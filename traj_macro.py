#!/usr/bin/env python3
"""
traj_macro.py — Batch single-molecule FRET pipeline.

For each matched client/chap image pair:
  1. Generate background from last frame (Gaussian blur, sigma=40)
  2. Beam-profile correction
  3. Interactive peak finding (adjust threshold, accept)
  4. Trajectory extraction (background-subtracted intensity per peak per frame)
  5. Client–chap colocalisation

Expected folder layout:
  Exp/
    Client/   client_1.tif  client_2.tif  ...
    Chap/     chap_1.tif    chap_2.tif    ...

Output created at Exp/Results/client_1/ etc.:
  client_1_background_corrected.tif
  chap_1_background_corrected.tif
  client_1_peaks_preview.png
  chap_1_peaks_preview.png
  client_1_results.csv
  chap_1_results.csv
  client_1_colocalisation.csv
  Client_trajectories/
    client_1_trajectories.csv
  Chap_trajectories/
    chap_1_trajectories.csv

Dependencies:
    pip install tifffile numpy scipy pandas matplotlib
"""

import re
import sys
import warnings
from pathlib import Path
from tqdm import tqdm

import numpy as np
import pandas as pd
import tifffile
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, TextBox
from scipy.ndimage import convolve, gaussian_filter
from tkinter import filedialog, Tk

warnings.filterwarnings('ignore')


# ─────────────────────────────────────────────────────────────────────────────
# CONFIG — edit to match your acquisition settings
# ─────────────────────────────────────────────────────────────────────────────

CLIENT_ELECTRONIC_OFFSET = 700    # camera electronic offset, client channel
CHAP_ELECTRONIC_OFFSET   = 500    # camera electronic offset, chap channel
BACKGROUND_BLUR_SIGMA    = 40     # Gaussian sigma (pixels) for background from last frame
CLIENT_MAX_FRAMES        = 20     # frames used for client max-intensity projection
PEAK_INNER_RADIUS        = 1      # discoidal filter inner radius (pixels)
PEAK_OUTER_RADIUS        = 3      # discoidal filter outer radius (pixels)
PEAK_THRESHOLD_SIGMA     = 6.0    # threshold = mean + N * std of filtered image
PEAK_MIN_DISTANCE        = 8      # minimum separation between peaks (pixels)
TRAJ_INNER_RADIUS        = 2      # trajectory measurement inner box radius → 5×5 px
TRAJ_OUTER_RADIUS        = 4      # trajectory measurement outer box radius → 9×9 px
COLOC_MAX_DISTANCE       = 3.0    # max pixel distance to call two peaks colocalised


# ─────────────────────────────────────────────────────────────────────────────
# Background estimation
# ─────────────────────────────────────────────────────────────────────────────

def make_background(stack: np.ndarray) -> np.ndarray:
    """
    Gaussian blur of the per-pixel minimum across all frames.
    Molecules blink/bleach so they're suppressed in the minimum,
    leaving a cleaner estimate of the beam profile than a single frame.
    """
    return gaussian_filter(stack.min(axis=0).astype(float), sigma=BACKGROUND_BLUR_SIGMA)


# ─────────────────────────────────────────────────────────────────────────────
# Beam profile correction  (BeamProfileCorrection.java)
# Formula: |(pixel − offset) / ((bg − offset) / (max_bg − offset))|
# ─────────────────────────────────────────────────────────────────────────────

def beam_profile_correction(
    stack: np.ndarray,
    background: np.ndarray,
    electronic_offset: float,
) -> np.ndarray:
    bg = background.astype(float)
    valid = bg[bg < 65535]
    max_bg = valid.max() if valid.size > 0 else bg.max()

    bg_norm = (bg - electronic_offset) / (max_bg - electronic_offset)
    bg_norm = np.where(np.abs(bg_norm) < 1e-10, 1e-10, bg_norm)

    s = stack.astype(float)
    if s.ndim == 3:
        bg_norm = bg_norm[np.newaxis]   # broadcast over frames

    return np.abs((s - electronic_offset) / bg_norm)


# ─────────────────────────────────────────────────────────────────────────────
# Peak finder  (PeakFinder.java + DiscoidalAveragingFilter.java)
# ─────────────────────────────────────────────────────────────────────────────

def _disk_kernel(radius: int) -> np.ndarray:
    y, x = np.ogrid[-radius: radius + 1, -radius: radius + 1]
    mask = x ** 2 + y ** 2 <= radius ** 2
    return mask.astype(float) / mask.sum()


def discoidal_filter(image: np.ndarray) -> np.ndarray:
    """
    Inner-disk mean minus outer-annulus mean.
    Enhances point-like fluorescent spots against background.
    """
    inner_k = _disk_kernel(PEAK_INNER_RADIUS)
    r = PEAK_OUTER_RADIUS
    y, x = np.ogrid[-r: r + 1, -r: r + 1]
    annulus = (x**2 + y**2 <= r**2) & ~(x**2 + y**2 <= PEAK_INNER_RADIUS**2)
    annulus_k = annulus.astype(float) / annulus.sum()
    f = image.astype(float)
    return convolve(f, inner_k, mode='reflect') - convolve(f, annulus_k, mode='reflect')


def find_peaks(
    image: np.ndarray,
    threshold_sigma: float = PEAK_THRESHOLD_SIGMA,
    threshold_value: float = 0.0,
    min_distance: int = PEAK_MIN_DISTANCE,
) -> list:
    """
    Returns list of (x, y) peak coordinates.
    Applies discoidal filter, thresholds at mean + N*σ, then iteratively
    picks the global maximum and blanks a radius around it.
    """
    filtered = discoidal_filter(image)
    t = threshold_value if threshold_value > 0 else (filtered.mean() + threshold_sigma * filtered.std())

    working = filtered.copy()
    peaks = []

    while working.max() >= t:
        y, x = np.unravel_index(working.argmax(), working.shape)
        peaks.append((x, y))
        rr = np.arange(max(0, y - min_distance), min(working.shape[0], y + min_distance + 1))
        cc = np.arange(max(0, x - min_distance), min(working.shape[1], x + min_distance + 1))
        gy, gx = np.meshgrid(rr, cc, indexing='ij')
        mask = (gy - y) ** 2 + (gx - x) ** 2 <= min_distance ** 2
        working[gy[mask], gx[mask]] = -np.inf

    return peaks


# ─────────────────────────────────────────────────────────────────────────────
# Interactive peak review
# ─────────────────────────────────────────────────────────────────────────────

def _display(image: np.ndarray, saturated: float = 0.35) -> np.ndarray:
    lo, hi = np.percentile(image, [saturated, 100 - saturated])
    return np.clip((image - lo) / max(hi - lo, 1e-6), 0, 1)


def review_peaks(image: np.ndarray, initial_peaks: list, title: str = "") -> list:
    """
    Shows the Z-projection with detected peaks overlaid.
    Adjust threshold σ and minimum distance, click Re-run to update,
    click Accept when happy.
    """
    state = {'peaks': list(initial_peaks)}

    fig, ax = plt.subplots(figsize=(10, 8))
    plt.subplots_adjust(bottom=0.28)
    ax.imshow(_display(image), cmap='gray')
    ax.set_title(f"{title}  —  {len(state['peaks'])} peaks", fontsize=11)

    offsets = np.array([[p[0], p[1]] for p in state['peaks']]) if state['peaks'] else np.empty((0, 2))
    pts = ax.scatter(*(offsets.T if offsets.size else ([], [])),
                     s=60, facecolors='none', edgecolors='cyan', linewidths=1)

    ax_sigma  = plt.axes([0.12, 0.16, 0.18, 0.05])
    ax_dist   = plt.axes([0.38, 0.16, 0.18, 0.05])
    ax_rerun  = plt.axes([0.12, 0.07, 0.20, 0.07])
    ax_accept = plt.axes([0.65, 0.07, 0.22, 0.07])

    tb_sigma  = TextBox(ax_sigma, 'Threshold σ', initial=str(PEAK_THRESHOLD_SIGMA))
    tb_dist   = TextBox(ax_dist,  'Min dist px',  initial=str(PEAK_MIN_DISTANCE))
    btn_rerun  = Button(ax_rerun,  'Re-run')
    btn_accept = Button(ax_accept, 'Accept', color='lightgreen')

    def rerun(_):
        try:
            sig  = float(tb_sigma.text)
            dist = int(float(tb_dist.text))
        except ValueError:
            return
        new_peaks = find_peaks(image, threshold_sigma=sig, min_distance=dist)
        state['peaks'] = new_peaks
        new_offsets = np.array([[p[0], p[1]] for p in new_peaks]) if new_peaks else np.empty((0, 2))
        pts.set_offsets(new_offsets)
        ax.set_title(f"{title}  —  {len(new_peaks)} peaks", fontsize=11)
        fig.canvas.draw_idle()

    btn_rerun.on_clicked(rerun)
    btn_accept.on_clicked(lambda _: plt.close(fig))
    plt.show()
    return state['peaks']


# ─────────────────────────────────────────────────────────────────────────────
# Trajectory measurement  (peak_intensity.txt)
# mean(inner box) − mean(outer box) per peak per frame
# ─────────────────────────────────────────────────────────────────────────────

def measure_trajectories(stack: np.ndarray, peaks: list, on_peak=None) -> pd.DataFrame:
    """
    Returns DataFrame: rows = frames, columns = Mean_0, Mean_1, …
    on_peak: optional callback(i, traj_array) called after each peak is computed.
    """
    n_frames, H, W = stack.shape
    ri, ro = TRAJ_INNER_RADIUS, TRAJ_OUTER_RADIUS
    data = np.zeros((n_frames, len(peaks)))

    for i, (cx, cy) in enumerate(tqdm(peaks, desc="    Trajectories", unit="peak")):
        cx, cy = int(round(cx)), int(round(cy))
        traj = np.zeros(n_frames)
        for f in range(n_frames):
            img = stack[f].astype(float)
            yi = slice(max(0, cy - ri), min(H, cy + ri + 1))
            xi = slice(max(0, cx - ri), min(W, cx + ri + 1))
            yo = slice(max(0, cy - ro), min(H, cy + ro + 1))
            xo = slice(max(0, cx - ro), min(W, cx + ro + 1))
            traj[f] = img[yi, xi].mean() - img[yo, xo].mean()
        data[:, i] = traj
        if on_peak is not None:
            on_peak(i, traj)

    return pd.DataFrame(data, columns=[f'Mean_{i}' for i in range(len(peaks))])


# ─────────────────────────────────────────────────────────────────────────────
# Colocalisation  (count_colocalized_peaks_Andrew.js)
# Greedy nearest-neighbour, closest pairs matched first
# ─────────────────────────────────────────────────────────────────────────────

def colocalize(table1: pd.DataFrame, table2: pd.DataFrame) -> pd.DataFrame:
    """
    Adds X2, Y2, distance columns to a copy of table1.
    distance = -1 means no match within COLOC_MAX_DISTANCE.
    """
    t1 = table1.copy().reset_index(drop=True)
    t2 = table2.reset_index(drop=True)

    candidates = []
    for i, r1 in t1.iterrows():
        for j, r2 in t2.iterrows():
            d = np.hypot(r1['X'] - r2['X'], r1['Y'] - r2['Y'])
            if d < COLOC_MAX_DISTANCE:
                candidates.append([i, j, d])

    candidates.sort(key=lambda c: c[2])   # closest first
    t1['X2'] = np.nan
    t1['Y2'] = np.nan
    t1['distance'] = -1.0

    while candidates:
        i, j, d = candidates.pop(0)
        if i == -1:
            continue
        t1.at[i, 'X2'] = t2.at[j, 'X']
        t1.at[i, 'Y2'] = t2.at[j, 'Y']
        t1.at[i, 'distance'] = d
        for c in candidates:
            if c[1] == j:
                c[0] = -1   # chap peak already matched; invalidate competing candidates

    return t1


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def load_stack(path: Path) -> np.ndarray:
    # Page-by-page avoids broken OME metadata parsing in tifffile
    with tifffile.TiffFile(str(path)) as tif:
        pages = tif.pages
        arr = np.stack([p.asarray() for p in tqdm(pages, desc=f"    Loading {path.name}", unit="frame")])
    print(f"    → {arr.shape[0]} frames, {arr.shape[1]}×{arr.shape[2]} px")
    return arr


def save_tiff(arr: np.ndarray, path: Path) -> None:
    if arr.dtype in (np.float32, np.float64):
        arr = np.clip(arr, 0, 65535).astype(np.uint16)
    tifffile.imwrite(str(path), arr)


def z_project(stack: np.ndarray, mode: str = 'max') -> np.ndarray:
    return stack.max(axis=0) if mode == 'max' else stack.mean(axis=0)


def peaks_to_df(peaks: list) -> pd.DataFrame:
    if not peaks:
        return pd.DataFrame(columns=['X', 'Y', 'Slice'])
    return pd.DataFrame({'X': [p[0] for p in peaks],
                         'Y': [p[1] for p in peaks],
                         'Slice': 1})


def save_peak_preview(image: np.ndarray, peaks: list, path: Path, title: str = "") -> None:
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(_display(image), cmap='gray')
    if peaks:
        xs, ys = zip(*peaks)
        ax.scatter(xs, ys, s=60, facecolors='none', edgecolors='cyan', linewidths=1)
    ax.set_title(f"{title}  ({len(peaks)} peaks)", fontsize=10)
    ax.axis('off')
    fig.savefig(str(path), dpi=150, bbox_inches='tight')
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Image pairing
# ─────────────────────────────────────────────────────────────────────────────

def pair_images(client_folder: Path, chap_folder: Path) -> list:
    """
    Match client_N.tif to chap_N.tif by the trailing integer in the filename.
    """
    def idx(p: Path) -> int:
        nums = re.findall(r'\d+', p.stem)
        return int(nums[-1]) if nums else -1

    clients = {idx(p): p for p in sorted(client_folder.glob('*.tif'))}
    chaps   = {idx(p): p for p in sorted(chap_folder.glob('*.tif'))}

    paired = []
    for n in sorted(clients):
        if n in chaps:
            paired.append((clients[n], chaps[n]))
        else:
            print(f"  Warning: no chap match for {clients[n].name}, skipping.")
    for n in chaps:
        if n not in clients:
            print(f"  Warning: no client match for {chaps[n].name}, skipping.")

    return paired


# ─────────────────────────────────────────────────────────────────────────────
# Process one pair
# ─────────────────────────────────────────────────────────────────────────────

def process_pair(client_path: Path, chap_path: Path, results_root: Path) -> None:
    name = client_path.stem
    out  = results_root / name
    (out / "Client_trajectories").mkdir(parents=True, exist_ok=True)
    (out / "Chap_trajectories").mkdir(parents=True, exist_ok=True)

    print(f"\n{'─'*60}")
    print(f"  {client_path.name}  ←→  {chap_path.name}")

    # ── Client ────────────────────────────────────────────────────────────
    print("  [Client] Loading and correcting…")
    client_stack = load_stack(client_path)
    client_bg    = make_background(client_stack)
    client_corr  = beam_profile_correction(client_stack, client_bg, CLIENT_ELECTRONIC_OFFSET)
    client_corr  = np.clip(client_corr, 0, 65535).astype(np.uint16)
    save_tiff(client_corr, out / f"{client_path.stem}_background_corrected.tif")
    del client_stack  # free memory — corrected version is saved to disk

    n_max        = min(CLIENT_MAX_FRAMES, client_corr.shape[0])
    client_proj  = z_project(client_corr[:n_max], 'max')
    del client_corr

    print("  [Client] Peak finding — review window opening…")
    client_peaks = find_peaks(client_proj)
    client_peaks = review_peaks(client_proj, client_peaks, title=f"Client  {name}")
    save_peak_preview(client_proj, client_peaks,
                      out / f"{client_path.stem}_peaks_preview.png", title=f"Client {name}")
    client_df    = peaks_to_df(client_peaks)
    client_df.to_csv(out / f"{client_path.stem}_results.csv", index=False)
    print(f"  [Client] {len(client_peaks)} peaks accepted.")

    # ── Chap ──────────────────────────────────────────────────────────────
    print("  [Chap] Loading and correcting…")
    chap_stack = load_stack(chap_path)
    chap_bg    = make_background(chap_stack)
    chap_corr  = beam_profile_correction(chap_stack, chap_bg, CHAP_ELECTRONIC_OFFSET)
    chap_corr  = np.clip(chap_corr, 0, 65535).astype(np.uint16)
    save_tiff(chap_corr, out / f"{chap_path.stem}_background_corrected.tif")
    del chap_stack
    del chap_corr

    chap_proj  = z_project(np.stack([p.asarray() for p in
                     tqdm(tifffile.TiffFile(str(out / f"{chap_path.stem}_background_corrected.tif")).pages,
                          desc="    Reloading chap for projection", unit="frame", leave=False)]),
                     'average')

    print("  [Chap] Peak finding — review window opening…")
    chap_peaks = find_peaks(chap_proj)
    chap_peaks = review_peaks(chap_proj, chap_peaks, title=f"Chap  {name}")
    save_peak_preview(chap_proj, chap_peaks,
                      out / f"{chap_path.stem}_peaks_preview.png", title=f"Chap {name}")
    chap_df    = peaks_to_df(chap_peaks)
    chap_df.to_csv(out / f"{chap_path.stem}_results.csv", index=False)
    print(f"  [Chap] {len(chap_peaks)} peaks accepted.")

    # ── Colocalisation ─────────────────────────────────────────────────────
    print("  Colocalisation…")
    coloc_df = colocalize(client_df, chap_df)
    coloc_df.to_csv(out / f"{client_path.stem}_colocalisation.csv", index=False)
    n_coloc = (coloc_df['distance'] >= 0).sum()
    print(f"  → {n_coloc} / {len(client_peaks)} client peaks colocalised with chap.")

    return out


def run_trajectories(out_dir: Path, live_ax=None) -> None:
    """Phase 2: trajectory extraction for one pair, using saved corrected TIFFs and CSVs."""

    corrected = {p.stem.replace('_background_corrected', ''): p
                 for p in out_dir.glob('*_background_corrected.tif')}
    csvs      = {p.stem.replace('_results', ''): p
                 for p in out_dir.glob('*_results.csv')}

    for channel, traj_folder in [('client', 'Client_trajectories'),
                                  ('chap',   'Chap_trajectories')]:
        key = next((k for k in corrected if channel in k), None)
        if key is None:
            print(f"  [{channel}] No corrected TIFF found — skipping.")
            continue

        print(f"\n  [{channel}] Loading corrected stack…")
        stack  = load_stack(corrected[key])
        coords = pd.read_csv(csvs[key])
        n_peaks = len(coords)
        n_frames = stack.shape[0]
        print(f"  [{channel}] {n_peaks} peaks × {n_frames} frames")

        # Clear the live window for this channel
        if live_ax is not None:
            live_ax.clear()
            _style_traj_ax(live_ax)
            live_ax.set_title(f"{out_dir.name}  —  {channel}  (0 / {n_peaks} peaks)",
                              color='white', fontsize=10)
            live_ax.figure.canvas.draw()
            live_ax.figure.canvas.flush_events()

        mean_line = None

        def on_peak(i, traj):
            nonlocal mean_line
            if live_ax is None:
                return
            live_ax.plot(traj, color='steelblue', alpha=0.2, linewidth=0.5)
            # update mean line
            if mean_line is not None:
                mean_line.remove()
            current_data = np.array([live_ax.lines[j].get_ydata()
                                     for j in range(len(live_ax.lines))])
            mean_line, = live_ax.plot(current_data.mean(axis=0),
                                      color='white', linewidth=1.5)
            live_ax.set_title(f"{out_dir.name}  —  {channel}  ({i+1} / {n_peaks} peaks)",
                              color='white', fontsize=10)
            live_ax.figure.canvas.draw()
            live_ax.figure.canvas.flush_events()

        traj_df = measure_trajectories(stack, list(zip(coords['X'], coords['Y'])),
                                       on_peak=on_peak if live_ax is not None else None)
        out_path = out_dir / traj_folder / f"{key}_trajectories.csv"
        traj_df.to_csv(out_path, index=False)
        print(f"  [{channel}] Saved → {out_path.name}")


def _style_traj_ax(ax) -> None:
    ax.set_facecolor('#1a1a1a')
    ax.figure.patch.set_facecolor('#1a1a1a')
    ax.tick_params(colors='white')
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_color('white')
    ax.set_xlabel('Frame', color='white')
    ax.set_ylabel('Intensity (a.u.)', color='white')


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    root = Tk()
    root.withdraw()

    print("Select Client image folder…")
    client_folder = Path(filedialog.askdirectory(title="Select Client folder  (e.g. Exp/Client)"))
    if not str(client_folder):
        print("Cancelled.")
        return

    print("Select Chap image folder…")
    chap_folder = Path(filedialog.askdirectory(title="Select Chap folder  (e.g. Exp/Chap)"))
    if not str(chap_folder):
        print("Cancelled.")
        return

    root.destroy()
    plt.switch_backend('TkAgg')

    results_root = client_folder.parent / "Results"
    results_root.mkdir(exist_ok=True)

    print(f"\nClient : {client_folder}")
    print(f"Chap   : {chap_folder}")
    print(f"Output : {results_root}")

    pairs = pair_images(client_folder, chap_folder)
    if not pairs:
        print("\nNo matched pairs found. Check filenames contain matching integers.")
        return

    print(f"\nFound {len(pairs)} pair(s):")
    for c, h in pairs:
        print(f"  {c.name}  ←→  {h.name}")

    # ── Phase 1: spot finding for ALL pairs ───────────────────────────────
    print(f"\n{'═'*60}")
    print("PHASE 1 — Spot finding (all pairs)")
    print(f"{'═'*60}")
    out_dirs = []
    for client_path, chap_path in pairs:
        out_dir = process_pair(client_path, chap_path, results_root)
        out_dirs.append(out_dir)

    # ── Phase 2: trajectory extraction ────────────────────────────────────
    print(f"\n{'═'*60}")
    print("PHASE 2 — Trajectory extraction")
    print(f"{'═'*60}")
    show_plots = input("\nShow live trajectory window? (y/n): ").strip().lower() == 'y'

    live_ax = None
    if show_plots:
        plt.ion()
        fig, live_ax = plt.subplots(figsize=(11, 4))
        _style_traj_ax(live_ax)
        fig.tight_layout()
        plt.show(block=False)

    for out_dir in tqdm(out_dirs, desc="Pairs", unit="pair"):
        print(f"\n  {out_dir.name}")
        run_trajectories(out_dir, live_ax=live_ax)

    if show_plots:
        plt.ioff()
        plt.show(block=True)   # keep window open after processing finishes

    print(f"\n{'═'*60}")
    print(f"Done. All results in:  {results_root}")


if __name__ == "__main__":
    main()
