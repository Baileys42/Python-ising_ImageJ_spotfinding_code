# traj_macro.py

Python port of a Fiji/ImageJ single-molecule FRET analysis pipeline. Processes batches of client (488 nm) and chaperone (647 nm) TIRF microscopy image pairs through background correction, spot finding, trajectory extraction, and colocalisation analysis.

---

## Scripts

| Script | Purpose |
|---|---|
| `traj_macro.py` | Main pipeline — background correction, spot finding, trajectories, colocalisation |
| `trajectory_extractor.py` | Standalone trajectory extraction from an existing Results folder |
| `movie_vis.py` | Quick TIFF stack viewer (Jupyter) |
| `troubleshoot.py` | Diagnose OME-TIFF loading issues |

---

## Requirements

```
pip install tifffile numpy scipy pandas matplotlib tqdm
```

---

## Folder Structure

**Input:**

```
Exp/
  Client/
    client_1.tif
    client_2.tif
  Chap/
    chap_1.tif
    chap_2.tif
```

Images are paired by the trailing integer in their filename (`client_1` ↔ `chap_1`).

**Output (auto-generated):**

```
Exp/
  Results/
    client_1/
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
```

---

## Usage

```
python traj_macro.py
```

Two folder dialogs open — select your **Client** folder then your **Chap** folder. The script runs in two phases:

**Phase 1 — Spot finding (all pairs)**
All interactive peak review windows are shown back-to-back so you can review everything in one go before any computation starts.

**Phase 2 — Trajectory extraction**
Prompted with `Show live trajectory window? (y/n)`. If yes, a single window shows each molecule's trace in white with the running mean in blue as processing happens in real time.

---

## CONFIG

All tunable parameters are at the top of `traj_macro.py`:

| Parameter | Default | Description |
|---|---|---|
| `CLIENT_ELECTRONIC_OFFSET` | `700` | Camera electronic offset, client channel |
| `CHAP_ELECTRONIC_OFFSET` | `500` | Camera electronic offset, chap channel |
| `BACKGROUND_BLUR_SIGMA` | `40` | Gaussian sigma (px) for beam profile estimation |
| `CLIENT_MAX_FRAMES` | `20` | Frames used for client max-intensity projection |
| `PEAK_THRESHOLD_SIGMA` | `6.0` | Peak threshold: mean + N × σ of filtered image |
| `PEAK_MIN_DISTANCE` | `8` | Minimum distance between peaks (px) |
| `TRAJ_INNER_RADIUS` | `2` | Inner box radius for intensity measurement (5×5 px) |
| `TRAJ_OUTER_RADIUS` | `4` | Outer box radius for background estimation (9×9 px) |
| `COLOC_MAX_DISTANCE` | `3.0` | Max pixel distance to call two peaks colocalised |

---

## Pipeline Detail

### Background correction

Per-pixel minimum projection across all frames is Gaussian blurred (σ = 40 px) to estimate the beam profile. Molecules blink/bleach so they are suppressed in the minimum, leaving a clean illumination profile. Correction formula (per `BeamProfileCorrection.java`):

```
output = |(pixel − offset) / ((bg − offset) / (max_bg − offset))|
```

### Spot finding

A discoidal averaging filter (inner disk mean minus outer annulus mean) enhances point-like fluorescent spots. Peaks are found by iteratively picking the global maximum and blanking a minimum-distance radius around it. An interactive matplotlib window lets you adjust threshold σ and minimum distance before accepting.

- **Client** — max intensity projection of first 20 frames
- **Chap** — average intensity projection of all frames

### Trajectory extraction

For each detected peak and each frame:

```
intensity = mean(5×5 box) − mean(9×9 box)
```

Output CSV: rows = frames, columns = `Mean_0`, `Mean_1`, … (one column per peak).

### Colocalisation

Greedy nearest-neighbour matching between client and chap peak tables. Closest pairs are matched first. A match is called if two peaks are within `COLOC_MAX_DISTANCE` pixels. Output adds `X2`, `Y2`, `distance` columns to the client results table (`-1` = no match).

---

## Notes

- Input files must be OME-TIFF or standard multi-page TIFF stacks. The loader reads page-by-page to work around a known `tifffile` bug where OME metadata incorrectly reports frame count.
- Run `troubleshoot.py` on any TIFF that behaves unexpectedly to diagnose loading issues.
