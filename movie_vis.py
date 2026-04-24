# %% Imports
%matplotlib widget
import ipywidgets as widgets
from IPython.display import display
from pathlib import Path
import numpy as np
import tifffile
import matplotlib.pyplot as plt

# %% Config
BASE_FPS      = 10      # assumed acquisition frame rate
DEFAULT_SPEED = 5.0     # playback speed multiplier
MAX_MOVIES    = 12
NROWS, NCOLS  = 6, 2

# %% Load movies — set your folder path here
folder = Path(r"C:\path\to\your\movies")

paths = sorted(list(folder.glob('*.tif')) + list(folder.glob('*.tiff')))[:MAX_MOVIES]
print(f"Found {len(paths)} movie(s) in {folder}")

movies = []
for p in paths:
    print(f"  Loading {p.name}…")
    with tifffile.TiffFile(str(p)) as tif:
        arr = tif.asarray()
    print(f"    raw shape: {arr.shape}, dtype: {arr.dtype}")

    # Squeeze any length-1 axes except the last two (H, W)
    arr = np.squeeze(arr)

    # Ensure (T, H, W) — if 2D it's a single frame, if >3D take first channel
    if arr.ndim == 2:
        arr = arr[np.newaxis]
    elif arr.ndim == 4:
        arr = arr[:, 0]     # drop channel axis, keep frames

    print(f"    interpreted as: {arr.shape[0]} frames, {arr.shape[1]}×{arr.shape[2]} px")

    lo, hi = np.percentile(arr, [0.5, 99.5])
    arr_f = np.clip((arr.astype(np.float32) - lo) / max(float(hi - lo), 1.0), 0.0, 1.0)
    movies.append({'name': p.stem, 'frames': arr_f, 'n': arr_f.shape[0]})

max_frames = max(m['n'] for m in movies)
print(f"Loaded. Max frames: {max_frames}")

# %% Display
n = len(movies)

fig, axes = plt.subplots(NROWS, NCOLS,
                          figsize=(NCOLS * 5, NROWS * 2.4),
                          facecolor='#1a1a1a')
fig.subplots_adjust(left=0.01, right=0.99, top=0.97,
                    bottom=0.01, hspace=0.12, wspace=0.04)

ims = []
for i, ax in enumerate(axes.flat):
    ax.set_facecolor('black')
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    if i < n:
        ax.set_title(movies[i]['name'], fontsize=7, color='white', pad=2)
        im = ax.imshow(movies[i]['frames'][0], cmap='gray',
                       vmin=0, vmax=1, aspect='equal', interpolation='nearest')
        ims.append(im)
    else:
        ax.set_visible(False)

frame_text = fig.text(0.5, 0.002, f"Frame 1 / {max_frames}",
                       ha='center', color='white', fontsize=8)
plt.show()

# ── Widgets ───────────────────────────────────────────────────────────────────
play = widgets.Play(
    value=0, min=0, max=max_frames - 1, step=1,
    interval=int(1000 / (BASE_FPS * DEFAULT_SPEED)),
    description='', layout=widgets.Layout(width='60px'),
)
frame_slider = widgets.IntSlider(
    min=0, max=max_frames - 1, value=0,
    description='Frame',
    layout=widgets.Layout(width='500px'),
    style={'description_width': '45px'},
)
speed_slider = widgets.FloatSlider(
    min=0.5, max=20.0, step=0.5, value=DEFAULT_SPEED,
    description='Speed ×',
    layout=widgets.Layout(width='400px'),
    style={'description_width': '55px'},
)

widgets.jslink((play, 'value'), (frame_slider, 'value'))

def on_frame_change(change):
    f = change['new']
    for i, movie in enumerate(movies):
        ims[i].set_data(movie['frames'][f % movie['n']])
    frame_text.set_text(f"Frame {f + 1} / {max_frames}")
    fig.canvas.draw_idle()

def on_speed_change(change):
    play.interval = int(1000 / (BASE_FPS * change['new']))

play.observe(on_frame_change, names='value')
speed_slider.observe(on_speed_change, names='value')

display(widgets.VBox([
    widgets.HBox([play, frame_slider]),
    speed_slider,
]))
