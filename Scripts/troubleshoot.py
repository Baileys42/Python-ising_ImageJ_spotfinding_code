"""
troubleshoot.py
Diagnoses TIFF loading issues — paste the path to one of your TIFF stacks below.
"""

from pathlib import Path
import numpy as np
import tifffile

TIFF_PATH = Path(r"C:\path\to\one\of\your\tiff\files.tif")

# ── Method 1: tifffile.imread (simplest) ─────────────────────────────────────
print("=" * 55)
print("Method 1: tifffile.imread")
arr1 = tifffile.imread(str(TIFF_PATH))
print(f"  shape : {arr1.shape}")
print(f"  dtype : {arr1.dtype}")
print(f"  ndim  : {arr1.ndim}")

# ── Method 2: TiffFile.asarray ────────────────────────────────────────────────
print("\nMethod 2: TiffFile.asarray")
with tifffile.TiffFile(str(TIFF_PATH)) as tif:
    arr2 = tif.asarray()
    print(f"  shape       : {arr2.shape}")
    print(f"  dtype       : {arr2.dtype}")
    print(f"  n_pages     : {len(tif.pages)}")
    print(f"  n_series    : {len(tif.series)}")
    if tif.series:
        s = tif.series[0]
        print(f"  series shape: {s.shape}")
        print(f"  series axes : {s.axes}")
    print(f"  is_imagej   : {tif.is_imagej}")
    print(f"  is_ome      : {tif.is_ome}")

# ── Method 3: read page by page ───────────────────────────────────────────────
print("\nMethod 3: page-by-page")
with tifffile.TiffFile(str(TIFF_PATH)) as tif:
    pages = [p.asarray() for p in tif.pages]
    arr3 = np.stack(pages)
    print(f"  n_pages read : {len(pages)}")
    print(f"  stacked shape: {arr3.shape}")

# ── After squeeze ─────────────────────────────────────────────────────────────
print("\nAfter np.squeeze on Method 2:")
squeezed = np.squeeze(arr2)
print(f"  shape: {squeezed.shape}")

print("\n" + "=" * 55)
print("CONCLUSION:")
print(f"  Use page-by-page → {arr3.shape[0]} frames, {arr3.shape[1]}×{arr3.shape[2]} px")
if arr3.shape[0] != arr2.shape[0]:
    print("\n  *** OME metadata is broken — asarray() collapses frames.")
    print("  *** Always use page-by-page loading for these files.")
