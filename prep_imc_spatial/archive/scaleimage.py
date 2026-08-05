import numpy as np
import pandas as pd
import tifffile
from PIL import Image
import tkinter as tk
from tkinter import filedialog

# ── Prompt user to select TIFF file ──────────────────────────────────────────
root = tk.Tk()
root.withdraw()
TIFF_PATH = filedialog.askopenfilename(
    title="Select OME-TIFF image",
    filetypes=[("TIFF files", "*.tiff *.tif"), ("All files", "*.*")]
)
if not TIFF_PATH:
    raise SystemExit("No file selected.")
print(f"Selected: {TIFF_PATH}")

# ── Prompt user to enter output CSV name ─────────────────────────────────────
csv_name = input("Enter a name for the output CSV (without .csv extension): ").strip()
if not csv_name:
    raise SystemExit("No filename entered.")
OUTPUT_DIR  = "/Users/jeanette.johnson/PDACAntigenPresentationABMs/PhysiCell/config/ics/substrate"
OUTPUT_PATH = f"{OUTPUT_DIR}/{csv_name}.csv"

# ── User-defined parameters ───────────────────────────────────────────────────
PERCENTILE  = 95    # clip ECM values above this percentile of nonzero values
DX          = 20.0  # voxel size in µm (set in PhysiCell XML)
DY          = 20.0  # voxel size in µm (set in PhysiCell XML)
# ─────────────────────────────────────────────────────────────────────────────

# ── Step 2: Read image and extract pixel size from metadata ──────────────────
with tifffile.TiffFile(TIFF_PATH) as tif:
    page = tif.pages[0]
    orig_w = page.tags["ImageWidth"].value
    orig_h = page.tags["ImageLength"].value
    xres_num, xres_den = page.tags["XResolution"].value
    yres_num, yres_den = page.tags["YResolution"].value
    res_unit = page.tags["ResolutionUnit"].value  # 2 = inch, 3 = cm
    pix_per_unit_x = xres_num / xres_den
    pix_per_unit_y = yres_num / yres_den
    if res_unit == 2:    # inches -> µm
        um_per_pix_x = 25400 / pix_per_unit_x
        um_per_pix_y = 25400 / pix_per_unit_y
    elif res_unit == 3:  # cm -> µm
        um_per_pix_x = 10000 / pix_per_unit_x
        um_per_pix_y = 10000 / pix_per_unit_y
    else:
        raise ValueError(f"Unknown ResolutionUnit: {res_unit}")

phys_w = orig_w * um_per_pix_x
phys_h = orig_h * um_per_pix_y
resize_w = round(phys_w / DX)
resize_h = round(phys_h / DY)
x_min, x_max = -phys_w / 2, phys_w / 2
y_min, y_max = -phys_h / 2, phys_h / 2

print(f"Original image: {orig_w} x {orig_h} pixels")
print(f"Pixel size: {um_per_pix_x:.4f} x {um_per_pix_y:.4f} µm/pixel")
print(f"Physical size: {phys_w:.1f} x {phys_h:.1f} µm")
print(f"\n── Copy these into PhysiCell XML <domain> ───────────────")
print(f"  <x_min>{x_min}</x_min>")
print(f"  <x_max>{x_max}</x_max>")
print(f"  <y_min>{y_min}</y_min>")
print(f"  <y_max>{y_max}</y_max>")
print(f"  <dx>{DX}</dx>")
print(f"  <dy>{DY}</dy>")
print(f"  <x_nodes>{resize_w}</x_nodes>")
print(f"  <y_nodes>{resize_h}</y_nodes>")
print(f"────────────────────────────────────────────────────────\n")

# ── Step 3: Resize image to voxel grid dimensions ────────────────────────────
print(f"Resizing image to: {resize_w} x {resize_h} voxels")

# ── Step 4: Read, resize, and transpose image (matching MATLAB pipeline) ─────
I = Image.open(TIFF_PATH)
B = np.array(I.resize((resize_w, resize_h), Image.BILINEAR))

# ── Step 5: Build coordinate grid and flatten (column-major, matching MATLAB) ─
rows, cols = B.shape
X, Y = np.meshgrid(np.arange(1, cols + 1), np.arange(1, rows + 1))
x_col = X.flatten(order='F')
y_col = Y.flatten(order='F')
z_col = B.flatten(order='F').astype(float)
image = pd.DataFrame({'x': x_col, 'y': y_col, 'ecm': z_col})

# ── Step 6: Scale ECM values to [0, 10] using nonzero 95th percentile ────────
nonzero_vals = image.loc[image["ecm"] > 0, "ecm"]
upper = np.percentile(nonzero_vals, PERCENTILE)
ecm_clipped = np.clip(image["ecm"].values, 0, upper)
image["ecm"] = (ecm_clipped / upper) * 10
print(f"ECM scaled: clipped at {PERCENTILE}th percentile ({upper:.2f}), range now [0, 10]")

# ── Step 7: Map pixel coordinates to PhysiCell voxel centers ─────────────────
# Map pixel index (1-indexed) directly to voxel center: x_min + (i - 0.5) * dx
# This guarantees spacing is exactly dx/dy — no voxel collisions after rotation.
spatial_base_coords = image[["x", "y"]].copy()
spatial_base_coords["x"] = x_min + (image["x"] - 0.5) * DX
spatial_base_coords["y"] = y_min + (image["y"] - 0.5) * DY
print(f"Coordinates mapped to voxel centers: x=[{spatial_base_coords['x'].min():.1f}, {spatial_base_coords['x'].max():.1f}], y=[{spatial_base_coords['y'].min():.1f}, {spatial_base_coords['y'].max():.1f}]")
spatial_base_coords["z"] = 0.0
spatial_base_coords["ecm"] = image["ecm"]

# ── Step 8: Flip around x axis (negate y) ────────────────────────────────────
spatial_base_coords["y"] = -spatial_base_coords["y"]

# ── Step 9: Save ──────────────────────────────────────────────────────────────
spatial_base_coords[["x", "y", "z", "ecm"]].to_csv(OUTPUT_PATH, index=False)
print(f"Saved to: {OUTPUT_PATH}")
