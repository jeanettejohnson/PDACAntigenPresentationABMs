"""
assemble_initial_conditions.py
-------------------------------
Converts a QuPath annotation CSV + optional duct annotation GeoJSON into a
PhysiCell ICS cell positions file, saved to PhysiCell/config/ics/JHH_IMC/.

Steps:
  1. Select annotation coordinates CSV (from QuPath export)
  2. Optionally select a duct annotation GeoJSON (exported from QuPath)
  3. Enter image dimensions (default 1100x1100) and hex radius for duct grid
  4. Type mapping applied automatically
  5. Hex grid generated from GeoJSON, SAME coordinate transform applied to both
  6. Merged CSV saved with same stem as annotation file

IMPORTANT: both cell and duct coordinates are transformed using the same image
dimensions, guaranteeing they end up in the same coordinate space.
"""

import tkinter as tk
from tkinter import filedialog, simpledialog, messagebox
from pathlib import Path
import json
import math
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, shape

# ── Type mapping ──────────────────────────────────────────────────────────────
# Edit this dict to add new QuPath label -> PhysiCell type mappings
TYPE_MAP = {
    'CAF: Other':                       'CAF',
    'CAF: HLA-DR':                      'apCAF',
    'ductal: Other':                    'epithelial_normal',
    'ductal: HLA-DR':                   'epithelial_normal',
    'tumor_epithelial: Other':          'epithelial_tumor_class1',
    'tumor_epithelial: HLA-DR':         'epithelial_tumor_class1_class2',
    'tumor_mesenchymal: Other':         'mesenchymal_tumor_class1',
    'tumor_mesenchymal: HLA-DR':        'mesenchymal_tumor_class1_class2',
    'CD4 T cell: Other':                'CD4_Tcell',
    'CD4 T cell: HLA-DR':              'CD4_Tcell',
    'CD4 T cell: FOXP3':               'Treg',
    'CD4 T cell: FOXP3: HLA-DR':       'Treg',
    'CD8 T cell: Other':               'CD8_Tcell',
    'CD8 T cell: HLA-DR':              'CD8_Tcell',
    'Myeloid: Other':                  'macrophage',
    'Myeloid: HLA-DR':                 'macrophage',
    'B cell: Other':                   'B cell',
    'B cell: HLA-DR':                  'B cell',
    'CD57: Other':                     'CD8_Tcell',
    'CD57: HLA-DR':                    'CD8_Tcell',
    'Other: HLA-DR':                   'other_tissue',
    'Other':                           'other_tissue',
}

IMMOVABLE = {'duct_filler', 'other_tissue'}
HEX_RADIUS_DEFAULT = 8.0
OUTDIR = Path(__file__).parent / 'PhysiCell/config/ics/JHH_IMC'

# ── Hex grid helpers (from hex_mask_workflow.py) ──────────────────────────────
def make_hexagon(cx, cy, radius):
    pts = [(cx + radius * math.cos(math.radians(60*k+30)),
            cy + radius * math.sin(math.radians(60*k+30))) for k in range(6)]
    return Polygon(pts)

def build_hex_grid(bounds, radius):
    minx, miny, maxx, maxy = bounds
    dx = math.sqrt(3) * radius
    dy = 1.5 * radius
    hexes, row = [], 0
    y = miny - 2 * radius
    while y <= maxy + 2 * radius:
        x_offset = dx / 2 if row % 2 else 0
        x = minx - 2 * dx + x_offset
        while x <= maxx + 2 * dx:
            hexes.append(make_hexagon(x, y, radius))
            x += dx
        y += dy
        row += 1
    hex_gdf = gpd.GeoDataFrame({'geometry': hexes}, geometry='geometry')
    hex_gdf['hex_id'] = np.arange(len(hex_gdf))
    return hex_gdf

def hex_centers_from_geojson(geojson_path, radius):
    with open(geojson_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    features = data.get('features', [])
    records = [{'geometry': shape(feat['geometry'])} for feat in features]
    ann_gdf = gpd.GeoDataFrame(records, geometry='geometry')
    hex_gdf = build_hex_grid(ann_gdf.total_bounds, radius=radius)
    joined = gpd.sjoin(hex_gdf, ann_gdf, how='inner', predicate='intersects')
    masked = joined.drop_duplicates(subset=['hex_id'])
    centers = masked.geometry.centroid
    return pd.DataFrame({'x': centers.x, 'y': centers.y, 'z': 0.0, 'type': 'duct_filler'})

# ── GUI setup ─────────────────────────────────────────────────────────────────
root = tk.Tk()
root.withdraw()

# ── Step 1: Select QuPath annotation TXT file ─────────────────────────────────
print("Step 1: Select QuPath annotation TXT file...")
ann_path = filedialog.askopenfilename(
    title="Select QuPath annotation TXT file",
    filetypes=[("Text files", "*.txt"), ("CSV files", "*.csv"), ("All files", "*.*")]
)
if not ann_path:
    raise SystemExit("No annotation file selected.")
ann_path = Path(ann_path)
print(f"  Selected: {ann_path.name}")

# ── Step 2: Optionally select duct GeoJSON ────────────────────────────────────
print("Step 2: Select duct annotation GeoJSON from QuPath (cancel to skip)...")
geojson_path = filedialog.askopenfilename(
    title="Select duct annotation GeoJSON (cancel to skip)",
    filetypes=[("GeoJSON files", "*.geojson"), ("All files", "*.*")]
)
if geojson_path:
    geojson_path = Path(geojson_path)
    print(f"  Selected: {geojson_path.name}")
else:
    geojson_path = None
    print("  No GeoJSON selected — skipping duct_filler.")

# ── Step 3: Image dimensions and hex radius ───────────────────────────────────
dims = simpledialog.askstring(
    "Image dimensions",
    "Enter image dimensions as width,height in pixels\n(default: 1100,1100)",
    initialvalue="1100,1100"
)
if not dims:
    raise SystemExit("No dimensions entered.")
try:
    img_w, img_h = [float(v.strip()) for v in dims.split(',')]
except ValueError:
    raise SystemExit(f"Could not parse dimensions: {dims}")
print(f"  Image: {img_w:.0f} x {img_h:.0f} px")

if geojson_path is not None:
    hex_r = simpledialog.askfloat(
        "Hex radius",
        "Enter hex radius in pixels (default: 8.0)",
        initialvalue=HEX_RADIUS_DEFAULT
    )
    if hex_r is None:
        hex_r = HEX_RADIUS_DEFAULT
    print(f"  Hex radius: {hex_r}")

# ── Step 4: Load annotation file (TXT or CSV) and map types ──────────────────
if ann_path.suffix.lower() == '.txt':
    raw = pd.read_table(ann_path)
    cells = raw[['Classification', 'Centroid X µm', 'Centroid Y µm']].copy()
    cells.columns = ['type', 'x', 'y']
else:
    cells = pd.read_csv(ann_path)
print(f"\nLoaded {len(cells)} cells from annotation file.")

unmapped = set(cells['type'].unique()) - set(TYPE_MAP.keys())
if unmapped:
    print(f"\n⚠️  Unmapped cell types found: {unmapped}")
    response = messagebox.askyesno(
        "Unmapped types",
        f"The following types have no mapping and will become 'UNKNOWNTYPE':\n\n"
        + '\n'.join(sorted(unmapped))
        + "\n\nContinue?"
    )
    if not response:
        raise SystemExit("Aborted by user.")
    for t in unmapped:
        TYPE_MAP[t] = 'UNKNOWNTYPE'

cells['type'] = cells['type'].map(TYPE_MAP)

# ── Step 5: Coordinate transform on cells ─────────────────────────────────────
cells['x'] = cells['x'] - img_w / 2
cells['y'] = img_h / 2 - cells['y']
cells['z'] = 0.0
cells['is_movable'] = cells['type'].apply(lambda t: 0 if t in IMMOVABLE else 1)
cells = cells[['x', 'y', 'z', 'type', 'is_movable']]

# ── Step 6: Generate and transform duct_filler from GeoJSON ──────────────────
if geojson_path is not None:
    print(f"\nGenerating duct_filler hex grid from {geojson_path.name}...")
    ducts = hex_centers_from_geojson(geojson_path, radius=hex_r)
    print(f"  Generated {len(ducts)} hex centers.")
    ducts['x'] = ducts['x'] - img_w / 2
    ducts['y'] = img_h / 2 - ducts['y']
    ducts['is_movable'] = 0
    ducts = ducts[['x', 'y', 'z', 'type', 'is_movable']]
    combined = pd.concat([cells, ducts], ignore_index=True)
    print(f"  Transformed and appended {len(ducts)} duct_filler cells.")
else:
    combined = cells

# ── Step 7: Report ────────────────────────────────────────────────────────────
print("\n── Cell type counts ──────────────────────────────────────────")
print(combined['type'].value_counts().to_string())
print(f"\nTotal cells : {len(combined)}")
print(f"x range     : [{combined.x.min():.1f}, {combined.x.max():.1f}]")
print(f"y range     : [{combined.y.min():.1f}, {combined.y.max():.1f}]")

# ── Step 8: Save ──────────────────────────────────────────────────────────────
OUTDIR.mkdir(parents=True, exist_ok=True)
out_stem = ann_path.stem.replace('_for_physicell', '').replace('_summary', '')
out_path = OUTDIR / f"{out_stem}.csv"
combined.to_csv(out_path, index=False)
print(f"\n✓ Saved to: {out_path}")

