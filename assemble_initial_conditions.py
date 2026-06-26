"""
assemble_initial_conditions.py
-------------------------------
Converts a QuPath detection TXT + optional duct GeoJSON into a PhysiCell ICS
CSV, saved to PhysiCell/config/ics/JHH_IMC/.

Run interactively (GUI file-picker) or import process_roi() for batch use.
"""

import re
import json
import math
import numpy as np
import pandas as pd
import geopandas as gpd
from pathlib import Path
from shapely.geometry import Polygon, shape

# ── Type mapping ──────────────────────────────────────────────────────────────
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

IMMOVABLE          = {'duct_filler', 'other_tissue'}
WINSOR_LO          = 5
WINSOR_HI          = 95
WINSOR_SKIP_TYPES  = {'duct_filler'}
WINSOR_MIN_CELLS   = 10
DUCT_FILLER_VOLUME = 1247.0   # half the original volume (~13 µm diameter)
HEX_RADIUS_DEFAULT = 4.0    # half the spacing → ~4× more cells in the same area
OUTDIR             = Path(__file__).parent / 'PhysiCell/config/ics/JHH_IMC'

# Induced types whose ICS volume should match their parent type's distribution.
PARENT_VOLUME_TYPE = {
    'Treg':                              'CD4_Tcell',
    'CD8_exhausted':                     'CD8_Tcell',
    'apCAF':                             'CAF',
    'epithelial_tumor_class1_class2':    'epithelial_tumor_class1',
    'epithelial_tumor_class2':           'epithelial_tumor_class1',
    'mesenchymal_tumor_class1_class2':   'mesenchymal_tumor_class1',
    'mesenchymal_tumor_class2':          'mesenchymal_tumor_class1',
}

# ── Hex grid helpers ──────────────────────────────────────────────────────────
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

def _lookup_type(c):
    """Try classification as-is, then reversed (handles QuPath classifier order changes)."""
    if c in TYPE_MAP:
        return TYPE_MAP[c]
    reversed_c = ': '.join(p.strip() for p in reversed(c.split(':')))
    return TYPE_MAP.get(reversed_c)

# ── Core processing ───────────────────────────────────────────────────────────
def process_roi(ann_path, geojson_path=None, img_w=1100.0, img_h=1100.0,
                hex_r=HEX_RADIUS_DEFAULT, batch=False):
    """
    Convert one ROI's detection TXT (+ optional duct GeoJSON) to a PhysiCell
    ICS CSV.

    Parameters
    ----------
    ann_path     : Path  – QuPath detection TXT file
    geojson_path : Path or None – duct annotation GeoJSON
    img_w, img_h : float – image dimensions in µm (used for coordinate centring)
    hex_r        : float – hex radius for duct grid (µm)
    batch        : bool  – if True, unmapped types become UNKNOWNTYPE automatically
                           instead of prompting

    Returns
    -------
    Path – the saved ICS CSV
    """
    ann_path = Path(ann_path)

    # Load detections
    if ann_path.suffix.lower() == '.txt':
        raw   = pd.read_table(ann_path)
        cells = raw[['Classification', 'Centroid X µm', 'Centroid Y µm',
                      'Max diameter µm']].copy()
        cells.columns = ['type', 'x', 'y', 'max_diam']
    else:
        cells = pd.read_csv(ann_path)
    print(f"  Loaded {len(cells)} cells from {ann_path.name}")

    # Map QuPath classifications → PhysiCell types
    unmapped = {c for c in cells['type'].unique() if _lookup_type(c) is None}
    if unmapped:
        print(f"  ⚠  Unmapped types: {sorted(unmapped)}")
        if not batch:
            from tkinter import messagebox
            response = messagebox.askyesno(
                "Unmapped types",
                f"The following types have no mapping and will become 'UNKNOWNTYPE':\n\n"
                + '\n'.join(sorted(unmapped)) + "\n\nContinue?"
            )
            if not response:
                raise SystemExit("Aborted by user.")
        for t in unmapped:
            TYPE_MAP[t] = 'UNKNOWNTYPE'

    cells['type'] = cells['type'].map(_lookup_type)

    # Coordinate transform: centre at (0,0), flip y
    cells['x'] = cells['x'] - img_w / 2
    cells['y'] = img_h / 2 - cells['y']
    cells['z'] = 0.0
    cells['is_movable'] = cells['type'].apply(lambda t: 0 if t in IMMOVABLE else 1)
    if 'max_diam' in cells.columns:
        cells['volume'] = (np.pi / 6.0) * cells['max_diam'] ** 3
    cells = cells[['x', 'y', 'z', 'type', 'is_movable', 'volume']]

    # Duct filler from GeoJSON
    if geojson_path is not None:
        geojson_path = Path(geojson_path)
        print(f"  Generating duct_filler hex grid from {geojson_path.name}...")
        ducts = hex_centers_from_geojson(geojson_path, radius=hex_r)
        ducts['x'] = ducts['x'] - img_w / 2
        ducts['y'] = img_h / 2 - ducts['y']
        ducts['is_movable'] = 0
        ducts['volume']     = DUCT_FILLER_VOLUME
        ducts = ducts[['x', 'y', 'z', 'type', 'is_movable', 'volume']]
        combined = pd.concat([cells, ducts], ignore_index=True)
        print(f"  Appended {len(ducts)} duct_filler cells.")
    else:
        combined = cells

    # Override induced-type volumes by sampling parent-type distribution
    rng = np.random.default_rng(seed=42)
    for child_type, parent_type in PARENT_VOLUME_TYPE.items():
        parent_vols = combined.loc[combined['type'] == parent_type, 'volume'].values
        mask = combined['type'] == child_type
        if mask.any() and len(parent_vols) > 0:
            n = int(mask.sum())
            combined.loc[mask, 'volume'] = rng.choice(parent_vols, size=n, replace=True)
            print(f"  {child_type}: volumes sampled from {parent_type} "
                  f"(n={len(parent_vols)} donors, {n} assigned)")
        elif mask.any():
            print(f"  {child_type}: no {parent_type} in this ROI — keeping QuPath volume")

    # Winsorize volumes per cell type
    n_clipped = 0
    for ct, group in combined.groupby('type'):
        if ct in WINSOR_SKIP_TYPES or len(group) < WINSOR_MIN_CELLS:
            continue
        lo = np.percentile(group['volume'], WINSOR_LO)
        hi = np.percentile(group['volume'], WINSOR_HI)
        mask = combined['type'] == ct
        before = combined.loc[mask, 'volume'].copy()
        combined.loc[mask, 'volume'] = combined.loc[mask, 'volume'].clip(lo, hi)
        n_clipped += (combined.loc[mask, 'volume'] != before).sum()
    print(f"  Winsorized {n_clipped} volumes (5th–95th pct per type)")

    # Report
    print(f"  Cell type counts:\n{combined['type'].value_counts().to_string()}")
    print(f"  Total: {len(combined)}  "
          f"x[{combined.x.min():.0f},{combined.x.max():.0f}]  "
          f"y[{combined.y.min():.0f},{combined.y.max():.0f}]")

    # Save
    OUTDIR.mkdir(parents=True, exist_ok=True)
    m = re.match(r'(JHH\d+R?)_?ROI0*(\d+)', ann_path.name)
    out_stem = f"{m.group(1)}ROI{m.group(2)}" if m else \
               ann_path.stem.replace('_for_physicell', '').replace('_summary', '')
    out_path = OUTDIR / f"{out_stem}.csv"
    combined.to_csv(out_path, index=False)
    print(f"  ✓ Saved → {out_path}")
    return out_path


# ── Interactive (GUI) entry point ─────────────────────────────────────────────
if __name__ == '__main__':
    import tkinter as tk
    from tkinter import filedialog, simpledialog

    root = tk.Tk()
    root.withdraw()

    print("Step 1: Select QuPath detection TXT file...")
    ann_path = filedialog.askopenfilename(
        title="Select QuPath detection TXT file",
        filetypes=[("Text files", "*.txt"), ("CSV files", "*.csv"), ("All files", "*.*")]
    )
    if not ann_path:
        raise SystemExit("No file selected.")
    print(f"  Selected: {Path(ann_path).name}")

    print("Step 2: Select duct annotation GeoJSON (cancel to skip)...")
    geojson_path = filedialog.askopenfilename(
        title="Select duct GeoJSON (cancel to skip)",
        filetypes=[("GeoJSON files", "*.geojson"), ("All files", "*.*")]
    )
    geojson_path = Path(geojson_path) if geojson_path else None
    print(f"  {'Selected: ' + geojson_path.name if geojson_path else 'Skipping duct_filler.'}")

    dims = simpledialog.askstring(
        "Image dimensions",
        "Enter image width,height in µm\n(default: 1100,1100)",
        initialvalue="1100,1100"
    )
    if not dims:
        raise SystemExit("No dimensions entered.")
    img_w, img_h = [float(v.strip()) for v in dims.split(',')]
    print(f"  Image: {img_w:.0f} x {img_h:.0f} µm")

    hex_r = HEX_RADIUS_DEFAULT
    if geojson_path is not None:
        hex_r = simpledialog.askfloat(
            "Hex radius", "Enter hex radius in µm (default: 8.0)",
            initialvalue=HEX_RADIUS_DEFAULT
        ) or HEX_RADIUS_DEFAULT

    process_roi(ann_path, geojson_path, img_w, img_h, hex_r, batch=False)
