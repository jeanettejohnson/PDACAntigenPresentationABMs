"""
batch_assemble_ics.py
---------------------
Runs process_roi() for every detection TXT that has a matching scaled substrate
CSV.  Image dimensions are derived from the substrate voxel grid so no manual
entry is needed.  Duct GeoJSON is used when present.

Matching is done by normalising filenames to a canonical ROI key
(e.g. JHH317ROI1) using the regex (JHH\d+R?)_?ROI0*(\d+).

Directory layout expected:
  qupath_detections/          *_detections.txt
  config/ics/ductcoordinates/ *.geojson
  config/ics/substrates/      *_ecm_scaled.csv   (domain reference)
  config/ics/JHH_IMC/         output ICS CSVs    (written by process_roi)
"""

import re
import sys
import pandas as pd
from pathlib import Path
from assemble_initial_conditions import process_roi, HEX_RADIUS_DEFAULT

BASE        = Path(__file__).parent
DET_DIR     = BASE / "qupath_detections"
GEOJSON_DIR = BASE / "PhysiCell/user_projects/antigen_presentation/config/ics/ductcoordinates"
SUB_DIR     = BASE / "PhysiCell/user_projects/antigen_presentation/config/ics/substrates"
DX          = 20.0

ROI_RE = re.compile(r'(JHH\d+R?)_?ROI0*(\d+)')

def roi_key(path):
    """Return canonical key like 'JHH317ROI1' from any matching filename."""
    m = ROI_RE.search(path.name)
    return f"{m.group(1)}ROI{m.group(2)}" if m else None

def dims_from_substrate(sub_path):
    """Derive image width/height in µm from the substrate voxel grid."""
    df = pd.read_csv(sub_path, usecols=['x', 'y'])
    img_w = (df.x.max() - df.x.min()) + DX
    img_h = (df.y.max() - df.y.min()) + DX
    return img_w, img_h

# ── Index files by ROI key ────────────────────────────────────────────────────
detections = {roi_key(p): p for p in DET_DIR.glob("*_detections.txt") if roi_key(p)}
geojsons   = {roi_key(p): p for p in GEOJSON_DIR.glob("*.geojson")    if roi_key(p)}
substrates = {roi_key(p): p for p in SUB_DIR.glob("*_ecm_scaled.csv") if roi_key(p)}

# ── Run ───────────────────────────────────────────────────────────────────────
candidates = sorted(detections.keys() & substrates.keys())
skipped    = sorted(detections.keys() - substrates.keys())

if skipped:
    print(f"No substrate found for: {', '.join(skipped)} — skipping\n")

print(f"Processing {len(candidates)} ROI(s)...\n")
errors = []

for key in candidates:
    det_path = detections[key]
    sub_path = substrates[key]
    geo_path = geojsons.get(key)          # optional
    img_w, img_h = dims_from_substrate(sub_path)

    print(f"── {key} ──────────────────────────────────────────")
    print(f"  detections : {det_path.name}")
    print(f"  substrate  : {sub_path.name}  ({img_w:.0f} x {img_h:.0f} µm)")
    print(f"  duct GeoJSON: {geo_path.name if geo_path else '(none)'}")

    try:
        process_roi(det_path, geo_path, img_w, img_h,
                    hex_r=HEX_RADIUS_DEFAULT, batch=True)
    except Exception as e:
        print(f"  ERROR: {e}")
        errors.append(key)
    print()

print(f"Done. {len(candidates) - len(errors)}/{len(candidates)} succeeded.")
if errors:
    print(f"Failed: {', '.join(errors)}")
    sys.exit(1)

# ── Generate per-ROI PhysiCell config XMLs ────────────────────────────────────
import subprocess
print("\n── Generating per-ROI PhysiCell configs ──────────────────────────────")
result = subprocess.run(
    [sys.executable, BASE / "generate_roi_configs.py"],
    check=False
)
if result.returncode != 0:
    print("  WARNING: generate_roi_configs.py exited with an error.")
