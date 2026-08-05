"""
winsorize_ics_volumes.py
-------------------------
One-off post-processing: clips cell volumes in each ICS CSV to the 5th–95th
percentile of that cell type's distribution, per ROI.

Operates on: PhysiCell/config/ics/JHH_IMC/*.csv  (canonical ROI files only)
Overwrites in place.
"""

import re
import numpy as np
import pandas as pd
from pathlib import Path

ICS_DIR      = Path(__file__).parent / "PhysiCell/config/ics/JHH_IMC"
LO_PCT       = 5
HI_PCT       = 95
SKIP_TYPES   = {"duct_filler"}   # uniform volume by design — skip
SKIP_SUFFIXES = {"rectangle", "withTreg", "withductfiller", "ductfiller",
                 "withTregs", "not_movable"}

def is_canonical(path: Path) -> bool:
    return not any(s in path.stem for s in SKIP_SUFFIXES)

csvs = sorted(p for p in ICS_DIR.glob("*.csv") if is_canonical(p))
print(f"Found {len(csvs)} canonical ICS files.\n")

total_clipped = 0

for csv_path in csvs:
    df = pd.read_csv(csv_path)
    if "volume" not in df.columns:
        print(f"  SKIP {csv_path.name} — no volume column")
        continue

    n_clipped_file = 0
    for ct, group in df.groupby("type"):
        if ct in SKIP_TYPES:
            continue
        vols = group["volume"].values
        if len(vols) < 10:
            continue   # too few cells to trust percentiles

        lo = np.percentile(vols, LO_PCT)
        hi = np.percentile(vols, HI_PCT)
        mask = df["type"] == ct
        before = df.loc[mask, "volume"].copy()
        df.loc[mask, "volume"] = df.loc[mask, "volume"].clip(lower=lo, upper=hi)
        n_clipped = (df.loc[mask, "volume"] != before).sum()
        if n_clipped:
            n_clipped_file += n_clipped
            print(f"  {csv_path.stem}  {ct:40s}  "
                  f"[{lo:.1f}, {hi:.1f}]  clipped {n_clipped}")

    df.to_csv(csv_path, index=False)
    total_clipped += n_clipped_file
    if n_clipped_file == 0:
        print(f"  {csv_path.stem}  (no clipping needed)")

print(f"\nDone. {total_clipped} volumes clipped across {len(csvs)} files.")
