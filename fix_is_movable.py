"""
fix_is_movable.py
-----------------
Ensures is_movable=0 for duct_filler and other_tissue in every CSV in
PhysiCell/config/ics/JHH_IMC/. Run this before launching any simulation.
"""

from pathlib import Path
import pandas as pd

IMMOVABLE = {'duct_filler', 'other_tissue'}
ICS_DIR = Path(__file__).parent / 'PhysiCell/config/ics/JHH_IMC'

csv_files = sorted(ICS_DIR.glob('*.csv'))
print(f"Found {len(csv_files)} ICS files in {ICS_DIR}\n")

for path in csv_files:
    df = pd.read_csv(path)
    df['is_movable'] = df['type'].apply(lambda t: 0 if t in IMMOVABLE else 1)
    df.to_csv(path, index=False)
    n_immovable = (df['is_movable'] == 0).sum()
    print(f"  ✓ {path.name}: {n_immovable} immovable cells fixed")

print("\nDone — all ICS files updated.")
