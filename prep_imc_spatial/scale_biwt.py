import numpy as np
import pandas as pd
from pathlib import Path

PERCENTILE   = 95   # clip substrate values above this percentile of nonzero values
SUBSTRATES_DIR = Path(__file__).parent / "PhysiCell/user_projects/antigen_presentation/config/ics/substrates"

csvs = sorted(p for p in SUBSTRATES_DIR.glob("*.csv") if not p.stem.endswith("_scaled"))
if not csvs:
    raise SystemExit(f"No unscaled CSVs found in {SUBSTRATES_DIR}")

print(f"Found {len(csvs)} file(s) to scale in {SUBSTRATES_DIR}\n")

for input_path in csvs:
    print(f"Processing: {input_path.name}")
    df = pd.read_csv(input_path)

    n_before = len(df)
    df = df.drop_duplicates(subset=['x', 'y', 'z'], keep='first')
    n_dupes = n_before - len(df)
    if n_dupes > 0:
        print(f"  Dropped {n_dupes} duplicate voxel(s)")

    substrate_cols = [c for c in df.columns if c not in ('x', 'y', 'z')]
    for col in substrate_cols:
        vals = df[col].values.astype(float)
        nonzero_vals = vals[vals > 0]
        if len(nonzero_vals) == 0:
            print(f"  {col}: all zeros — skipping")
            continue
        upper = np.percentile(nonzero_vals, PERCENTILE)
        df[col] = (np.clip(vals, 0, upper) / upper) * 10
        print(f"  {col}: clipped at {PERCENTILE}th percentile ({upper:.3f}), scaled to [0, 10]")
        print(f"         min={df[col].min():.3f}, median={df[col].median():.3f}, max={df[col].max():.3f}")

    output_path = input_path.with_stem(input_path.stem + "_scaled")
    df.to_csv(output_path, index=False)
    print(f"  Saved: {output_path.name}\n")

print("Done.")
