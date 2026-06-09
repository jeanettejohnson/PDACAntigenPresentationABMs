import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import filedialog

PERCENTILE = 95  # clip ECM values above this percentile of nonzero values

# ── Prompt user to select BIWT CSV ───────────────────────────────────────────
root = tk.Tk()
root.withdraw()
input_path = filedialog.askopenfilename(
    title="Select BIWT substrate CSV",
    filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
)
if not input_path:
    raise SystemExit("No file selected.")
print(f"Selected: {input_path}")

# ── Load and scale ────────────────────────────────────────────────────────────
df = pd.read_csv(input_path)

ecm_col = [c for c in df.columns if c not in ['x', 'y', 'z']]
print(f"Substrate columns to scale: {ecm_col}")

for col in ecm_col:
    vals = df[col].values.astype(float)
    nonzero_vals = vals[vals > 0]
    upper = np.percentile(nonzero_vals, PERCENTILE)
    df[col] = (np.clip(vals, 0, upper) / upper) * 10
    print(f"  {col}: clipped at {PERCENTILE}th percentile ({upper:.3f}), scaled to [0, 10]")
    print(f"         min={df[col].min():.3f}, median={df[col].median():.3f}, max={df[col].max():.3f}")

# ── Save to same location with _scaled suffix ─────────────────────────────────
output_path = input_path.replace(".csv", "_scaled.csv")
df.to_csv(output_path, index=False)
print(f"\nSaved to: {output_path}")
