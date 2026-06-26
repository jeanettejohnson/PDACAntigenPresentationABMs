"""
plot_cell_size_by_type.py
--------------------------
Plots max diameter distributions per cell type across all QuPath detection
files, using the same type mapping as assemble_initial_conditions.py.
Output: cell_size_by_type.png
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path
from assemble_initial_conditions import TYPE_MAP, _lookup_type

DET_DIR = Path(__file__).parent / "qupath_detections"
OUT     = Path(__file__).parent / "analysis/cell_size_by_type.png"

# ── Load all detections ───────────────────────────────────────────────────────
dfs = []
for p in sorted(DET_DIR.glob("*_detections.txt")):
    df = pd.read_table(p, usecols=["Classification", "Max diameter µm"])
    df.columns = ["classification", "max_diam"]
    dfs.append(df)

if not dfs:
    raise SystemExit(f"No detection files found in {DET_DIR}")

all_cells = pd.concat(dfs, ignore_index=True)
print(f"Loaded {len(all_cells):,} cells from {len(dfs)} files")

# Map to PhysiCell types; drop unmapped
all_cells["type"] = all_cells["classification"].map(_lookup_type)
n_before = len(all_cells)
all_cells = all_cells.dropna(subset=["type"])
all_cells = all_cells[all_cells["type"] != "UNKNOWNTYPE"]
print(f"Retained {len(all_cells):,} cells after type mapping (dropped {n_before - len(all_cells):,})")

# ── Order types by median diameter ───────────────────────────────────────────
order = (all_cells.groupby("type")["max_diam"]
         .median()
         .sort_values(ascending=False)
         .index.tolist())

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(12, 5))

# PhysiCell paint_by_number colors, indexed by cell type definition order
PHYSICELL_COLORS = {
    "epithelial_normal":               "palegreen",
    "mesenchymal_normal":              "lightcyan",
    "CAF":                             "yellow",
    "epithelial_tumor":                "green",
    "mesenchymal_tumor":               "blue",
    "other_tissue":                    "magenta",
    "CD4_Tcell":                       "darkorange",
    "CD8_Tcell":                       "maroon",
    "Treg":                            "plum",
    "CD8_exhausted":                   "lightcoral",
    "B cell":                          "papayawhip",
    "macrophage":                      "lightpink",
    "epithelial_tumor_class1":         "chartreuse",
    "epithelial_tumor_class1_class2":  "darkolivegreen",
    "epithelial_tumor_class2":         "seagreen",
    "mesenchymal_tumor_class1":        "lightskyblue",
    "mesenchymal_tumor_class1_class2": "dodgerblue",
    "mesenchymal_tumor_class2":        "royalblue",
    "apCAF":                           "grey",
    "PDAC_unclassified":               "lightgrey",  # white → lightgrey for visibility
    "duct_filler":                     "tan",
}

data_by_type = [all_cells.loc[all_cells["type"] == t, "max_diam"].values for t in order]
counts       = [len(d) for d in data_by_type]
face_colors  = [PHYSICELL_COLORS.get(t, "#4C9BE8") for t in order]

bp = ax.boxplot(data_by_type, positions=range(len(order)),
                patch_artist=True, showfliers=False,
                medianprops=dict(color="black", linewidth=1.5),
                whiskerprops=dict(linewidth=1),
                capprops=dict(linewidth=1))

for patch, color in zip(bp["boxes"], face_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.85)

ax.set_xticks(range(len(order)))
ax.set_xticklabels(
    [f"{t}\n(n={c:,})" for t, c in zip(order, counts)],
    rotation=40, ha="right", fontsize=8
)
ax.set_ylabel("Max diameter (µm)")
ax.set_title("Cell size by type across all ROIs")
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.grid(axis="y", linewidth=0.5, alpha=0.5)
ax.spines[["top", "right"]].set_visible(False)

fig.tight_layout()
fig.savefig(OUT, dpi=150)
print(f"Saved → {OUT}")
