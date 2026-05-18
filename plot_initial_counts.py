"""
Stacked bar chart of initial cell counts for each simulation,
colored with the exact PhysiCell cell-type colors from legend.svg.

Reads initial_cells.mat from each completed simulation, counts cells
by type using row 5 (cell_type ID), and labels bars by sample_id
from assignmentsummary_HTAN_singlecell.csv.

Writes: analysis/initial_cell_counts.png
"""

import os
import glob
import numpy as np
import pandas as pd
import scipy.io
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

ROOT = os.path.dirname(__file__)
SIM_DIR = os.path.expanduser(
    "~/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine"
    "/AntigenPresentationSimulations/HTANSingleCell"
)
CSV_PATH = os.path.join(ROOT, "assignmentsummary_HTAN_singlecell.csv")
OUT_PATH = os.path.join(ROOT, "analysis", "initial_cell_counts.png")

# PhysiCell cell type ID → name (from initial.xml cellular_information/cell_types)
CELL_TYPE_NAMES = {
    0:  "epithelial_normal",
    1:  "mesenchymal_normal",
    2:  "CAF",
    3:  "epithelial_tumor",
    4:  "mesenchymal_tumor",
    5:  "other_tissue",
    6:  "CD4_Tcell",
    7:  "CD8_Tcell",
    8:  "Treg",
    9:  "CD8_exhausted",
    10: "B cell",
    11: "macrophage",
    12: "epithelial_tumor_class1",
    13: "epithelial_tumor_class1_class2",
    14: "epithelial_tumor_class2",
    15: "mesenchymal_tumor_class1",
    16: "mesenchymal_tumor_class1_class2",
    17: "mesenchymal_tumor_class2",
    18: "apCAF",
    19: "PDAC_unclassified",
}

# Colors from PhysiCell legend.svg (exact CSS color names)
CELL_COLORS = {
    "epithelial_normal":              "palegreen",
    "mesenchymal_normal":             "lightcyan",
    "CAF":                            "yellow",
    "epithelial_tumor":               "green",
    "mesenchymal_tumor":              "blue",
    "other_tissue":                   "magenta",
    "CD4_Tcell":                      "darkorange",
    "CD8_Tcell":                      "maroon",
    "Treg":                           "plum",
    "CD8_exhausted":                  "lightcoral",
    "B cell":                         "papayawhip",
    "macrophage":                     "lightpink",
    "epithelial_tumor_class1":        "chartreuse",
    "epithelial_tumor_class1_class2": "darkolivegreen",
    "epithelial_tumor_class2":        "seagreen",
    "mesenchymal_tumor_class1":       "lightskyblue",
    "mesenchymal_tumor_class1_class2":"dodgerblue",
    "mesenchymal_tumor_class2":       "royalblue",
    "apCAF":                          "grey",
    "PDAC_unclassified":              "white",
}

# Cell types that are nearly white — give them a visible edge
LIGHT_TYPES = {"B cell", "papayawhip", "PDAC_unclassified", "epithelial_normal",
               "mesenchymal_normal", "CAF"}

def sample_order():
    """Return list of sample_ids in CSV row order."""
    df = pd.read_csv(CSV_PATH)
    return list(df["sample_id"])

def count_cells(mat_path):
    """Return dict: cell_type_name → count from a PhysiCell initial_cells.mat."""
    mat = scipy.io.loadmat(mat_path)
    cell_type_ids = mat["cells"][5, :].astype(int)
    counts = {}
    for type_id, name in CELL_TYPE_NAMES.items():
        n = int(np.sum(cell_type_ids == type_id))
        if n > 0:
            counts[name] = n
    return counts

def main():
    ordered_samples = sample_order()

    # Collect data from all sims that have initial_cells.mat, in CSV order
    records = []
    for sample in ordered_samples:
        mat_path = os.path.join(SIM_DIR, sample, "output", "initial_cells.mat")
        if not os.path.exists(mat_path):
            continue
        counts = count_cells(mat_path)
        records.append({"sample": sample, **counts})

    if not records:
        print("No simulations with initial_cells.mat found.")
        return

    df = pd.DataFrame(records).fillna(0).set_index("sample")

    # Only plot cell types that appear in at least one simulation
    present_types = [t for t in CELL_TYPE_NAMES.values() if t in df.columns and df[t].sum() > 0]
    df = df[present_types]

    # Sort rows by sim_id order (already sorted above)
    n_sims = len(df)
    fig_height = max(6, n_sims * 0.45)
    fig, ax = plt.subplots(figsize=(12, fig_height))

    left = np.zeros(n_sims)
    y_pos = np.arange(n_sims)

    for cell_type in present_types:
        vals = df[cell_type].values
        color = CELL_COLORS.get(cell_type, "gray")
        edge = "black" if cell_type in LIGHT_TYPES or color in ("white", "yellow", "palegreen", "lightcyan", "papayawhip") else "none"
        linewidth = 0.4 if edge == "black" else 0
        ax.barh(y_pos, vals, left=left, color=color,
                edgecolor=edge, linewidth=linewidth, label=cell_type)
        left += vals

    ax.set_yticks(y_pos)
    ax.set_yticklabels(df.index, fontsize=8)
    ax.set_xlabel("Initial cell count")
    ax.set_title("Initial cell counts by simulation (PhysiCell colors)")
    ax.invert_yaxis()

    # Legend — one entry per cell type present, in display order
    handles = []
    for cell_type in present_types:
        color = CELL_COLORS.get(cell_type, "gray")
        edge = "black" if color in ("white", "yellow", "palegreen", "lightcyan", "papayawhip") else "none"
        handles.append(mpatches.Patch(facecolor=color, edgecolor=edge, label=cell_type))

    ax.legend(handles=handles, loc="lower right", fontsize=7,
              ncol=2, framealpha=0.9)

    plt.tight_layout()
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    plt.savefig(OUT_PATH, dpi=150, bbox_inches="tight")
    print(f"Saved: {OUT_PATH}")

if __name__ == "__main__":
    main()
