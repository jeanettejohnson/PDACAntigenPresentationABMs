"""
Paired initial/final antigen presentation chart for all HTANSingleCell simulations.

For each simulation, classifies tumor cells by antigen presentation status:
  class I+II  — epithelial/mesenchymal_tumor_class1_class2 (IDs 13, 16)
  class I     — epithelial/mesenchymal_tumor_class1        (IDs 12, 15)
  class II    — epithelial/mesenchymal_tumor_class2        (IDs 14, 17)
  none        — epithelial/mesenchymal_tumor + PDAC_unclassified (IDs 3, 4, 19)

Shows proportion of each class at initial and final timepoints.
Writes: analysis/antigen_presentation_initial_vs_final.png
"""

import os
import re
import glob
import numpy as np
import pandas as pd
import scipy.io
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

ROOT    = os.path.dirname(os.path.dirname(__file__))
SIM_DIR = os.path.expanduser(
    "~/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine"
    "/AntigenPresentationSimulations/HTANSingleCell"
)
CSV_PATH = os.path.join(ROOT, "assignmentsummary_HTAN_singlecell.csv")
OUT_PATH = os.path.join(ROOT, "figures", "antigen_presentation_initial_vs_final.png")

# Tumor cell type IDs → antigen class
TUMOR_CLASS = {
    3:  "none",       # epithelial_tumor
    4:  "none",       # mesenchymal_tumor
    12: "class I",    # epithelial_tumor_class1
    13: "class I+II", # epithelial_tumor_class1_class2
    14: "class II",   # epithelial_tumor_class2
    15: "class I",    # mesenchymal_tumor_class1
    16: "class I+II", # mesenchymal_tumor_class1_class2
    17: "class II",   # mesenchymal_tumor_class2
    19: "none",       # PDAC_unclassified
}

AG_CLASSES = ["class I+II", "class I", "class II", "none"]
AG_COLORS  = {
    "class I+II": "#1f78b4",
    "class I":    "#a6cee3",
    "class II":   "#e31a1c",
    "none":       "#d9d9d9",
}

def proportions_from_mat(mat_path):
    """Return dict of ag_class → proportion of tumor cells, or None if no tumor cells."""
    mat = scipy.io.loadmat(mat_path)
    type_ids = mat["cells"][5, :].astype(int)
    tumor_ids = type_ids[np.isin(type_ids, list(TUMOR_CLASS.keys()))]
    total = len(tumor_ids)
    if total == 0:
        return None
    counts = {cls: 0 for cls in AG_CLASSES}
    for tid in tumor_ids:
        counts[TUMOR_CLASS[tid]] += 1
    return {cls: counts[cls] / total for cls in AG_CLASSES}

def patient_of(sample):
    m = re.match(r"(HT\d+P?\d*)", sample)
    return m.group(1) if m else sample

def main():
    ordered_samples = list(pd.read_csv(CSV_PATH)["sample_id"])

    records = []
    for sample in ordered_samples:
        out_dir = os.path.join(SIM_DIR, sample, "output")
        init_mat = os.path.join(out_dir, "initial_cells.mat")
        final_mat = os.path.join(out_dir, "final_cells.mat")
        if not os.path.exists(init_mat):
            continue
        p_init  = proportions_from_mat(init_mat)
        p_final = proportions_from_mat(final_mat) if os.path.exists(final_mat) else None
        records.append({"sample": sample, "initial": p_init, "final": p_final})

    n = len(records)
    if n == 0:
        print("No data found.")
        return

    # Layout: two columns — initial (left) and final (right)
    bar_h  = 0.35
    gap    = 0.15   # gap between initial/final pair
    stride = 1.0    # row spacing

    fig_height = max(8, n * 0.55)
    fig, axes = plt.subplots(1, 2, figsize=(14, fig_height), sharey=True)
    ax_init, ax_final = axes

    y_ticks = []
    y_labels = []

    for i, rec in enumerate(records):
        y = (n - 1 - i) * stride
        y_ticks.append(y)
        y_labels.append(rec["sample"])

        for ax, props, title in [
            (ax_init,  rec["initial"], "Initial"),
            (ax_final, rec["final"],   "Final"),
        ]:
            if props is None:
                ax.barh(y, 1.0, height=bar_h, color="#d9d9d9", alpha=0.3)
                continue
            left = 0.0
            for cls in AG_CLASSES:
                val = props[cls]
                if val > 0:
                    ax.barh(y, val, left=left, height=bar_h,
                            color=AG_COLORS[cls], edgecolor="none")
                    left += val

    for ax, title in [(ax_init, "Initial"), (ax_final, "Final")]:
        ax.set_xlim(0, 1)
        ax.set_xlabel("Proportion of tumor cells")
        ax.set_title(title, fontsize=13, fontweight="bold")
        ax.axvline(0, color="black", linewidth=0.5)
        ax.tick_params(axis="x", labelsize=8)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xticklabels(["0", "25%", "50%", "75%", "100%"])

    ax_init.set_yticks(y_ticks)
    ax_init.set_yticklabels(y_labels, fontsize=7.5)

    # Legend
    handles = [mpatches.Patch(facecolor=AG_COLORS[cls], label=cls) for cls in AG_CLASSES]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=9,
               title="Antigen presentation", title_fontsize=9,
               bbox_to_anchor=(0.5, 0.0), framealpha=0.9)

    fig.suptitle("Tumor antigen presentation: initial vs final timepoint", fontsize=13)
    plt.tight_layout(rect=[0, 0.04, 1, 0.97])

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    plt.savefig(OUT_PATH, dpi=150, bbox_inches="tight")
    print(f"Saved: {OUT_PATH}")

if __name__ == "__main__":
    main()
