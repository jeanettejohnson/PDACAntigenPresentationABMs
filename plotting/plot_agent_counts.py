"""
Bar plot of agent_assignment counts aggregated across all CSVs in RDataObjects.
Writes: analysis/agent_assignment_counts.png
"""

import os
import glob
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_DIR     = os.path.expanduser(
    "~/OneDrive - University of Maryland School of Medicine"
    "/HTANDATA/HTANWUSTL/RDataObjects"
)
ANALYSIS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "figures")

# ── colour groups ─────────────────────────────────────────────────────────────
COLORS = {
    # T cells
    "CD8_Tcell":                            "#8B0000",
    "CD4_Tcell":                            "#FFA500",
    "Treg":                                 "#C8A2C8",
    "CD8_exhausted":                        "#FA8072",
    # Tumor
    "PDAC":                                 "#00CED1",
    "PDAC_unclassified":                    "#48D1CC",
    "epithelial":                           "#DA70D6",
    "epithelial_class1":                    "#C060C0",
    "epithelial_class2":                    "#B040B0",
    "epithelial_class1_class2":             "#9020A0",
    "mesenchymal":                          "#7FFF00",
    "mesenchymal_class1":                   "#66CC00",
    "mesenchymal_class2":                   "#4D9900",
    "mesenchymal_class1_class2":            "#336600",
    "epithelial_mesenchymal":               "#7C3AED",
    "epithelial_mesenchymal_class1":        "#6A30D0",
    "epithelial_mesenchymal_class2":        "#5820B0",
    "epithelial_mesenchymal_class1_class2": "#461090",
    "PanIN_early":                          "#b0d4f1",
    "PanIN_late":                           "#5b9bd5",
    "ADM":                                  "#2e75b6",
    "Acinar":                               "#1f4e79",
    "Acinar_REG+":                          "#16365c",
    "Duct_like_1":                          "#9dc3e6",
    "Duct_like_2":                          "#6baed6",
    # Stroma / CAF
    "CAF":                                  "#f4a460",
    "iCAF":                                 "#e8891a",
    "myCAF":                                "#c46e10",
    "apCAF":                                "#a05208",
    "CXCR4+_iCAF":                          "#d4a96a",
    "CD133+_iCAF":                          "#bf8040",
    # Immune (non-T)
    "B":                                    "#74c476",
    "Plasma":                               "#41ab5d",
    "NK":                                   "#238b45",
    "Macrophage":                           "#005a32",
    "Monocyte":                             "#a8ddb5",
    "Mast":                                 "#ccebc5",
    "Neutrophil":                           "#e5f5e0",
    "cDC1":                                 "#78c679",
    "cDC2":                                 "#31a354",
    "other_immune":                         "#bdbdbd",
    # Other tissue
    "Endothelial":                          "#9ecae1",
    "Islet":                                "#4292c6",
    "Erythrocyte":                          "#fc9272",
    "Tuft":                                 "#fdae6b",
    "other_tissue":                         "#d9d9d9",
}

DEFAULT_COLOR = "#aaaaaa"


def main():
    os.makedirs(ANALYSIS_DIR, exist_ok=True)

    # ── aggregate counts ──────────────────────────────────────────────────
    counts: dict[str, int] = {}
    csv_files = [
        f for f in glob.glob(os.path.join(DATA_DIR, "*.csv"))
        if os.path.basename(f).count("_") == 1
    ]
    print(f"Reading {len(csv_files)} CSVs (case_tissue format only) ...")

    for f in csv_files:
        with open(f) as fh:
            header = fh.readline()
        if "agent_assignment" not in header:
            continue
        df = pd.read_csv(f, usecols=["agent_assignment"])
        for val, n in df["agent_assignment"].value_counts().items():
            counts[val] = counts.get(val, 0) + int(n)

    df_counts = (
        pd.Series(counts, name="count")
          .sort_values(ascending=True)
          .reset_index()
          .rename(columns={"index": "agent"})
    )
    df_counts.columns = ["agent", "count"]

    # ── plot ──────────────────────────────────────────────────────────────
    n      = len(df_counts)
    fig_h  = max(6, n * 0.32 + 1.5)
    fig, ax = plt.subplots(figsize=(9, fig_h))

    bar_colors = [COLORS.get(a, DEFAULT_COLOR) for a in df_counts["agent"]]
    bars = ax.barh(df_counts["agent"], df_counts["count"],
                   color=bar_colors, edgecolor="white", linewidth=0.4)

    for bar, val in zip(bars, df_counts["count"]):
        ax.text(bar.get_width() + df_counts["count"].max() * 0.005,
                bar.get_y() + bar.get_height() / 2,
                f"{val:,}", va="center", ha="left", fontsize=7)

    ax.set_xlabel("Total cell count (all samples)", fontsize=9)
    ax.set_title("Agent Assignment Counts (case_tissue samples only)",
                 fontsize=11, fontweight="bold", pad=10)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax.tick_params(axis="y", labelsize=8)
    ax.tick_params(axis="x", labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(0, df_counts["count"].max() * 1.12)

    fig.tight_layout()
    out = os.path.join(ANALYSIS_DIR, "agent_assignment_counts_by_tissue.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
