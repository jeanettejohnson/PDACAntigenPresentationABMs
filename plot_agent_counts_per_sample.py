"""
Stacked horizontal bar chart of agent_assignment counts, one row per sample.
Only includes CSVs with the case_tissue filename format (exactly one underscore).
Writes: analysis/agent_assignment_per_sample.png
"""

import os
import glob
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

DATA_DIR     = os.path.expanduser(
    "~/OneDrive - University of Maryland School of Medicine"
    "/HTANDATA/HTANWUSTL/RDataObjects"
)
ANALYSIS_DIR = os.path.join(os.path.dirname(__file__), "analysis")

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

# Patient ID = everything before the first underscore
def patient_of(sample: str) -> str:
    return sample.split("_")[0]


def main():
    os.makedirs(ANALYSIS_DIR, exist_ok=True)

    csv_files = sorted([
        f for f in glob.glob(os.path.join(DATA_DIR, "*.csv"))
        if os.path.basename(f).count("_") == 1
    ])
    print(f"Reading {len(csv_files)} CSVs (case_tissue format only) ...")

    # ── build per-sample counts table ─────────────────────────────────────
    rows = []
    for f in csv_files:
        sample = os.path.splitext(os.path.basename(f))[0]
        with open(f) as fh:
            header = fh.readline()
        if "agent_assignment" not in header:
            continue
        df = pd.read_csv(f, usecols=["agent_assignment"])
        vc = df["agent_assignment"].value_counts()
        for agent, n in vc.items():
            rows.append({"sample": sample, "agent": agent, "count": int(n)})

    df_long = pd.DataFrame(rows)

    # order samples: group by patient, alphabetical within patient
    samples_sorted = sorted(df_long["sample"].unique(),
                            key=lambda s: (patient_of(s), s))

    # order agents by total count descending (most abundant on the left)
    agent_order = (
        df_long.groupby("agent")["count"].sum()
               .sort_values(ascending=False)
               .index.tolist()
    )

    # pivot to wide: rows = samples, cols = agents
    df_wide = (
        df_long.pivot(index="sample", columns="agent", values="count")
               .fillna(0)
               .reindex(index=samples_sorted, columns=agent_order, fill_value=0)
    )

    total_cells = df_wide.sum(axis=1)

    # ── figure ─────────────────────────────────────────────────────────────
    n_samples   = len(samples_sorted)
    fig_h       = max(8, n_samples * 0.38 + 2.5)
    fig, ax     = plt.subplots(figsize=(13, fig_h))

    bar_colors  = [COLORS.get(a, DEFAULT_COLOR) for a in agent_order]

    lefts = [0.0] * n_samples
    y_pos = list(range(n_samples))

    for agent, color in zip(agent_order, bar_colors):
        vals = df_wide[agent].values
        bars = ax.barh(y_pos, vals, left=lefts,
                       color=color, edgecolor="white", linewidth=0.3, height=0.72)
        # label segments that are wide enough
        for bar, left, val in zip(bars, lefts, vals):
            seg_w = bar.get_width()
            row_total = total_cells[samples_sorted[y_pos.index(bar.get_y() + bar.get_height() / 2 - 0.36 + 0.36)]] if False else None
            if seg_w > total_cells.max() * 0.035 and val > 0:
                ax.text(left + seg_w / 2,
                        bar.get_y() + bar.get_height() / 2,
                        f"{int(val):,}",
                        ha="center", va="center", fontsize=4.5, color="black")
        lefts = [l + v for l, v in zip(lefts, vals)]

    # total cell count label at the right end of each bar
    for i, sample in enumerate(samples_sorted):
        ax.text(total_cells[sample] + total_cells.max() * 0.003,
                i, f"{int(total_cells[sample]):,}",
                va="center", ha="left", fontsize=6, color="#444444")

    # patient group separators
    current_patient = None
    for i, sample in enumerate(samples_sorted):
        pat = patient_of(sample)
        if pat != current_patient:
            if current_patient is not None:
                ax.axhline(i - 0.5, color="#999999", linewidth=0.7, linestyle="--")
            current_patient = pat

    ax.set_yticks(y_pos)
    ax.set_yticklabels(samples_sorted, fontsize=7)
    ax.set_xlabel("Cell count", fontsize=9)
    ax.set_title("Agent Assignment Counts per Sample (case_tissue samples only)",
                 fontsize=11, fontweight="bold", pad=10)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax.tick_params(axis="y", labelsize=7)
    ax.tick_params(axis="x", labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(0, total_cells.max() * 1.13)
    ax.set_ylim(-0.6, n_samples - 0.4)

    # ── legend ─────────────────────────────────────────────────────────────
    legend_patches = [
        mpatches.Patch(color=COLORS.get(a, DEFAULT_COLOR), label=a)
        for a in agent_order
    ]
    n_cols = min(6, len(legend_patches))
    fig.legend(handles=legend_patches, loc="upper center",
               ncol=n_cols, fontsize=6.5, frameon=False,
               bbox_to_anchor=(0.5, 1.0),
               title="Agent type", title_fontsize=7.5)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(ANALYSIS_DIR, "agent_assignment_per_sample.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
