"""
Visualization: antigen presentation class proportions, initial vs final,
for all HT* simulation runs.

Reads:  analysis/antigen_class_proportions.csv
Writes: analysis/antigen_class_paired_bars.png
"""

import os
import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.gridspec import GridSpec

CSV_PATH = os.path.join(os.path.dirname(__file__), "analysis", "antigen_class_proportions.csv")
OUT_PATH = os.path.join(os.path.dirname(__file__), "analysis", "antigen_class_paired_bars.png")

AG_CLASS_ORDER  = ["class I+II", "class I", "class II", "none"]
AG_CLASS_COLORS = {
    "class I+II": "#1f78b4",
    "class I":    "#a6cee3",
    "class II":   "#e31a1c",
    "none":       "#d9d9d9",
}

PATIENT_COLORS = {
    "HT056":   "#8dd3c7",
    "HT060":   "#ffffb3",
    "HT061":   "#bebada",
    "HT061P1": "#fb8072",
    "HT064":   "#80b1d3",
    "HT071":   "#fdb462",
    "HT085":   "#b3de69",
    "HT115":   "#fccde5",
}


def patient_of(run: str) -> str:
    m = re.match(r"(HT\d+P?\d*)", run)
    return m.group(1) if m else run


def short_label(run: str) -> str:
    """Strip patient prefix to keep y-axis labels readable."""
    pat = patient_of(run)
    suffix = run[len(pat):]
    return suffix.lstrip("_") or run


def pivot_run(df: pd.DataFrame, run: str, timepoint: str) -> dict:
    sub = df[(df["run"] == run) & (df["timepoint"] == timepoint)]
    return {row["ag_class"]: row["proportion"] for _, row in sub.iterrows()}


def draw_stacked_bar(ax, y, props: dict, height=0.35, alpha=1.0, hatch=None):
    left = 0.0
    for cls in AG_CLASS_ORDER:
        val = props.get(cls, 0.0)
        ax.barh(y, val, left=left, height=height,
                color=AG_CLASS_COLORS[cls], alpha=alpha,
                edgecolor="white", linewidth=0.4,
                hatch=hatch)
        if val > 0.07:
            ax.text(left + val / 2, y, f"{val:.0%}",
                    ha="center", va="center", fontsize=5.5,
                    color="black", fontweight="normal")
        left += val


def main():
    df = pd.read_csv(CSV_PATH)

    # ── organise runs by patient ───────────────────────────────────────────
    runs = df["run"].unique().tolist()
    runs_sorted = sorted(runs, key=lambda r: (patient_of(r), r))

    patients_ordered = []
    seen = set()
    for r in runs_sorted:
        p = patient_of(r)
        if p not in seen:
            patients_ordered.append(p)
            seen.add(p)

    patient_groups = {p: [r for r in runs_sorted if patient_of(r) == p]
                      for p in patients_ordered}

    # build y-positions: small gap between runs, larger gap between patients
    y_positions = {}
    y_labels    = {}
    y           = 0
    group_spans = {}  # patient → (y_min, y_max)

    for pat in patients_ordered:
        group_runs = patient_groups[pat]
        y_start = y
        for run in group_runs:
            y_positions[run] = y
            y_labels[run]    = short_label(run)
            y += 1
        group_spans[pat] = (y_start, y - 1)
        y += 0.7   # gap between patient groups

    total_height = y

    # ── figure ────────────────────────────────────────────────────────────
    fig_h = max(8, total_height * 0.55 + 2.5)
    fig = plt.figure(figsize=(13, fig_h))

    gs = GridSpec(
        1, 3,
        figure=fig,
        width_ratios=[0.03, 1, 1],
        wspace=0.04,
        left=0.18, right=0.97, top=0.93, bottom=0.07,
    )

    ax_pat  = fig.add_subplot(gs[0])   # patient colour strip
    ax_init = fig.add_subplot(gs[1])   # initial bars
    ax_fin  = fig.add_subplot(gs[2])   # final bars

    # ── draw bars ─────────────────────────────────────────────────────────
    for run in runs_sorted:
        y_pos = y_positions[run]
        props_init  = pivot_run(df, run, "initial")
        props_final = pivot_run(df, run, "final")
        draw_stacked_bar(ax_init, y_pos, props_init,  alpha=1.0)
        draw_stacked_bar(ax_fin,  y_pos, props_final, alpha=1.0)

    # ── patient colour strip + group labels ───────────────────────────────
    for pat, (y0, y1) in group_spans.items():
        mid = (y0 + y1) / 2
        color = PATIENT_COLORS.get(pat, "#cccccc")
        ax_pat.barh(mid, 1, left=0, height=(y1 - y0 + 0.8),
                    color=color, edgecolor="none")
        ax_pat.text(0.5, mid, pat, ha="center", va="center",
                    fontsize=7, fontweight="bold", rotation=90)

    ax_pat.set_xlim(0, 1)
    ax_pat.set_ylim(-0.5, total_height)
    ax_pat.axis("off")

    # ── y-axis labels & shared formatting ─────────────────────────────────
    y_vals  = [y_positions[r] for r in runs_sorted]
    y_texts = [y_labels[r]    for r in runs_sorted]

    for ax, title in [(ax_init, "Initial"), (ax_fin, "Final")]:
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.5, total_height)
        ax.set_yticks(y_vals)
        ax.set_yticklabels(y_texts if ax is ax_init else [],
                           fontsize=7)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
        ax.set_xticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=7)
        ax.set_title(title, fontsize=10, fontweight="bold", pad=6)
        ax.tick_params(axis="y", length=0)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.grid(axis="x", color="#eeeeee", linewidth=0.6, zorder=0)

        # patient group separators
        for pat, (y0, y1) in group_spans.items():
            sep_y = y1 + 0.85
            ax.axhline(sep_y, color="#cccccc", linewidth=0.6, linestyle="--")

    # ── legend ────────────────────────────────────────────────────────────
    legend_patches = [
        mpatches.Patch(color=AG_CLASS_COLORS[c], label=c) for c in AG_CLASS_ORDER
    ]
    fig.legend(
        handles=legend_patches,
        loc="upper center",
        ncol=4,
        fontsize=8,
        frameon=False,
        title="Antigen presentation class  (proportion of all tumor cells)",
        title_fontsize=8.5,
        bbox_to_anchor=(0.58, 0.99),
    )

    fig.suptitle(
        "Antigen Presentation Class Composition — Initial vs Final Timepoint",
        fontsize=12, fontweight="bold", y=1.02,
    )

    fig.savefig(OUT_PATH, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUT_PATH}")


if __name__ == "__main__":
    main()
