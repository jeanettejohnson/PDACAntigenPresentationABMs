"""
Compare antigen presentation class proportions at initial vs final timepoints
across all HT* simulation runs.

Only loads initial.xml and final.xml per run (fast — skips all middle timesteps).

Outputs in analysis/:
  antigen_class_proportions.csv   : tidy table of proportions per run/timepoint
  antigen_class_comparison.png    : grid of stacked bar plots, one panel per run
  antigen_class_aggregate.png     : aggregate view across all runs
"""

import os
import glob
import time
import traceback
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pcdl

OUTPUTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "PhysiCell", "outputs")
ANALYSIS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "figures")

# ── Antigen class definitions ─────────────────────────────────────────────────
# Each tumor cell type is assigned: morphology (epi/mes/unclassified) + ag class

TUMOR_CELL_CLASSES = {
    "epithelial_tumor":               ("epithelial",    "none"),
    "epithelial_tumor_class1":        ("epithelial",    "class I"),
    "epithelial_tumor_class1_class2": ("epithelial",    "class I+II"),
    "epithelial_tumor_class2":        ("epithelial",    "class II"),
    "mesenchymal_tumor":              ("mesenchymal",   "none"),
    "mesenchymal_tumor_class1":       ("mesenchymal",   "class I"),
    "mesenchymal_tumor_class1_class2":("mesenchymal",   "class I+II"),
    "mesenchymal_tumor_class2":       ("mesenchymal",   "class II"),
    "PDAC_unclassified":              ("unclassified",  "none"),
}

AG_CLASS_ORDER = ["class I+II", "class I", "class II", "none"]
MORPH_ORDER    = ["epithelial", "mesenchymal"]

AG_CLASS_COLORS = {
    "class I+II": "#00CED1",   # cyan
    "class I":    "#DA70D6",   # orchid
    "class II":   "#7FFF00",   # chartreuse
    "none":       "#d9d9d9",   # neutral gray
}


def load_timepoint(xml_path: str) -> pd.DataFrame:
    mcds = pcdl.TimeStep(xml_path, microenv=False, graph=False, verbose=False)
    df = mcds.get_cell_df()[["cell_type"]].copy()
    df["morphology"] = df["cell_type"].map(
        lambda c: TUMOR_CELL_CLASSES.get(c, (None, None))[0]
    )
    df["ag_class"] = df["cell_type"].map(
        lambda c: TUMOR_CELL_CLASSES.get(c, (None, None))[1]
    )
    return df


def compute_proportions(df: pd.DataFrame) -> pd.DataFrame:
    tumor = df[df["morphology"].notna()].copy()
    if tumor.empty:
        return pd.DataFrame()
    total = len(tumor)
    rows = []
    for ag_class in AG_CLASS_ORDER:
        n = (tumor["ag_class"] == ag_class).sum()
        rows.append({"ag_class": ag_class, "count": n, "proportion": n / total})
    return pd.DataFrame(rows)


def compute_morph_proportions(df: pd.DataFrame) -> pd.DataFrame:
    tumor = df[df["morphology"].notna()].copy()
    if tumor.empty:
        return pd.DataFrame()
    total = len(tumor)
    rows = []
    for morph in MORPH_ORDER:
        sub = tumor[tumor["morphology"] == morph]
        morph_total = len(sub)
        for ag_class in AG_CLASS_ORDER:
            n = (sub["ag_class"] == ag_class).sum()
            rows.append({
                "morphology": morph,
                "ag_class": ag_class,
                "count": n,
                "proportion_of_tumor": n / total,
                "proportion_of_morph": n / morph_total if morph_total > 0 else 0,
            })
    return pd.DataFrame(rows)


def get_ht_folders():
    return sorted(glob.glob(os.path.join(OUTPUTS_DIR, "HT*")))


BAR_X    = {"initial": 0.22, "final": 0.78}  # bar centres in axes coords
BAR_W    = 0.28


def _ribbon(ax, x0_r, x1_l, y0_bot, y0_top, y1_bot, y1_top, color):
    """Bezier ribbon connecting one class segment in the initial bar to the final bar."""
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path

    if y0_top <= y0_bot and y1_top <= y1_bot:
        return

    mx = (x0_r + x1_l) / 2
    verts = [
        (x0_r, y0_bot),
        (mx,   y0_bot),
        (mx,   y1_bot),
        (x1_l, y1_bot),
        (x1_l, y1_top),
        (mx,   y1_top),
        (mx,   y0_top),
        (x0_r, y0_top),
        (x0_r, y0_bot),
    ]
    codes = [
        Path.MOVETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.LINETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.CLOSEPOLY,
    ]
    ax.add_patch(PathPatch(Path(verts, codes),
                           facecolor=color, alpha=0.25, edgecolor="none", zorder=1))


def plot_run_panel(ax, props_initial, props_final, run_name, use_counts=False):
    """Stacked vertical bars with Bezier ribbons connecting each class."""
    segments = {}
    n_init  = int(props_initial["count"].sum())
    n_final = int(props_final["count"].sum())
    y_max   = max(n_init, n_final, 1)

    for x_ctr, tp in [(BAR_X["initial"], "initial"), (BAR_X["final"], "final")]:
        props = props_initial if tp == "initial" else props_final
        bottom = 0.0
        for ag_class in AG_CLASS_ORDER:
            row = props[props["ag_class"] == ag_class]
            if use_counts:
                val = float(row["count"].values[0]) if not row.empty else 0.0
            else:
                val = float(row["proportion"].values[0]) if not row.empty else 0.0
            top = bottom + val
            ax.bar(x_ctr, val, bottom=bottom, width=BAR_W,
                   color=AG_CLASS_COLORS[ag_class],
                   edgecolor="white", linewidth=0.5, zorder=2)
            threshold = (y_max * 0.07) if use_counts else 0.07
            if val > threshold:
                label_text = f"{int(val):,}" if use_counts else f"{val:.0%}"
                ax.text(x_ctr, bottom + val / 2, label_text,
                        ha="center", va="center", fontsize=6.5,
                        color="black", zorder=3)
            segments.setdefault(ag_class, {})[tp] = (bottom, top)
            bottom = top

    # draw ribbons
    x0_r = BAR_X["initial"] + BAR_W / 2
    x1_l = BAR_X["final"]   - BAR_W / 2
    for ag_class in AG_CLASS_ORDER:
        seg = segments.get(ag_class, {})
        if "initial" in seg and "final" in seg:
            y0_bot, y0_top = seg["initial"]
            y1_bot, y1_top = seg["final"]
            _ribbon(ax, x0_r, x1_l, y0_bot, y0_top, y1_bot, y1_top,
                    AG_CLASS_COLORS[ag_class])

    ax.set_xlim(0, 1)
    ax.set_xticks([BAR_X["initial"], BAR_X["final"]])
    ax.set_xticklabels(["initial", "final"], fontsize=7)

    if use_counts:
        ax.set_ylim(0, y_max * 1.08)
        ax.yaxis.set_major_locator(plt.MaxNLocator(3, integer=True))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{int(v):,}"))
        ax.tick_params(axis="y", labelsize=6)
        ax.set_ylabel("cell count", fontsize=5.5, labelpad=2)
    else:
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.5, 1])
        ax.set_yticklabels(["0%", "50%", "100%"], fontsize=6)

    ax.set_title(run_name, fontsize=7.5, fontweight="bold", pad=3, loc="left")
    ax.tick_params(axis="both", length=2)
    ax.spines[["top", "right"]].set_visible(False)


def parse_run(run_name: str) -> tuple[str, str, str]:
    """Return (patient, tissue, short_label) from a run name."""
    parts = run_name.split("_")
    patient = parts[0]
    tissue  = parts[1] if len(parts) > 1 else ""
    label   = "_".join(parts[1:]) if len(parts) > 1 else run_name
    return patient, tissue, label


def save_run_grid(all_run_data: list[dict], use_counts=False):
    from collections import defaultdict
    from matplotlib.gridspec import GridSpec

    # ── group and sort by patient → tissue ────────────────────────────────
    patient_runs: dict[str, list[tuple[str, dict]]] = defaultdict(list)
    for rd in all_run_data:
        patient, tissue, label = parse_run(rd["run_name"])
        patient_runs[patient].append((tissue, label, rd))

    patients = sorted(patient_runs.keys())
    for p in patients:
        patient_runs[p].sort(key=lambda x: x[0])   # sort by tissue

    max_cols  = max(len(v) for v in patient_runs.values())
    n_patients = len(patients)

    PANEL_W   = 3.0
    PANEL_H   = 2.6
    LABEL_W   = 0.7   # width of patient-label column (in panel units)

    fig = plt.figure(figsize=(LABEL_W + max_cols * PANEL_W,
                               n_patients * PANEL_H))

    gs = GridSpec(
        n_patients, max_cols + 1,
        figure=fig,
        width_ratios=[LABEL_W / PANEL_W] + [1] * max_cols,
        hspace=0.55,
        wspace=0.35,
        left=0.02, right=0.98,
        top=0.93,  bottom=0.04,
    )

    for row_idx, patient in enumerate(patients):
        # patient label in the first column
        ax_lbl = fig.add_subplot(gs[row_idx, 0])
        ax_lbl.text(0.5, 0.5, patient,
                    ha="center", va="center",
                    fontsize=9, fontweight="bold",
                    transform=ax_lbl.transAxes)
        ax_lbl.axis("off")

        runs_in_patient = patient_runs[patient]
        for col_idx, (tissue, label, rd) in enumerate(runs_in_patient):
            ax = fig.add_subplot(gs[row_idx, col_idx + 1])
            plot_run_panel(ax, rd["props_initial"], rd["props_final"], label,
                           use_counts=use_counts)

        # hide unused columns in this row
        for col_idx in range(len(runs_in_patient), max_cols):
            fig.add_subplot(gs[row_idx, col_idx + 1]).axis("off")

    # ── legend & title ────────────────────────────────────────────────────
    legend_patches = [
        mpatches.Patch(color=AG_CLASS_COLORS[c], label=c) for c in AG_CLASS_ORDER
    ]
    fig.legend(handles=legend_patches, loc="upper center", ncol=4,
               fontsize=8, frameon=False,
               title="Antigen presentation class",
               title_fontsize=8.5,
               bbox_to_anchor=(0.55, 0.99))

    fig.suptitle(
        "Tumor Antigen Presentation Dynamics",
        fontsize=11, fontweight="bold", y=1.0,
    )

    filename = "antigen_class_counts.png" if use_counts else "antigen_class_comparison.png"
    out = os.path.join(ANALYSIS_DIR, filename)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def save_aggregate_plot(all_run_data: list[dict]):
    """Dot plot showing change in antigen-presenting fraction (class I + class I+II + class II) per run."""
    rows = []
    for rd in all_run_data:
        for tp, props in [("initial", rd["props_initial"]), ("final", rd["props_final"])]:
            presenting = props[props["ag_class"] != "none"]["proportion"].sum()
            rows.append({"run": rd["run_name"], "timepoint": tp, "presenting_fraction": presenting})

    df = pd.DataFrame(rows)
    runs = [rd["run_name"] for rd in all_run_data]

    fig, ax = plt.subplots(figsize=(10, max(4, len(runs) * 0.35 + 1)))

    for run in runs:
        sub = df[df["run"] == run]
        init_val = sub[sub["timepoint"] == "initial"]["presenting_fraction"].values[0]
        final_val = sub[sub["timepoint"] == "final"]["presenting_fraction"].values[0]
        y = runs.index(run)
        ax.plot([init_val, final_val], [y, y], color="#aaaaaa", lw=1, zorder=1)
        ax.scatter(init_val, y, color="#4dac26", s=50, zorder=2, label="initial" if y == 0 else "")
        ax.scatter(final_val, y, color="#d01c8b", s=50, zorder=2, label="final" if y == 0 else "",
                   marker="D")

    ax.set_yticks(range(len(runs)))
    ax.set_yticklabels(runs, fontsize=8)
    ax.set_xlabel("Fraction of tumor cells with antigen presentation (any class)", fontsize=9)
    ax.set_xlim(-0.02, 1.02)
    ax.axvline(0.5, color="#cccccc", lw=0.8, ls="--")
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title(
        "Change in Antigen-Presenting Tumor Cell Fraction\nInitial → Final",
        fontsize=11, fontweight="bold",
    )
    fig.tight_layout()

    out = os.path.join(ANALYSIS_DIR, "antigen_class_aggregate.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def main():
    total_start = time.time()
    os.makedirs(ANALYSIS_DIR, exist_ok=True)
    folders = get_ht_folders()
    total = len(folders)
    print(f"Found {total} HT* simulation folders.")
    print(f"Loading only initial.xml and final.xml per run.\n")

    all_run_data = []
    tidy_rows = []

    for idx, folder in enumerate(folders, start=1):
        run_name = os.path.basename(folder)
        initial_xml = os.path.join(folder, "initial.xml")
        final_xml   = os.path.join(folder, "final.xml")
        print(f"[{idx}/{total}] {run_name}")

        if not os.path.exists(initial_xml) or not os.path.exists(final_xml):
            print(f"  SKIP — missing initial.xml or final.xml")
            continue

        try:
            print(f"  Loading initial.xml ...")
            df_init = load_timepoint(initial_xml)
            props_init = compute_proportions(df_init)
            morph_init = compute_morph_proportions(df_init)

            print(f"  Loading final.xml ...")
            df_final = load_timepoint(final_xml)
            props_final = compute_proportions(df_final)
            morph_final = compute_morph_proportions(df_final)

            n_tumor_init  = df_init["morphology"].notna().sum()
            n_tumor_final = df_final["morphology"].notna().sum()
            print(f"  Tumor cells — initial: {n_tumor_init}, final: {n_tumor_final}")

            for tp, props, n_tumor in [
                ("initial", props_init, n_tumor_init),
                ("final",   props_final, n_tumor_final),
            ]:
                for _, row in props.iterrows():
                    tidy_rows.append({
                        "run":        run_name,
                        "timepoint":  tp,
                        "ag_class":   row["ag_class"],
                        "count":      row["count"],
                        "proportion": row["proportion"],
                        "n_tumor":    n_tumor,
                    })

            all_run_data.append({
                "run_name":      run_name,
                "props_initial": props_init,
                "props_final":   props_final,
                "morph_initial": morph_init,
                "morph_final":   morph_final,
            })

            elapsed = time.time() - total_start
            print(f"  Done ({elapsed:.1f}s elapsed)")

        except Exception as e:
            print(f"  ERROR: {e}")
            traceback.print_exc()
            continue

    if not all_run_data:
        print("No runs loaded successfully — exiting.")
        return

    # ── Save tidy CSV ──────────────────────────────────────────────────────
    df_tidy = pd.DataFrame(tidy_rows)
    csv_path = os.path.join(ANALYSIS_DIR, "antigen_class_proportions.csv")
    df_tidy.to_csv(csv_path, index=False)
    print(f"\nSaved tidy table: {csv_path}")
    print(df_tidy.to_string(index=False))

    # ── Plots ─────────────────────────────────────────────────────────────
    wellmixed_data = [rd for rd in all_run_data if "wellmixed" in rd["run_name"]]
    print(f"\nGenerating proportions grid ({len(wellmixed_data)} wellmixed runs) ...")
    save_run_grid(wellmixed_data, use_counts=False)
    print(f"Generating counts grid ({len(wellmixed_data)} wellmixed runs) ...")
    save_run_grid(wellmixed_data, use_counts=True)

    print("Generating aggregate slope plot ...")
    save_aggregate_plot(wellmixed_data)

    total_elapsed = time.time() - total_start
    print(f"\nAll done in {total_elapsed:.1f}s. Outputs in {ANALYSIS_DIR}/")


if __name__ == "__main__":
    main()
