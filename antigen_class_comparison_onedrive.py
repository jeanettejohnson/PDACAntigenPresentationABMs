"""
Compare antigen presentation class proportions at initial vs final timepoints
for all HTANSingleCell simulations in OneDrive.

Reads initial.xml and final.xml from each simulation's output/ subfolder.

Outputs in analysis/:
  antigen_class_comparison_htan.png  : grid of stacked bar plots, one panel per run
  antigen_class_aggregate_htan.png   : aggregate slope plot
"""

import os
import re
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

SIM_DIR      = os.path.expanduser(
    "~/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine"
    "/AntigenPresentationSimulations/HTANSingleCell"
)
ROOT         = os.path.dirname(__file__)
ANALYSIS_DIR = os.path.join(ROOT, "analysis")
CSV_PATH     = os.path.join(ROOT, "assignmentsummary_HTAN_singlecell.csv")

TCELL_TYPES = {
    "CD8_Tcell":      "#6B0000",   # maroon
    "CD4_Tcell":      "#FF8C00",   # darkorange
    "Treg":           "#C8A2C8",   # plum
    "CD8_exhausted":  "#F08080",   # lightcoral
}
TCELL_ORDER = ["CD8_Tcell", "CD4_Tcell", "Treg", "CD8_exhausted"]

TUMOR_CELL_CLASSES = {
    "epithelial_tumor":                ("epithelial",   "none"),
    "epithelial_tumor_class1":         ("epithelial",   "class I"),
    "epithelial_tumor_class1_class2":  ("epithelial",   "class I+II"),
    "epithelial_tumor_class2":         ("epithelial",   "class II"),
    "mesenchymal_tumor":               ("mesenchymal",  "none"),
    "mesenchymal_tumor_class1":        ("mesenchymal",  "class I"),
    "mesenchymal_tumor_class1_class2": ("mesenchymal",  "class I+II"),
    "mesenchymal_tumor_class2":        ("mesenchymal",  "class II"),
    "PDAC_unclassified":               ("unclassified", "none"),
}

AG_CLASS_ORDER  = ["class I+II", "class I", "class II", "none"]
AG_CLASS_COLORS = {
    "class I+II": "#00CED1",
    "class I":    "#DA70D6",
    "class II":   "#7FFF00",
    "none":       "#d9d9d9",
}

BAR_X = {"initial": 0.22, "final": 0.78}
BAR_W = 0.28


def get_sim_folders():
    ordered = list(pd.read_csv(CSV_PATH)["sample_id"])
    folders = []
    for s in ordered:
        p = os.path.join(SIM_DIR, s)
        if os.path.isdir(p):
            folders.append(p)
    return folders


def load_timepoint(xml_path):
    mcds = pcdl.TimeStep(xml_path, microenv=False, graph=False, verbose=False)
    df = mcds.get_cell_df()[["cell_type"]].copy()
    df["morphology"] = df["cell_type"].map(
        lambda c: TUMOR_CELL_CLASSES.get(c, (None, None))[0]
    )
    df["ag_class"] = df["cell_type"].map(
        lambda c: TUMOR_CELL_CLASSES.get(c, (None, None))[1]
    )
    return df


def compute_tcell_counts(df):
    """Return dict: tcell_type → count."""
    return {t: int((df["cell_type"] == t).sum()) for t in TCELL_ORDER}


def compute_proportions(df):
    tumor = df[df["morphology"].notna()].copy()
    if tumor.empty:
        return pd.DataFrame()
    total = len(tumor)
    rows = []
    for ag_class in AG_CLASS_ORDER:
        n = (tumor["ag_class"] == ag_class).sum()
        rows.append({"ag_class": ag_class, "count": n, "proportion": n / total})
    return pd.DataFrame(rows)


def _ribbon(ax, x0_r, x1_l, y0_bot, y0_top, y1_bot, y1_top, color):
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path
    if y0_top <= y0_bot and y1_top <= y1_bot:
        return
    mx = (x0_r + x1_l) / 2
    verts = [
        (x0_r, y0_bot), (mx, y0_bot), (mx, y1_bot), (x1_l, y1_bot),
        (x1_l, y1_top), (mx, y1_top), (mx, y0_top), (x0_r, y0_top),
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


def plot_tcell_panel(ax, tc_init, tc_final, run_name):
    """Grouped bars: initial and final T cell counts, stacked by T cell type."""
    for x_ctr, tc in [(BAR_X["initial"], tc_init), (BAR_X["final"], tc_final)]:
        bottom = 0
        for ttype in TCELL_ORDER:
            val = tc.get(ttype, 0)
            ax.bar(x_ctr, val, bottom=bottom, width=BAR_W,
                   color=TCELL_TYPES[ttype], edgecolor="white", linewidth=0.5, zorder=2)
            bottom += val

    y_max = max(sum(tc_init.values()), sum(tc_final.values()), 1)
    ax.set_xlim(0, 1)
    ax.set_xticks([BAR_X["initial"], BAR_X["final"]])
    ax.set_xticklabels(["init", "final"], fontsize=6)
    ax.set_ylim(0, y_max * 1.12)
    ax.yaxis.set_major_locator(plt.MaxNLocator(3, integer=True))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.tick_params(axis="y", labelsize=5.5)
    ax.set_ylabel("T cells", fontsize=5.5, labelpad=2)
    ax.tick_params(axis="both", length=2)
    ax.spines[["top", "right"]].set_visible(False)


def plot_run_panel(ax, props_initial, props_final, run_name, use_counts=False):
    segments = {}
    n_init  = int(props_initial["count"].sum()) if not props_initial.empty else 0
    n_final = int(props_final["count"].sum())   if not props_final.empty  else 0
    y_max   = max(n_init, n_final, 1)

    for x_ctr, tp, props in [
        (BAR_X["initial"], "initial", props_initial),
        (BAR_X["final"],   "final",   props_final),
    ]:
        if props.empty:
            ax.text(x_ctr, 0.5, "no tumor\ncells",
                    ha="center", va="center", fontsize=6, color="#aaaaaa",
                    transform=ax.transAxes)
            continue
        bottom = 0.0
        for ag_class in AG_CLASS_ORDER:
            row = props[props["ag_class"] == ag_class]
            val = float(row["count"].values[0] if use_counts else row["proportion"].values[0]) \
                  if not row.empty else 0.0
            top = bottom + val
            ax.bar(x_ctr, val, bottom=bottom, width=BAR_W,
                   color=AG_CLASS_COLORS[ag_class],
                   edgecolor="white", linewidth=0.5, zorder=2)
            threshold = (y_max * 0.07) if use_counts else 0.07
            if val > threshold:
                label = f"{int(val):,}" if use_counts else f"{val:.0%}"
                ax.text(x_ctr, bottom + val / 2, label,
                        ha="center", va="center", fontsize=6.5, color="black", zorder=3)
            segments.setdefault(ag_class, {})[tp] = (bottom, top)
            bottom = top

    x0_r = BAR_X["initial"] + BAR_W / 2
    x1_l = BAR_X["final"]   - BAR_W / 2
    for ag_class in AG_CLASS_ORDER:
        seg = segments.get(ag_class, {})
        if "initial" in seg and "final" in seg:
            y0b, y0t = seg["initial"]
            y1b, y1t = seg["final"]
            _ribbon(ax, x0_r, x1_l, y0b, y0t, y1b, y1t, AG_CLASS_COLORS[ag_class])

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


def parse_run(run_name):
    m = re.match(r"(HT\d+P?\d*)", run_name)
    patient = m.group(1) if m else run_name
    label = run_name[len(patient):].lstrip("_") or run_name
    return patient, label


def save_run_grid(all_run_data, use_counts=False):
    from collections import defaultdict
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

    patient_runs = defaultdict(list)
    for rd in all_run_data:
        patient, label = parse_run(rd["run_name"])
        patient_runs[patient].append((label, rd))

    patients   = list(dict.fromkeys(parse_run(rd["run_name"])[0] for rd in all_run_data))
    max_cols   = max(len(v) for v in patient_runs.values())
    n_patients = len(patients)
    PANEL_W, PANEL_H, LABEL_W = 3.0, 4.0, 0.7

    fig = plt.figure(figsize=(LABEL_W + max_cols * PANEL_W, n_patients * PANEL_H))
    gs  = GridSpec(
        n_patients, max_cols + 1, figure=fig,
        width_ratios=[LABEL_W / PANEL_W] + [1] * max_cols,
        hspace=0.6, wspace=0.35,
        left=0.02, right=0.98, top=0.93, bottom=0.04,
    )

    for row_idx, patient in enumerate(patients):
        ax_lbl = fig.add_subplot(gs[row_idx, 0])
        ax_lbl.text(0.5, 0.5, patient, ha="center", va="center",
                    fontsize=9, fontweight="bold", transform=ax_lbl.transAxes)
        ax_lbl.axis("off")

        for col_idx, (label, rd) in enumerate(patient_runs[patient]):
            # nested 2-row GridSpec: top = tumor (65%), bottom = T cells (35%)
            inner = GridSpecFromSubplotSpec(
                2, 1, subplot_spec=gs[row_idx, col_idx + 1],
                height_ratios=[2, 1], hspace=0.45,
            )
            ax_tumor  = fig.add_subplot(inner[0])
            ax_tcells = fig.add_subplot(inner[1])

            plot_run_panel(ax_tumor, rd["props_initial"], rd["props_final"],
                           label, use_counts=use_counts)
            plot_tcell_panel(ax_tcells, rd["tc_initial"], rd["tc_final"], label)

        for col_idx in range(len(patient_runs[patient]), max_cols):
            fig.add_subplot(gs[row_idx, col_idx + 1]).axis("off")

    tumor_patches = [
        mpatches.Patch(color=AG_CLASS_COLORS[c], label=c) for c in AG_CLASS_ORDER
    ]
    tcell_patches = [
        mpatches.Patch(color=TCELL_TYPES[t], label=t.replace("_", " ")) for t in TCELL_ORDER
    ]
    fig.legend(handles=tumor_patches, loc="upper left", ncol=4,
               fontsize=7.5, frameon=False,
               title="Antigen presentation class", title_fontsize=8,
               bbox_to_anchor=(0.08, 0.995))
    fig.legend(handles=tcell_patches, loc="upper right", ncol=4,
               fontsize=7.5, frameon=False,
               title="T cell type", title_fontsize=8,
               bbox_to_anchor=(0.92, 0.995))
    fig.suptitle("Tumor Antigen Presentation & T Cell Dynamics — HTANSingleCell",
                 fontsize=11, fontweight="bold", y=1.0)

    suffix = "counts" if use_counts else "proportions"
    out = os.path.join(ANALYSIS_DIR, f"antigen_class_comparison_htan_{suffix}.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def save_aggregate_plot(all_run_data):
    rows = []
    for rd in all_run_data:
        for tp, props in [("initial", rd["props_initial"]), ("final", rd["props_final"])]:
            if props.empty:
                presenting = 0.0
            else:
                presenting = props[props["ag_class"] != "none"]["proportion"].sum()
            rows.append({"run": rd["run_name"], "timepoint": tp,
                         "presenting_fraction": presenting})

    df   = pd.DataFrame(rows)
    runs = [rd["run_name"] for rd in all_run_data]

    fig, ax = plt.subplots(figsize=(10, max(4, len(runs) * 0.35 + 1)))
    for y, run in enumerate(runs):
        sub = df[df["run"] == run]
        iv  = sub[sub["timepoint"] == "initial"]["presenting_fraction"].values[0]
        fv  = sub[sub["timepoint"] == "final"]["presenting_fraction"].values[0]
        ax.plot([iv, fv], [y, y], color="#aaaaaa", lw=1, zorder=1)
        ax.scatter(iv, y, color="#4dac26", s=50, zorder=2,
                   label="initial" if y == 0 else "")
        ax.scatter(fv, y, color="#d01c8b", s=50, zorder=2, marker="D",
                   label="final" if y == 0 else "")

    ax.set_yticks(range(len(runs)))
    ax.set_yticklabels(runs, fontsize=7)
    ax.set_xlabel("Fraction of tumor cells with antigen presentation (any class)")
    ax.set_xlim(-0.02, 1.02)
    ax.axvline(0.5, color="#cccccc", lw=0.8, ls="--")
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title("Change in Antigen-Presenting Tumor Cell Fraction\nInitial → Final",
                 fontsize=11, fontweight="bold")
    fig.tight_layout()

    out = os.path.join(ANALYSIS_DIR, "antigen_class_aggregate_htan.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def main():
    t0 = time.time()
    os.makedirs(ANALYSIS_DIR, exist_ok=True)
    folders = get_sim_folders()
    print(f"Found {len(folders)} simulations.\n")

    all_run_data = []
    for idx, folder in enumerate(folders, 1):
        run_name    = os.path.basename(folder)
        initial_xml = os.path.join(folder, "output", "initial.xml")
        final_xml   = os.path.join(folder, "output", "final.xml")
        print(f"[{idx}/{len(folders)}] {run_name}", end="  ", flush=True)

        if not os.path.exists(initial_xml):
            print("SKIP — missing initial.xml")
            continue

        try:
            df_init    = load_timepoint(initial_xml)
            props_init = compute_proportions(df_init)
            tc_init    = compute_tcell_counts(df_init)

            if os.path.exists(final_xml):
                df_final    = load_timepoint(final_xml)
                props_final = compute_proportions(df_final)
                tc_final    = compute_tcell_counts(df_final)
            else:
                df_final    = pd.DataFrame(columns=["cell_type"])
                props_final = pd.DataFrame()
                tc_final    = {t: 0 for t in TCELL_ORDER}

            n_i = df_init["morphology"].notna().sum()
            n_f = df_final["morphology"].notna().sum() if not df_final.empty else 0
            print(f"tumor cells: {n_i} → {n_f}  ({time.time()-t0:.0f}s)")

            all_run_data.append({
                "run_name":      run_name,
                "props_initial": props_init,
                "props_final":   props_final,
                "tc_initial":    tc_init,
                "tc_final":      tc_final,
            })
        except Exception as e:
            print(f"ERROR: {e}")
            traceback.print_exc()

    if not all_run_data:
        print("No runs loaded — exiting.")
        return

    print(f"\nGenerating proportions grid ({len(all_run_data)} runs)...")
    save_run_grid(all_run_data, use_counts=False)
    print(f"Generating counts grid...")
    save_run_grid(all_run_data, use_counts=True)
    print(f"Generating aggregate slope plot...")
    save_aggregate_plot(all_run_data)
    print(f"\nAll done in {time.time()-t0:.1f}s.")


if __name__ == "__main__":
    main()
