"""
T cell subset composition: initial vs final, wellmixed HT* runs.
Same paired-bar + Bezier-ribbon format as antigen_class_comparison.py.

Reads:  PhysiCell/outputs/HT*wellmixed/initial.xml & final.xml
Writes: analysis/tcell_comparison.png
        analysis/tcell_proportions.csv
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
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from matplotlib.gridspec import GridSpec
from collections import defaultdict
import pcdl

OUTPUTS_DIR  = os.path.join(os.path.dirname(__file__), "PhysiCell", "outputs")
ANALYSIS_DIR = os.path.join(os.path.dirname(__file__), "analysis")

TCELL_TYPES  = ["CD8_Tcell", "CD4_Tcell", "Treg", "CD8_exhausted"]
TCELL_COLORS = {
    "CD8_Tcell":     "#8B0000",   # dark red
    "CD4_Tcell":     "#FFA500",   # orange
    "Treg":          "#C8A2C8",   # lilac
    "CD8_exhausted": "#FA8072",   # salmon
}

BAR_X = {"initial": 0.22, "final": 0.78}
BAR_W = 0.28


def get_wellmixed_folders():
    return sorted(glob.glob(os.path.join(OUTPUTS_DIR, "HT*wellmixed")))


def load_tcell_props(xml_path: str) -> tuple[pd.DataFrame, int]:
    mcds = pcdl.TimeStep(xml_path, microenv=False, graph=False, verbose=False)
    df   = mcds.get_cell_df()[["cell_type"]]
    tcells = df[df["cell_type"].isin(TCELL_TYPES)]
    n_total = len(tcells)
    rows = []
    for t in TCELL_TYPES:
        n = (tcells["cell_type"] == t).sum()
        rows.append({
            "cell_type":  t,
            "count":      int(n),
            "proportion": n / n_total if n_total > 0 else 0.0,
        })
    return pd.DataFrame(rows), n_total


def parse_run(run_name: str) -> tuple[str, str, str]:
    parts   = run_name.split("_")
    patient = parts[0]
    tissue  = parts[1] if len(parts) > 1 else ""
    label   = "_".join(parts[1:]) if len(parts) > 1 else run_name
    return patient, tissue, label


# ── drawing helpers ───────────────────────────────────────────────────────────

def _ribbon(ax, x0_r, x1_l, y0_bot, y0_top, y1_bot, y1_top, color):
    if y0_top <= y0_bot and y1_top <= y1_bot:
        return
    mx    = (x0_r + x1_l) / 2
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
                           facecolor=color, alpha=0.25,
                           edgecolor="none", zorder=1))


def plot_run_panel(ax, props_init, props_final, n_init, n_final, label):
    segments = {}
    y_max = max(n_init, n_final, 1)

    for x_ctr, tp, props in [
        (BAR_X["initial"], "initial", props_init),
        (BAR_X["final"],   "final",   props_final),
    ]:
        bottom = 0.0
        for ct in TCELL_TYPES:
            row = props[props["cell_type"] == ct]
            val = float(row["count"].values[0]) if not row.empty else 0.0
            top = bottom + val
            ax.bar(x_ctr, val, bottom=bottom, width=BAR_W,
                   color=TCELL_COLORS[ct],
                   edgecolor="white", linewidth=0.5, zorder=2)
            if val / y_max > 0.07:
                ax.text(x_ctr, bottom + val / 2, f"{int(val):,}",
                        ha="center", va="center",
                        fontsize=6, color="black", zorder=3)
            segments.setdefault(ct, {})[tp] = (bottom, top)
            bottom = top

    x0_r = BAR_X["initial"] + BAR_W / 2
    x1_l = BAR_X["final"]   - BAR_W / 2
    for ct in TCELL_TYPES:
        seg = segments.get(ct, {})
        if "initial" in seg and "final" in seg:
            y0_bot, y0_top = seg["initial"]
            y1_bot, y1_top = seg["final"]
            _ribbon(ax, x0_r, x1_l, y0_bot, y0_top, y1_bot, y1_top,
                    TCELL_COLORS[ct])

    ax.set_xlim(0, 1)
    ax.set_ylim(0, y_max * 1.08)
    ax.set_xticks([BAR_X["initial"], BAR_X["final"]])
    ax.set_xticklabels(["initial", "final"], fontsize=7)
    # y-axis: 3 evenly spaced ticks, formatted as integers
    ax.yaxis.set_major_locator(plt.MaxNLocator(3, integer=True))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(
        lambda v, _: f"{int(v):,}"
    ))
    ax.tick_params(axis="y", labelsize=6)
    ax.set_ylabel("cell count", fontsize=5.5, labelpad=2)
    ax.set_title(label, fontsize=7.5, fontweight="bold", pad=3, loc="left")
    ax.tick_params(axis="x", length=2)
    ax.spines[["top", "right"]].set_visible(False)


def save_grid(all_run_data: list[dict]):
    patient_runs = defaultdict(list)
    for rd in all_run_data:
        patient, tissue, label = parse_run(rd["run_name"])
        patient_runs[patient].append((tissue, label, rd))

    patients = sorted(patient_runs.keys())
    for p in patients:
        patient_runs[p].sort(key=lambda x: x[0])

    max_cols   = max(len(v) for v in patient_runs.values())
    n_patients = len(patients)
    PANEL_W, PANEL_H, LABEL_W = 3.0, 2.6, 0.7

    fig = plt.figure(figsize=(LABEL_W + max_cols * PANEL_W, n_patients * PANEL_H))
    gs  = GridSpec(
        n_patients, max_cols + 1, figure=fig,
        width_ratios=[LABEL_W / PANEL_W] + [1] * max_cols,
        hspace=0.55, wspace=0.35,
        left=0.02, right=0.98, top=0.93, bottom=0.04,
    )

    for row_idx, patient in enumerate(patients):
        ax_lbl = fig.add_subplot(gs[row_idx, 0])
        ax_lbl.text(0.5, 0.5, patient, ha="center", va="center",
                    fontsize=9, fontweight="bold", transform=ax_lbl.transAxes)
        ax_lbl.axis("off")

        runs_in_patient = patient_runs[patient]
        for col_idx, (tissue, label, rd) in enumerate(runs_in_patient):
            ax = fig.add_subplot(gs[row_idx, col_idx + 1])
            plot_run_panel(ax,
                           rd["props_initial"], rd["props_final"],
                           rd["n_initial"],     rd["n_final"],
                           label)

        for col_idx in range(len(runs_in_patient), max_cols):
            fig.add_subplot(gs[row_idx, col_idx + 1]).axis("off")

    legend_patches = [
        mpatches.Patch(color=TCELL_COLORS[t], label=t) for t in TCELL_TYPES
    ]
    fig.legend(handles=legend_patches, loc="upper center", ncol=4,
               fontsize=8, frameon=False,
               title="T cell subset", title_fontsize=8.5,
               bbox_to_anchor=(0.55, 0.99))

    fig.suptitle("T Cell Subset Dynamics: Initial vs Final",
                 fontsize=11, fontweight="bold", y=1.0)

    out = os.path.join(ANALYSIS_DIR, "tcell_comparison.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def main():
    total_start = time.time()
    os.makedirs(ANALYSIS_DIR, exist_ok=True)
    folders = get_wellmixed_folders()
    total   = len(folders)
    print(f"Found {total} HT*wellmixed folders.\n")

    all_run_data = []
    tidy_rows    = []

    for idx, folder in enumerate(folders, start=1):
        run_name    = os.path.basename(folder)
        initial_xml = os.path.join(folder, "initial.xml")
        final_xml   = os.path.join(folder, "final.xml")
        print(f"[{idx}/{total}] {run_name}")

        if not os.path.exists(initial_xml) or not os.path.exists(final_xml):
            print(f"  SKIP — missing initial.xml or final.xml")
            continue

        try:
            print(f"  Loading initial.xml ...")
            props_init, n_init = load_tcell_props(initial_xml)

            print(f"  Loading final.xml ...")
            props_final, n_final = load_tcell_props(final_xml)

            print(f"  T cells — initial: {n_init}, final: {n_final}")

            for tp, props, n in [("initial", props_init, n_init),
                                  ("final",   props_final, n_final)]:
                for _, row in props.iterrows():
                    tidy_rows.append({
                        "run":        run_name,
                        "timepoint":  tp,
                        "cell_type":  row["cell_type"],
                        "count":      row["count"],
                        "proportion": row["proportion"],
                        "n_tcell":    n,
                    })

            all_run_data.append({
                "run_name":     run_name,
                "props_initial": props_init,
                "props_final":   props_final,
                "n_initial":     n_init,
                "n_final":       n_final,
            })

            print(f"  Done ({time.time() - total_start:.1f}s elapsed)")

        except Exception as e:
            print(f"  ERROR: {e}")
            traceback.print_exc()

    if not all_run_data:
        print("No runs loaded — exiting.")
        return

    csv_path = os.path.join(ANALYSIS_DIR, "tcell_proportions.csv")
    pd.DataFrame(tidy_rows).to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")

    print("\nGenerating T cell grid plot ...")
    save_grid(all_run_data)

    total_elapsed = time.time() - total_start
    print(f"\nAll done in {total_elapsed:.1f}s.")


if __name__ == "__main__":
    main()
