"""
antigen_class_comparison_imc.py
---------------------------------
All cell types in one combined stacked bar per timepoint (initial → final),
with Bezier ribbons, for the JHH IMC simulations.

Reads:  PhysiCell/outputs/JHH*/
Writes: analysis/antigen_class_comparison_imc.png
        analysis/antigen_class_proportions_imc.csv
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
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from matplotlib.gridspec import GridSpec
import pcdl

from antigen_class_comparison import TUMOR_CELL_CLASSES

OUTPUTS_DIR  = os.path.join(os.path.dirname(os.path.dirname(__file__)), "PhysiCell", "outputs")
ANALYSIS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "figures")

ROI_RE = re.compile(r"(JHH\d+R?)ROI(\d+)")

# ── Cell type order (bottom → top) and colors ─────────────────────────────────
ALL_CELL_ORDER = [
    "CAF", "apCAF",
    "CD4_Tcell", "CD8_Tcell", "Treg", "CD8_exhausted", "macrophage", "B cell",
    "epithelial_normal", "mesenchymal_normal",
    "PDAC_unclassified",
    "epithelial_tumor", "epithelial_tumor_class1",
    "epithelial_tumor_class1_class2", "epithelial_tumor_class2",
    "mesenchymal_tumor", "mesenchymal_tumor_class1",
    "mesenchymal_tumor_class1_class2", "mesenchymal_tumor_class2",
]

CELL_COLORS = {
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
    "PDAC_unclassified":               "#cccccc",
}

BAR_X = {"initial": 0.22, "final": 0.78}
BAR_W = 0.30


# ── Data helpers ──────────────────────────────────────────────────────────────

def get_jhh_folders():
    return sorted(glob.glob(os.path.join(OUTPUTS_DIR, "JHH*")))


def first_and_last_xml(folder):
    xmls = sorted(glob.glob(os.path.join(folder, "output*.xml")))
    if not xmls:
        return None, None
    return xmls[0], xmls[-1]


def load_counts(xml_path):
    """Return Series: count of each cell type (excl. duct_filler)."""
    mcds = pcdl.TimeStep(xml_path, microenv=False, graph=False, verbose=False)
    df   = mcds.get_cell_df()[["cell_type"]]
    df   = df[df["cell_type"] != "duct_filler"]
    return df["cell_type"].value_counts()


# ── Ribbon ────────────────────────────────────────────────────────────────────

def _ribbon(ax, x0_r, x1_l, y0_bot, y0_top, y1_bot, y1_top, color):
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


# ── Panel ─────────────────────────────────────────────────────────────────────

def plot_combined_panel(ax, counts_init, counts_final, run_name):
    n_init  = int(counts_init.sum())
    n_final = int(counts_final.sum())
    y_max   = max(n_init, n_final, 1)
    segments = {}

    for x_ctr, tp, counts in [
        (BAR_X["initial"], "initial", counts_init),
        (BAR_X["final"],   "final",   counts_final),
    ]:
        bottom = 0
        for ct in ALL_CELL_ORDER:
            val = int(counts.get(ct, 0))
            if val == 0:
                segments.setdefault(ct, {})[tp] = (bottom, bottom)
                continue
            top = bottom + val
            ax.bar(x_ctr, val, bottom=bottom, width=BAR_W,
                   color=CELL_COLORS.get(ct, "#cccccc"),
                   edgecolor="white", linewidth=0.4, zorder=2)
            if val > y_max * 0.04:
                ax.text(x_ctr, bottom + val / 2, f"{val:,}",
                        ha="center", va="center", fontsize=5.5,
                        color="black", zorder=3)
            segments.setdefault(ct, {})[tp] = (bottom, top)
            bottom = top

    x0_r = BAR_X["initial"] + BAR_W / 2
    x1_l = BAR_X["final"]   - BAR_W / 2
    for ct in ALL_CELL_ORDER:
        seg = segments.get(ct, {})
        if "initial" in seg and "final" in seg:
            y0b, y0t = seg["initial"]
            y1b, y1t = seg["final"]
            _ribbon(ax, x0_r, x1_l, y0b, y0t, y1b, y1t,
                    CELL_COLORS.get(ct, "#cccccc"))

    ax.set_xlim(0, 1)
    ax.set_ylim(0, y_max * 1.08)
    ax.set_xticks([BAR_X["initial"], BAR_X["final"]])
    ax.set_xticklabels(
        [f"initial\n(n={n_init:,})", f"final\n(n={n_final:,})"],
        fontsize=7,
    )
    ax.yaxis.set_major_locator(plt.MaxNLocator(4, integer=True))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.tick_params(axis="y", labelsize=6.5, labelrotation=90)
    ax.set_ylabel("cell count", fontsize=6.5, labelpad=2)
    ax.set_title(run_name, fontsize=8, fontweight="bold", pad=4)
    ax.tick_params(axis="both", length=2)
    ax.spines[["top", "right"]].set_visible(False)


# ── Grid ──────────────────────────────────────────────────────────────────────

def parse_imc_run(run_name):
    m = ROI_RE.match(run_name)
    return (m.group(1), f"ROI{m.group(2)}") if m else (run_name, run_name)


def save_grid(all_run_data):
    from collections import defaultdict
    patient_runs = defaultdict(list)
    for rd in all_run_data:
        patient, roi = parse_imc_run(rd["run_name"])
        patient_runs[patient].append((roi, rd))
    patients = sorted(patient_runs.keys())
    for p in patients:
        patient_runs[p].sort(key=lambda x: x[0])

    # one row per patient, columns = ROIs
    max_cols   = max(len(v) for v in patient_runs.values())
    n_patients = len(patients)

    PANEL_W = 2.6
    PANEL_H = 3.2
    LABEL_W = 0.5

    fig = plt.figure(figsize=(LABEL_W + max_cols * PANEL_W, n_patients * PANEL_H + 1.5))
    gs  = GridSpec(
        n_patients, max_cols + 1,
        figure=fig,
        width_ratios=[LABEL_W / PANEL_W] + [1] * max_cols,
        hspace=0.55, wspace=0.30,
        left=0.04, right=0.98,
        top=0.90, bottom=0.10,
    )

    for row_idx, patient in enumerate(patients):
        ax_lbl = fig.add_subplot(gs[row_idx, 0])
        ax_lbl.text(0.5, 0.5, patient, ha="center", va="center",
                    fontsize=9, fontweight="bold", transform=ax_lbl.transAxes)
        ax_lbl.axis("off")

        for col_idx, (roi, rd) in enumerate(patient_runs[patient]):
            ax = fig.add_subplot(gs[row_idx, col_idx + 1])
            plot_combined_panel(ax, rd["counts_init"], rd["counts_final"], roi)

        for col_idx in range(len(patient_runs[patient]), max_cols):
            fig.add_subplot(gs[row_idx, col_idx + 1]).axis("off")

    # legend — only types that actually appear
    present = set()
    for rd in all_run_data:
        present |= set(rd["counts_init"].index) | set(rd["counts_final"].index)
    legend_patches = [
        mpatches.Patch(color=CELL_COLORS.get(ct, "#cccccc"), label=ct)
        for ct in ALL_CELL_ORDER if ct in present
    ]
    fig.legend(handles=legend_patches, loc="lower center",
               ncol=min(7, len(legend_patches)),
               fontsize=7, frameon=False,
               bbox_to_anchor=(0.52, 0.0))

    fig.suptitle(
        "Cell population changes: initial → final\nIMC-initialised PhysiCell simulations",
        fontsize=11, fontweight="bold",
    )

    out = os.path.join(ANALYSIS_DIR, "antigen_class_comparison_imc.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    t0 = time.time()
    os.makedirs(ANALYSIS_DIR, exist_ok=True)
    folders = get_jhh_folders()
    print(f"Found {len(folders)} JHH* simulation folders.\n")

    all_run_data = []
    tidy_rows    = []

    for idx, folder in enumerate(folders, start=1):
        run_name = os.path.basename(folder)
        init_xml, final_xml = first_and_last_xml(folder)
        print(f"[{idx}/{len(folders)}] {run_name}")

        if not init_xml or not final_xml or init_xml == final_xml:
            print("  SKIP — insufficient output XMLs")
            continue

        try:
            counts_init  = load_counts(init_xml)
            counts_final = load_counts(final_xml)

            print(f"  init n={counts_init.sum():,}  final n={counts_final.sum():,}  "
                  f"({time.time()-t0:.1f}s)")

            for tp, counts in [("initial", counts_init), ("final", counts_final)]:
                for ct, n in counts.items():
                    tidy_rows.append({
                        "run": run_name, "timepoint": tp,
                        "cell_type": ct, "count": int(n),
                    })

            all_run_data.append({
                "run_name":    run_name,
                "counts_init": counts_init,
                "counts_final": counts_final,
            })

        except Exception as e:
            print(f"  ERROR: {e}")
            traceback.print_exc()

    if not all_run_data:
        print("No runs loaded — exiting.")
        return

    csv_out = os.path.join(ANALYSIS_DIR, "antigen_class_proportions_imc.csv")
    pd.DataFrame(tidy_rows).to_csv(csv_out, index=False)
    print(f"\nSaved: {csv_out}")

    print(f"\nGenerating grid ({len(all_run_data)} runs)...")
    save_grid(all_run_data)

    print(f"\nAll done in {time.time()-t0:.1f}s.")


if __name__ == "__main__":
    main()
