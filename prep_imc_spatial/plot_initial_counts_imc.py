"""
plot_initial_counts_imc.py
---------------------------
Two clustered stacked bar charts of initial cell counts per ROI.
ROIs reordered by Ward clustering on cell type proportions.
Layout: legend → dendrogram → bars → cluster strip.

  1. duct_filler + other_tissue excluded  → initial_cell_counts_imc.png
  2. all types included                   → initial_cell_counts_imc_all_types.png

Reads:  PhysiCell/config/ics/JHH_IMC/*.csv  (canonical files only)
"""

import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from pathlib import Path

ICS_DIR    = Path(__file__).parent / "PhysiCell/config/ics/JHH_IMC"
N_CLUSTERS = 5

SKIP_SUFFIXES = {"rectangle", "withTreg", "withductfiller", "ductfiller",
                 "withTregs", "not_movable"}
ROI_RE = re.compile(r"(JHH\d+R?)ROI(\d+)")

CELL_ORDER = [
    "CAF", "apCAF",
    "CD4_Tcell", "CD8_Tcell", "Treg", "CD8_exhausted", "macrophage", "B cell",
    "epithelial_normal", "mesenchymal_normal",
    "PDAC_unclassified",
    "epithelial_tumor", "epithelial_tumor_class1",
    "epithelial_tumor_class1_class2", "epithelial_tumor_class2",
    "mesenchymal_tumor", "mesenchymal_tumor_class1",
    "mesenchymal_tumor_class1_class2", "mesenchymal_tumor_class2",
    "other_tissue", "duct_filler",
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
    "duct_filler":                     "tan",
}

LIGHT_COLORS = {"yellow", "palegreen", "lightcyan", "papayawhip", "#cccccc", "tan"}

CLUSTER_CMAP = cm.get_cmap("Set2", N_CLUSTERS)

def roi_sort_key(name):
    m = ROI_RE.match(name)
    return (m.group(1), int(m.group(2))) if m else (name, 0)

# ── Load counts ───────────────────────────────────────────────────────────────
records = {}
for p in sorted(ICS_DIR.glob("*.csv")):
    if any(s in p.stem for s in SKIP_SUFFIXES):
        continue
    m = ROI_RE.match(p.stem)
    if not m:
        continue
    roi = f"{m.group(1)}ROI{m.group(2)}"
    df  = pd.read_csv(p, usecols=["type"])
    records[roi] = df["type"].value_counts().to_dict()

if not records:
    raise SystemExit(f"No ICS files found in {ICS_DIR}")

rois = sorted(records.keys(), key=roi_sort_key)

def build_counts(exclude=None):
    order   = [t for t in CELL_ORDER if t not in (exclude or set())]
    df      = pd.DataFrame(records, index=order).fillna(0).T.loc[rois]
    present = [t for t in order if df[t].sum() > 0]
    return df[present], present

def cluster_order(counts_df):
    props  = counts_df.div(counts_df.sum(axis=1), axis=0).fillna(0)
    Z      = linkage(props.values, method="ward", metric="euclidean")
    dend   = dendrogram(Z, no_plot=True)
    ordered = [counts_df.index[i] for i in dend["leaves"]]
    labels  = fcluster(Z, N_CLUSTERS, criterion="maxclust")
    roi_cluster = {counts_df.index[i]: int(labels[i]) for i in range(len(counts_df))}
    return ordered, Z, roi_cluster

def plot_clustered(counts_df, present, ordered_rois, Z, roi_cluster, title, out_path,
                   show_dendrogram=True, show_cluster_strip=True):
    n     = len(ordered_rois)
    fig_w = max(10, n * 0.45)

    n_rows       = int(show_dendrogram) + 1 + int(show_cluster_strip)
    height_ratios = ([1.2] if show_dendrogram else []) + [4] + ([0.25] if show_cluster_strip else [])
    fig_h        = 11 if show_dendrogram else (9 if show_cluster_strip else 8)
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs  = GridSpec(n_rows, 1, figure=fig,
                   height_ratios=height_ratios,
                   hspace=0.02,
                   top=0.95, bottom=0.18, right=0.78)

    row = 0
    if show_dendrogram:
        ax_dend = fig.add_subplot(gs[row]); row += 1
    ax_bar = fig.add_subplot(gs[row]); row += 1
    if show_cluster_strip:
        ax_clust = fig.add_subplot(gs[row])

    # ── Dendrogram ────────────────────────────────────────────────────────────
    if show_dendrogram:
        dendrogram(Z, ax=ax_dend, link_color_func=lambda _: "0.4",
                   no_labels=True, above_threshold_color="0.4")
        ax_dend.set_axis_off()

    # ── Stacked bars ──────────────────────────────────────────────────────────
    x_pos  = np.arange(n)
    bottom = np.zeros(n)
    df_ord = counts_df.loc[ordered_rois]

    for ct in present:
        vals  = df_ord[ct].values
        color = CELL_COLORS.get(ct, "#cccccc")
        edge  = "black" if color in LIGHT_COLORS else "none"
        ax_bar.bar(x_pos, vals, bottom=bottom, color=color,
                   edgecolor=edge, linewidth=0.4 if edge == "black" else 0)
        bottom += vals

    ax_bar.set_xticks(x_pos)
    ax_bar.set_xticklabels([])
    ax_bar.set_ylabel("Initial cell count", fontsize=11)
    ax_bar.spines[["top", "right", "bottom"]].set_visible(False)
    ax_bar.tick_params(axis="x", which="both", length=0)

    # ── Cluster strip ─────────────────────────────────────────────────────────
    if show_cluster_strip:
        for i, roi in enumerate(ordered_rois):
            cid   = roi_cluster[roi]
            color = CLUSTER_CMAP(cid - 1)
            ax_clust.bar(i, 1, width=1.0, color=color, edgecolor="white", linewidth=0.4)
        ax_clust.set_xlim(ax_bar.get_xlim())
        ax_clust.set_ylim(0, 1)
        ax_clust.set_xticks(x_pos)
        ax_clust.set_xticklabels(ordered_rois, fontsize=10, rotation=45, ha="right")
        ax_clust.tick_params(axis="x", which="both", length=0)
        ax_clust.spines[["top", "left", "right", "bottom"]].set_visible(False)
        ax_clust.set_yticks([])
    else:
        ax_bar.set_xticks(x_pos)
        ax_bar.set_xticklabels(ordered_rois, fontsize=10, rotation=45, ha="right")

    # ── Legend above dendrogram ───────────────────────────────────────────────
    handles = [
        mpatches.Patch(facecolor=CELL_COLORS.get(ct, "#cccccc"),
                       edgecolor="black" if CELL_COLORS.get(ct) in LIGHT_COLORS else "none",
                       label=ct)
        for ct in present
    ]
    fig.legend(handles=handles, loc="center left",
               bbox_to_anchor=(0.79, 0.5),
               ncol=1, fontsize=13, frameon=False,
               handlelength=1.2)

    fig.suptitle(title, fontsize=13, fontweight="bold")

    out_path.parent.mkdir(exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")

VARIANTS = [
    # (show_dendrogram, show_cluster_strip, use_clustering, suffix)
    (True,  True,  True,  ""),
    (False, False, False, "_ordered"),
]

for show_dend, show_strip, do_cluster, suffix in VARIANTS:
    for excl, name in [
        ({"duct_filler", "other_tissue"}, "initial_cell_counts_imc"),
        ({"duct_filler"},                 "initial_cell_counts_imc_with_other"),
        (None,                            "initial_cell_counts_imc_all_types"),
    ]:
        counts_df, present = build_counts(exclude=excl)
        if do_cluster:
            ordered_rois, Z, roi_cluster = cluster_order(counts_df)
        else:
            ordered_rois = rois[:]
            _, Z, roi_cluster = cluster_order(counts_df)  # still need Z/roi_cluster as dummies

        titles = {
            "initial_cell_counts_imc":            "Agent types per IMC ROI",
            "initial_cell_counts_imc_with_other": "Agent types per IMC ROI including other tissue",
            "initial_cell_counts_imc_all_types":  "Agent types per IMC ROI all cell classes",
        }
        plot_clustered(counts_df, present, ordered_rois, Z, roi_cluster,
                       titles[name],
                       Path(__file__).parent / f"analysis/{name}{suffix}.png",
                       show_dendrogram=show_dend,
                       show_cluster_strip=show_strip)
