"""
plot_apcaf_vs_tcells.py
------------------------
Scatter: CAF-family type vs total T cell count per ROI initial condition.
T cells = CD4_Tcell + CD8_Tcell + Treg + CD8_exhausted

Reads:  PhysiCell/config/ics/JHH_IMC/*.csv
Writes: analysis/apcaf_vs_tcells_ic.png
        analysis/caf_vs_tcells_ic.png
"""

import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pathlib import Path
from adjustText import adjust_text

ICS_DIR = Path(__file__).parent.parent / "PhysiCell/config/ics/JHH_IMC"

SKIP_SUFFIXES = {"rectangle", "withTreg", "withductfiller", "ductfiller",
                 "withTregs", "not_movable"}

T_CELL_TYPES = {"CD4_Tcell", "CD8_Tcell", "Treg", "CD8_exhausted"}
EPI_TYPES    = {"epithelial_tumor", "epithelial_tumor_class1",
                "epithelial_tumor_class1_class2", "epithelial_tumor_class2"}
MES_TYPES    = {"mesenchymal_tumor", "mesenchymal_tumor_class1",
                "mesenchymal_tumor_class1_class2", "mesenchymal_tumor_class2"}

ROI_RE = re.compile(r"(JHH\d+R?)ROI(\d+)")

def roi_sort_key(name):
    m = ROI_RE.match(name)
    return (m.group(1), int(m.group(2))) if m else (name, 0)

# ── Load counts (once) ────────────────────────────────────────────────────────
rows = []
for p in sorted(ICS_DIR.glob("*.csv")):
    if any(s in p.stem for s in SKIP_SUFFIXES):
        continue
    m = ROI_RE.match(p.stem)
    if not m:
        continue
    roi     = f"{m.group(1)}ROI{m.group(2)}"
    patient = m.group(1)
    df      = pd.read_csv(p, usecols=["type"])
    vc      = df["type"].value_counts()
    rows.append({
        "roi":     roi,
        "patient": patient,
        "apCAF":   int(vc.get("apCAF", 0)),
        "CD8":     int(vc.get("CD8_Tcell", 0)),
        "CD4":     int(vc.get("CD4_Tcell", 0)),
        "CAF":     int(vc.get("CAF", 0)),
        "tcells":  int(sum(vc.get(t, 0) for t in T_CELL_TYPES)),
        "epi":              int(sum(vc.get(t, 0) for t in EPI_TYPES)),
        "mes":              int(sum(vc.get(t, 0) for t in MES_TYPES)),
        "epi_class1":       int(vc.get("epithelial_tumor_class1", 0)),
        "epi_class1_class2":int(vc.get("epithelial_tumor_class1_class2", 0)),
        "mes_class1":       int(vc.get("mesenchymal_tumor_class1", 0)),
        "mes_class1_class2":int(vc.get("mesenchymal_tumor_class1_class2", 0)),
        "class1":           int(vc.get("epithelial_tumor_class1", 0) + vc.get("mesenchymal_tumor_class1", 0)),
        "class1_class2":    int(vc.get("epithelial_tumor_class1_class2", 0) + vc.get("mesenchymal_tumor_class1_class2", 0)),
    })

data = pd.DataFrame(rows).sort_values("roi", key=lambda s: s.map(roi_sort_key))

patients = sorted(data["patient"].unique())
cmap     = cm.get_cmap("tab20", len(patients))
pal      = {p: cmap(i) for i, p in enumerate(patients)}


# ── Plotting function ─────────────────────────────────────────────────────────
LABELS = {
    "apCAF":  "apCAF count",
    "CAF":    "CAF count",
    "tcells": "Total T cell count",
    "epi":          "Total epithelial tumor cell count",
    "mes":          "Total mesenchymal tumor cell count",
    "class1":       "Class I expressing tumor cells",
    "class1_class2":"Class I & II expressing tumor cells",
    "CD8":          "CD8 T cell count",
    "CD4":          "CD4 T cell count",
}

def make_scatter(x_col, y_col, title, out_path):
    fig, ax = plt.subplots(figsize=(7, 6))

    texts = []
    for patient, grp in data.groupby("patient"):
        ax.scatter(grp[x_col], grp[y_col],
                   color=pal[patient], s=60, zorder=3,
                   edgecolors="white", linewidths=0.5, label=patient)
        for _, row in grp.iterrows():
            roi_id = "ROI" + ROI_RE.match(row["roi"]).group(2)
            texts.append(ax.text(row[x_col], row[y_col], roi_id,
                                 fontsize=6.5, color=pal[patient]))

    adjust_text(texts, ax=ax, expand=(1.15, 1.3), force_points=(0.4, 0.6))
    ax.legend(title="Case", fontsize=7.5, title_fontsize=8,
              loc="upper left", framealpha=0.85, ncol=2)

    m_fit, b_fit = np.polyfit(data[x_col], data[y_col], 1)
    x_range = np.array([data[x_col].min(), data[x_col].max()])
    ax.plot(x_range, m_fit * x_range + b_fit,
            color="0.4", linewidth=1.2, linestyle="--", zorder=1)

    r = np.corrcoef(data[x_col], data[y_col])[0, 1]
    ax.text(0.97, 0.97, f"r = {r:.2f}", transform=ax.transAxes,
            ha="right", va="top", fontsize=9, color="0.4")

    ax.set_xlabel(LABELS.get(x_col, x_col), fontsize=10)
    ax.set_ylabel(LABELS.get(y_col, y_col), fontsize=10)
    ax.set_title(title, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    out_path.parent.mkdir(exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


# ── Generate both plots ───────────────────────────────────────────────────────
make_scatter(
    x_col    = "apCAF",
    y_col    = "tcells",
    title    = "apCAF vs T cell abundance per IMC ROI",
    out_path = Path(__file__).parent.parent / "figures/apcaf_vs_tcells_ic.png",
)

make_scatter(
    x_col    = "CAF",
    y_col    = "tcells",
    title    = "CAF vs T cell abundance per IMC ROI",
    out_path = Path(__file__).parent.parent / "figures/caf_vs_tcells_ic.png",
)

make_scatter(
    x_col    = "CAF",
    y_col    = "apCAF",
    title    = "CAF vs apCAF abundance per IMC ROI",
    out_path = Path(__file__).parent.parent / "figures/caf_vs_apcaf_ic.png",
)

def make_class1_scatter(out_path):
    fig, ax = plt.subplots(figsize=(7, 6))

    series = [
        ("epi_class1", "epi_class1_class2", "o", "epithelial"),
        ("mes_class1", "mes_class1_class2", "^", "mesenchymal"),
    ]

    texts = []
    for x_col, y_col, marker, phenotype in series:
        for patient, grp in data.groupby("patient"):
            ax.scatter(grp[x_col], grp[y_col],
                       color=pal[patient], s=60, marker=marker, zorder=3,
                       edgecolors="white", linewidths=0.5)
            for _, row in grp.iterrows():
                texts.append(ax.text(row[x_col], row[y_col], row["roi"],
                                     fontsize=6.5, color=pal[patient]))

    adjust_text(texts, ax=ax, expand=(1.15, 1.3), force_points=(0.4, 0.6))

    # best-fit lines per phenotype
    for x_col, y_col, _, phenotype in series:
        m_fit, b_fit = np.polyfit(data[x_col], data[y_col], 1)
        x_range = np.array([data[x_col].min(), data[x_col].max()])
        ax.plot(x_range, m_fit * x_range + b_fit,
                color="0.4", linewidth=1.2, linestyle="--", zorder=1)

    # marker legend for phenotype
    for marker, label in [("o", "epithelial"), ("^", "mesenchymal")]:
        ax.scatter([], [], marker=marker, color="0.4", s=50, label=label)
    ax.legend(fontsize=8, framealpha=0.85, loc="upper left")

    ax.set_xlabel("Class I expressing tumor cells", fontsize=10)
    ax.set_ylabel("Class I & II expressing tumor cells", fontsize=10)
    ax.set_title("Tumor antigen presentation: class I vs class I+II per IMC ROI", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    out_path.parent.mkdir(exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


make_scatter(
    x_col    = "epi",
    y_col    = "mes",
    title    = "Epithelial vs mesenchymal tumor cell abundance per IMC ROI",
    out_path = Path(__file__).parent.parent / "figures/epi_vs_mes_tumor_ic.png",
)

make_scatter(
    x_col    = "apCAF",
    y_col    = "CD8",
    title    = "apCAF vs CD8 T cell abundance per IMC ROI",
    out_path = Path(__file__).parent.parent / "figures/apcaf_vs_cd8_ic.png",
)

make_scatter(
    x_col    = "apCAF",
    y_col    = "CD4",
    title    = "apCAF vs CD4 T cell abundance per IMC ROI",
    out_path = Path(__file__).parent.parent / "figures/apcaf_vs_cd4_ic.png",
)

make_scatter(
    x_col    = "class1",
    y_col    = "class1_class2",
    title    = "Class I vs class I+II tumor cell abundance per IMC ROI",
    out_path = Path(__file__).parent.parent / "figures/class1_vs_class1class2_ic.png",
)
