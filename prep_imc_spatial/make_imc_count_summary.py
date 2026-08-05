"""
make_imc_count_summary.py
-------------------------
Build a wide per-ROI cell-count table from the JHH IMC initial-condition files --
the IMC analog of assignmentsummary_HTAN_singlecell.csv.

Output: assignmentsummary_JHH_IMC.csv, one row per ROI (sample_id), one column
per PhysiCell cell type, values = cell counts.

Source files live in PhysiCell/config/ics/JHH_IMC/:
  - *.csv : 'type' column is already in PhysiCell type space -> counted directly
  - *.txt : raw QuPath export; 'Classification' mapped via TYPE_MAP, then counted

One canonical file is chosen per ROI (prefer _withTreg.csv > plain .csv > .txt).
other_tissue, duct_filler and blank/unmapped cells are dropped (no annulus patch
in the well-mixed cells.xml / not simulatable).
"""

import re
import sys
from pathlib import Path

import pandas as pd

# QuPath Classification -> PhysiCell type. Copied verbatim from
# assemble_initial_conditions.py so the .txt files map identically to the .csv files.
TYPE_MAP = {
    "CAF: Other": "CAF",
    "CAF: HLA-DR": "apCAF",
    "ductal: Other": "epithelial_normal",
    "ductal: HLA-DR": "epithelial_normal",
    "tumor_epithelial: Other": "epithelial_tumor_class1",
    "tumor_epithelial: HLA-DR": "epithelial_tumor_class1_class2",
    "tumor_mesenchymal: Other": "mesenchymal_tumor_class1",
    "tumor_mesenchymal: HLA-DR": "mesenchymal_tumor_class1_class2",
    "CD4 T cell: Other": "CD4_Tcell",
    "CD4 T cell: HLA-DR": "CD4_Tcell",
    "CD4 T cell: FOXP3": "Treg",
    "CD4 T cell: FOXP3: HLA-DR": "Treg",
    "CD8 T cell: Other": "CD8_Tcell",
    "CD8 T cell: HLA-DR": "CD8_Tcell",
    "Myeloid: Other": "macrophage",
    "Myeloid: HLA-DR": "macrophage",
    "B cell: Other": "B cell",
    "B cell: HLA-DR": "B cell",
    "CD57: Other": "CD8_Tcell",
    "CD57: HLA-DR": "CD8_Tcell",
    "Other: HLA-DR": "other_tissue",
    "Other": "other_tissue",
}

DROP = {"other_tissue", "duct_filler"}  # present in IMC, not placeable well-mixed
ROI_RE = re.compile(r"^(JHH\d+[A-Z]?ROI\d+)")  # canonical ROI key, e.g. JHH317ROI1 or JHH417RROI4

# Inputs resolve against the repo root; the summary is written next to this
# script, so the prep pipeline's outputs live with the code that makes them.
HERE = Path(__file__).parent
REPO = HERE.parent
IMC_DIR = REPO / "PhysiCell/user_projects/antigen_presentation/config/ics/JHH_IMC"
OUT_PATH = HERE / "assignmentsummary_JHH_IMC.csv"


def choose_canonical(files):
    """One file per ROI: prefer {roi}_withTreg.csv > {roi}.csv > {roi}.txt."""
    by_roi = {}
    for f in files:
        m = ROI_RE.match(f.stem)
        if m:
            by_roi.setdefault(m.group(1), {})[f.name] = f
    chosen = {}
    for roi, names in sorted(by_roi.items()):
        pick = (
            names.get(f"{roi}_withTreg.csv")
            or names.get(f"{roi}.csv")
            or names.get(f"{roi}.txt")
        )
        if pick is not None:
            chosen[roi] = pick
    return chosen


def types_for(path):
    """Series of PhysiCell type labels for one IC file."""
    if path.suffix.lower() == ".txt":
        cls = pd.read_csv(path, sep="\t", usecols=["Classification"])["Classification"]
        cls = cls.astype(str).str.strip()
        unmapped = sorted(set(cls) - set(TYPE_MAP))
        if unmapped:
            print(
                f"  ! {path.name}: unmapped classifications dropped: {unmapped}",
                file=sys.stderr,
            )
        return cls.map(TYPE_MAP)
    return pd.read_csv(path, usecols=["type"])["type"]


def main():
    if not IMC_DIR.is_dir():
        raise SystemExit(f"IMC directory not found: {IMC_DIR}")

    chosen = choose_canonical(sorted(IMC_DIR.glob("JHH*")))
    if not chosen:
        raise SystemExit(f"No IMC ROI files found in {IMC_DIR}")

    rows = {}
    for roi, path in chosen.items():
        t = types_for(path).dropna().astype(str).str.strip()
        t = t[(t != "") & (~t.isin(DROP))]
        rows[roi] = t.value_counts()
        print(f"{roi:14s} <- {path.name:42s} {int(t.size):6d} cells")

    summary = pd.DataFrame(rows).T.fillna(0).astype(int).sort_index().sort_index(axis=1)
    summary.insert(0, "total", summary.sum(axis=1))
    summary.index.name = "sample_id"
    summary.reset_index().to_csv(OUT_PATH, index=False)
    print(
        f"\nWrote {OUT_PATH}  ({summary.shape[0]} ROIs x {summary.shape[1] - 1} type cols)"
    )


if __name__ == "__main__":
    main()
