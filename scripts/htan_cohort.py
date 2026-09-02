#!/usr/bin/env python3
"""Read assignmentsummary_HTAN_singlecell.csv in the model's own vocabulary.

The CSV names cell types the way the HTAN annotation does -- Pattern2, Pattern7,
CD8_T_cytotoxic -- while every palette, figure and simulation output uses the
PhysiCell names. The mapping between them is not a matter of taste: it is the one
`slurm/archive/run_htan_singlecell_tme_geometries.jl` applied when it built the
initial conditions, so using anything else would label a figure differently from
the simulation it describes.

Pattern2 is the epithelial lineage and Pattern7 the mesenchymal one; cells
carrying both patterns count as epithelial, which is what the Julia does.
"""

from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent.parent
CSV = BASE / "assignmentsummary_HTAN_singlecell.csv"

#: The IMC summary already names its columns the way PhysiCell does, so it needs
#: no mapping -- only the HTAN one speaks Pattern2/Pattern7.
#:
#: Note the path. `assignmentsummary_JHH_IMC_test.csv` sits at the repository root
#: in every clone and holds three samples: it is a fixture, not the cohort. The
#: real one has 48 rows and lives under prep_imc_spatial/, which is what
#: slurm/run_imc_wellmixed.jl reads when it seeds the runs.
IMC_CSV = BASE / "prep_imc_spatial" / "assignmentsummary_JHH_IMC.csv"

#: PhysiCell cell type -> the CSV columns summed to produce it.
#: Transcribed from run_htan_singlecell_tme_geometries.jl lines 42-59.
COMPOSITION = {
    "CAF": ["CAF"],
    "apCAF": ["apCAF"],
    "CD4_Tcell": ["CD4_T"],
    "CD8_Tcell": ["CD8_T", "CD8_T_cytotoxic", "Proliferating_T"],
    "Treg": ["Treg"],
    "epithelial": ["Pattern2_Pattern7", "Pattern2"],
    "epithelial_class1": ["Pattern2_Pattern7_class_1", "Pattern2_class_1"],
    "epithelial_class1_class2": ["Pattern2_Pattern7_class_1_class_2",
                                 "Pattern2_class_1_class_2"],
    "epithelial_class2": ["Pattern2_Pattern7_class_2", "Pattern2_class_2"],
    "mesenchymal": ["Pattern7"],
    "mesenchymal_class1": ["Pattern7_class_1"],
    "mesenchymal_class1_class2": ["Pattern7_class_1_class_2"],
    "mesenchymal_class2": ["Pattern7_class_2"],
    "PDAC_unclassified": ["PDAC_unclassified"],
    # Seeded by htan_geometries but not by htan_wellmixed. Harmless to compute
    # either way -- the matcher asks only for the types a run actually seeded.
    "CD8_exhausted": ["Exhausted_T"],
}

#: Annotated but never seeded by any of the four simulation sets.
NOT_SEEDED = ["other_immune", "other_tissue", "PDAC"]


def load_imc(long=True):
    """Cell counts per IMC sample. Columns are already PhysiCell names."""
    raw = pd.read_csv(IMC_CSV)
    raw.columns = [c.strip() for c in raw.columns]
    keep = [c for c in raw.columns if c not in ("sample_id", "total")]
    wide = raw[["sample_id"] + keep].copy()
    wide["patient_id"] = wide["sample_id"].str.extract(r"^([A-Za-z]+\d+)")[0]
    if not long:
        return wide
    return wide.melt(id_vars=["sample_id", "patient_id"],
                     var_name="cell_type", value_name="count")


def load(long=True):
    """Cell counts per sample, in PhysiCell names.

    Returns long form by default -- one row per sample and cell type, which is
    what every pdacagviz function expects. Pass ``long=False`` for the wide table.
    """
    raw = pd.read_csv(CSV)
    raw.columns = [c.strip() for c in raw.columns]

    wide = pd.DataFrame({"sample_id": raw["sample_id"]})
    for cell_type, columns in COMPOSITION.items():
        missing = [c for c in columns if c not in raw.columns]
        if missing:
            raise KeyError(f"{cell_type}: {CSV.name} has no column(s) {missing}")
        wide[cell_type] = raw[columns].sum(axis=1).astype(int)

    wide["patient_id"] = wide["sample_id"].str.split("P").str[0]
    if not long:
        return wide

    return wide.melt(
        id_vars=["sample_id", "patient_id"],
        var_name="cell_type",
        value_name="count",
    )


if __name__ == "__main__":
    df = load()
    print(f"{df['sample_id'].nunique()} samples, {df['patient_id'].nunique()} patients, "
          f"{df['cell_type'].nunique()} cell types, {df['count'].sum():,} cells")
