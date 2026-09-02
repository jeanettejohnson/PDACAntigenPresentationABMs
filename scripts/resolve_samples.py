#!/usr/bin/env python3
"""Work out which sample each simulation was seeded from, from the run itself.

The model manager records how a simulation was configured, never which sample
the configuration came from. That link lived in `simulation_to_sample_id*.csv`
files generated alongside the runs -- which exist for three of the four
simulation sets, live in a checkout other than this one, and have to be
regenerated whenever simulations are added. htan_geometries has none at all.

Everything needed to reconstruct it is already in each clone, by two routes:

1. **Directly.** Where a run has its own initial-condition folder, the folder is
   named for the sample: imc_spatial's 48 simulations point at `ic_cells` rows
   called `JHH317ROI1` and so on. Nothing to infer.

2. **By composition.** Where every run shares one initial-condition folder and
   differs by variation -- htan_wellmixed, imc_wellmixed, htan_geometries --
   `ic_cell_variations.db` records the starting count of each cell type per
   variation. Those counts came from the assignment summary, one row per sample,
   so matching the count vector recovers the sample exactly. It is a lookup
   rather than a guess: the numbers are equal, not merely close.

Route 2 is checked against `simulation_to_sample_id.csv` where that file exists,
which is how it was validated in the first place.
"""

import argparse
import re
import sqlite3
import sys
from pathlib import Path

import pandas as pd

#: ic_cell_variations names cell types the way PhysiCell does; the assignment
#: summary names them the way pdacagviz's atlas palette does. Same alias map.
SIMULATION_TO_ATLAS = {
    "epithelial_tumor": "epithelial",
    "epithelial_tumor_class1": "epithelial_class1",
    "epithelial_tumor_class2": "epithelial_class2",
    "epithelial_tumor_class1_class2": "epithelial_class1_class2",
    "mesenchymal_tumor": "mesenchymal",
    "mesenchymal_tumor_class1": "mesenchymal_class1",
    "mesenchymal_tumor_class2": "mesenchymal_class2",
    "mesenchymal_tumor_class1_class2": "mesenchymal_class1_class2",
}


def _direct(base):
    """simulation -> sample, where each run has its own ic_cells folder."""
    with sqlite3.connect(base / "data" / "pcmm.db") as con:
        rows = con.execute(
            "SELECT s.simulation_id, ic.folder_name FROM simulations s "
            "JOIN ic_cells ic USING(ic_cell_id)"
        ).fetchall()
    names = {name for _, name in rows}
    # One shared folder means this route says nothing -- every run would get the
    # same name, which is the htan_wellmixed case that made the earlier version
    # of this plan wrong.
    if len(names) < 2:
        return {}
    return dict(rows)


def _variation_counts(base):
    """variation id -> {cell type: starting count}, from ic_cell_variations.db."""
    dbs = list((base / "data" / "inputs" / "ics" / "cells").glob("*/ic_cell_variations.db"))
    if not dbs:
        return {}
    frame = pd.read_sql("SELECT * FROM ic_cell_variations", sqlite3.connect(dbs[0]))
    # Only the count columns. htan_geometries' variations also carry inner_radius
    # and outer_radius per cell type -- 30 of its 44 columns -- and matching on
    # "name:" alone swallows those too, so radii get read as cell counts. That is
    # what made this set look like a synthetic sweep rather than 73 samples.
    columns = {c: re.search(r"name:([^/]+)/", c).group(1)
               for c in frame.columns if "name:" in c and c.endswith("/number")}
    if not columns:
        return {}
    counts = {}
    for row in frame.itertuples(index=False):
        record = frame.loc[frame.ic_cell_variation_id == row.ic_cell_variation_id]
        counts[int(row.ic_cell_variation_id)] = {
            name: int(record.iloc[0][col]) for col, name in columns.items()
        }
    return counts


def _by_composition(base, cohort):
    """simulation -> sample, by matching starting counts against the cohort."""
    variations = _variation_counts(base)
    if not variations or cohort is None:
        return {}

    # The run's own recipe decides which types matter. htan_wellmixed seeds 14,
    # htan_geometries 15, imc_wellmixed 8 -- reading it beats assuming it.
    # Variation names are PhysiCell's. The IMC summary uses those directly; the
    # HTAN one uses the atlas spellings. Try the names as they come, then aliased,
    # rather than picking one and being wrong for half the clones.
    raw = sorted(next(iter(variations.values())))
    for alias in (dict(), SIMULATION_TO_ATLAS):
        types = [alias.get(t, t) for t in raw]
        if all(t in cohort.columns for t in types):
            break
    else:
        return {}, [], sorted(variations)
    lookup = dict(zip(types, raw))

    # Exact count vectors, so a dict keyed on the tuple is the whole matcher.
    by_counts = {}
    for row in cohort.itertuples(index=False):
        key = tuple(int(getattr(row, t)) for t in types)
        by_counts.setdefault(key, []).append(row.sample_id)

    with sqlite3.connect(base / "data" / "pcmm.db") as con:
        sims = con.execute(
            "SELECT simulation_id, ic_cell_variation_id FROM simulations"
        ).fetchall()

    resolved, ambiguous, unmatched = {}, [], []
    for sim, variation in sims:
        counts = variations.get(int(variation))
        if counts is None:
            unmatched.append(sim)
            continue
        candidates = by_counts.get(tuple(counts[lookup[t]] for t in types), [])
        if len(candidates) == 1:
            resolved[sim] = candidates[0]
        elif candidates:
            ambiguous.append((sim, candidates))
        else:
            unmatched.append(sim)
    return resolved, ambiguous, unmatched


def _from_csv(base):
    """simulation -> sample, from a checked-in mapping if one is present.

    Last resort, not first: a CSV is generated once alongside a run and then
    drifts, whereas the database and the variation counts are what the run
    actually used. But imc_wellmixed's source composition is not in any clone --
    the IMC summary here is a three-sample test fixture against 48 runs -- so
    without this route those simulations have no identity at all.
    """
    for pattern in ("simulation_to_sample_id*.csv",):
        for path in sorted(Path(base).glob(pattern)):
            table = pd.read_csv(path)
            if {"simulation", "sample_id"} <= set(table.columns):
                return dict(zip(table["simulation"], table["sample_id"])), path.name
    return {}, None


def resolve(base, cohorts=None):
    """simulation id -> sample id, by whichever route the clone supports.

    Returns ``(mapping, report)``. The report names what could not be resolved
    rather than leaving it to be discovered downstream.
    """
    base = Path(base)
    direct = _direct(base)
    if direct:
        return direct, {"route": "ic_cells folder", "resolved": len(direct),
                        "ambiguous": [], "unmatched": []}

    best = ({}, [], [], None)
    for label, cohort in (cohorts or []):
        resolved, ambiguous, unmatched = _by_composition(base, cohort)
        if len(resolved) > len(best[0]):
            best = (resolved, ambiguous, unmatched, label)
    resolved, ambiguous, unmatched, label = best
    if resolved:
        return resolved, {"route": f"initial composition ({label})",
                          "resolved": len(resolved),
                          "ambiguous": ambiguous, "unmatched": unmatched}

    from_csv, name = _from_csv(base)
    if from_csv:
        return from_csv, {"route": f"csv ({name})", "resolved": len(from_csv),
                          "ambiguous": [], "unmatched": []}

    # Nothing resolved, and that may be correct rather than a failure:
    # htan_geometries sweeps synthetic cell-count conditions, so its runs are
    # identified by (condition, geometry) and were never seeded from a sample.
    return {}, {"route": "none -- runs may not be sample-derived",
                "resolved": 0, "ambiguous": [], "unmatched": []}


#: Compartments the geometry annuli are defined over. Confirmed against the
#: imm/caf/tum fields in run_htan_geometries.jl; everything unlisted is tumour.
COMPARTMENT = {
    "CAF": "caf", "apCAF": "caf",
    "CD4_Tcell": "immune", "CD8_Tcell": "immune",
    "Treg": "immune", "CD8_exhausted": "immune",
}


def sim_type(base):
    """Which of the four simulation sets this clone holds.

    Decided from what the runs were configured with, never from the directory
    name -- the same conversion runs in four checkouts and derived files get
    copied, so a path is not identity.
    """
    base = Path(base)
    with sqlite3.connect(base / "data" / "pcmm.db") as con:
        sims = pd.read_sql("SELECT * FROM simulations", con)
    if sims.ic_cell_id.nunique() > 1:
        return "imc_spatial"            # per-run initial conditions

    for db in (base / "data" / "inputs" / "configs").glob("*/config_variations.db"):
        cfg = pd.read_sql("SELECT * FROM config_variations", sqlite3.connect(db))
        if any("spatial_config_index" in c for c in cfg.columns):
            return "htan_geometries"    # annulus sweep

    seeded = set()
    for db in (base / "data" / "inputs" / "ics" / "cells").glob("*/ic_cell_variations.db"):
        table = pd.read_sql("SELECT * FROM ic_cell_variations", sqlite3.connect(db))
        seeded = {re.search(r"name:([^/]+)/", c).group(1)
                  for c in table.columns if c.endswith("/number")}
        break
    return "htan_wellmixed" if "Treg" in seeded else "imc_wellmixed"


def geometries(base):
    """variation id -> annulus layout, read from the radii rather than a label.

    Compartments run centre outward: "-" joins those sharing a ring, "_"
    separates rings. Derived from the data, so a new arrangement would be
    described correctly rather than mislabelled.
    """
    base = Path(base)
    for db in (base / "data" / "inputs" / "ics" / "cells").glob("*/ic_cell_variations.db"):
        table = pd.read_sql("SELECT * FROM ic_cell_variations",
                            sqlite3.connect(db)).set_index("ic_cell_variation_id")
        radii = {}
        for col in table.columns:
            if "radius" not in col:
                continue
            name = re.search(r"name:([^/]+)/", col).group(1)
            radii.setdefault(name, {})["in" if "inner" in col else "out"] = col
        if not radii:
            return {}
        out = {}
        for vid, row in table.iterrows():
            rings = {}
            for name, cols in radii.items():
                key = (row[cols["in"]], row[cols["out"]])
                rings.setdefault(key, set()).add(COMPARTMENT.get(name, "tumor"))
            out[int(vid)] = "_".join("-".join(sorted(c))
                                     for _, c in sorted(rings.items()))
        return out
    return {}


def load_cohorts(base):
    """Every assignment summary this clone carries, widest first.

    Both summaries sit in every clone, and which one a given simulation set was
    seeded from is not recorded anywhere -- so rather than guess from the folder
    name, the matcher tries each and keeps whichever resolves. A wrong guess
    cannot silently succeed: the counts are matched exactly, so a summary that
    did not seed these runs matches nothing at all.
    """
    base = Path(base)
    sys.path.insert(0, str(base))
    cohorts = []
    try:
        from scripts.htan_cohort import IMC_CSV, load, load_imc
        if (base / "assignmentsummary_HTAN_singlecell.csv").is_file():
            cohorts.append(("HTAN", load(long=False)))
        if IMC_CSV.is_file():
            cohorts.append(("IMC", load_imc(long=False)))
    except (ImportError, KeyError, FileNotFoundError):
        pass
    return cohorts


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo", default=Path(__file__).resolve().parent.parent)
    ap.add_argument("--check", metavar="CSV",
                    help="compare against a simulation_to_sample_id csv")
    args = ap.parse_args()

    base = Path(args.repo)
    mapping, report = resolve(base, load_cohorts(base))
    print(f"{base.parent.name}: resolved {report['resolved']} via {report['route']}")
    for sim, candidates in report["ambiguous"][:5]:
        print(f"  ambiguous: simulation {sim} matches {candidates}")
    if report["unmatched"]:
        print(f"  unmatched: {report['unmatched'][:8]}")

    if args.check:
        truth = pd.read_csv(args.check).set_index("simulation")["sample_id"].to_dict()
        shared = set(mapping) & set(truth)
        wrong = {s: (mapping[s], truth[s]) for s in shared if mapping[s] != truth[s]}
        print(f"  against {Path(args.check).name}: {len(shared)} shared, "
              f"{len(wrong)} disagree")
        for sim, (got, want) in list(wrong.items())[:5]:
            print(f"    simulation {sim}: resolved {got}, csv says {want}")
        return 1 if wrong else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
