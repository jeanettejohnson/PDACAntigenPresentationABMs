#!/usr/bin/env python3
"""Extract CD8 attack events from the derived series into a long-form table.

One row per (time, attacker, target). Nothing is aggregated here: what counts as
"an attack" is a modelling decision, and the alternatives -- duration-weighted
timepoints, distinct attacker/target engagements, first contacts, damage-weighted
totals -- are each one groupby away from this table but need different data if
baked in too early. The 88 GB of series is read once; every later question is
answered from a few MB of parquet.

    derived/attacks/<sim_id>-attacks.parquet    the events
    derived/attacks/<sim_id>-summary.parquet    denominators for normalising
    derived/attacks/all_attacks.parquet         every simulation, concatenated
    derived/attacks/all_summary.parquet

Attacks live only in the series. `total_attack_time` and
`attack_total_damage_delivered` are per-attacker cumulative totals, so a T cell
that attacked a class I and a class I+II target has one lumped number in the
snapshots with no way to split it. The breakdown by target exists only in the
per-timepoint `attack_target`, which names the cell being attacked right then.

Columns are read individually through h5py rather than by `read_h5ad`. The
series X is chunked (16384, 64), so five columns cost about a thirtieth of what
the whole matrix costs -- measured at 0.21 s against 6.9 s for one column of a
660k-row file.

Run in `physicell-analysis-260901`, which carries the parquet engine. It is not
in physicell-sim-260606 and is deliberately not added there: that environment's
figure stack is version-pinned, and adding anything re-solves it.

Writer and reader are the same environment, so parquet is readable where it is
written. `--csv` remains for handing the tables to something that has no parquet
engine at all -- a larger file with dtypes lost, so prefer parquet.
"""

import argparse
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
DERIVED = BASE / "derived"
OUT = DERIVED / "attacks"

#: Read from X. attack_target names the cell being attacked at that timepoint,
#: or -1 when the attacker is idle; the rest give effort and outcome.
X_COLUMNS = ("attack_target", "total_attack_time",
             "attack_total_damage_delivered", "damage", "dead")

ATTACKER = "CD8_Tcell"

IDENTITY = ("sim_type", "sim_db_id", "sample_id", "patient_id", "geometry", "sim_id")


def _decode(values):
    return [v.decode() if isinstance(v, bytes) else str(v) for v in values]


def _obs_column(handle, name):
    """One obs column as a numpy array of str, categorical or not."""
    node = handle["obs"][name]
    if isinstance(node, h5py.Group):          # categorical
        categories = np.asarray(_decode(node["categories"][:]), dtype=object)
        return categories[node["codes"][:]]
    return np.asarray(_decode(node[:]), dtype=object)


def split_cell_type(frame, column="target_cell_type"):
    """Lineage and antigen-presentation class from a tumour cell-type name.

    `epithelial_tumor_class1_class2` -> lineage epithelial, class class1_class2.
    Non-tumour targets keep the type and get no lineage or class, rather than
    being dropped: a CD8 attacking something unexpected is a finding, not noise.
    """
    kind = frame[column].astype(str)
    is_tumour = kind.str.contains("_tumor")
    lineage = kind.str.split("_tumor", n=1).str[0].where(is_tumour)
    rest = kind.str.split("_tumor", n=1).str[-1].str.lstrip("_")
    klass = rest.where(is_tumour & rest.ne(""), other=pd.NA)
    return lineage, klass


def extract(path):
    """Attack events and per-simulation totals for one series file."""
    with h5py.File(path, "r") as handle:
        var_names = _decode(handle["var"]["_index"][:])
        missing = [c for c in X_COLUMNS if c not in var_names]
        if missing:
            raise KeyError(f"{path.name} has no {', '.join(missing)}")
        values = {c: handle["X"][:, var_names.index(c)] for c in X_COLUMNS}

        obs_names = _decode(handle["obs"]["_index"][:])
        cell_type = _obs_column(handle, "cell_type")
        times = _obs_column(handle, "time")
        identity = {k: _obs_column(handle, k)[0] if k in handle["obs"] else ""
                    for k in IDENTITY}

    # obs_names are "<cell id>_<time>"; the id is what attack_target refers to.
    cell_id = np.array([n.rsplit("_", 1)[0] for n in obs_names], dtype=object)
    target = values["attack_target"].astype(np.int64)

    table = pd.DataFrame({
        "time": times,
        "cell_id": cell_id,
        "cell_type": cell_type,
        "attack_target": target,
        "total_attack_time": values["total_attack_time"],
        "attack_total_damage_delivered": values["attack_total_damage_delivered"],
        "damage": values["damage"],
        "dead": values["dead"],
    })

    # The target's type has to be read at the same timepoint: cells change type
    # and die, so a lookup that ignored time would attribute attacks to whatever
    # the target became later.
    by_time = table.set_index(["time", "cell_id"])
    attacks = table[(table.cell_type == ATTACKER) & (table.attack_target >= 0)].copy()
    key = pd.MultiIndex.from_arrays(
        [attacks["time"], attacks["attack_target"].astype(str)])
    attacks["target_cell_type"] = by_time["cell_type"].reindex(key).to_numpy()
    attacks["target_damage"] = by_time["damage"].reindex(key).to_numpy()
    attacks["target_dead"] = by_time["dead"].reindex(key).to_numpy()

    attacks = attacks.rename(columns={
        "cell_id": "attacker_id",
        "attack_target": "target_id",
        "total_attack_time": "attacker_total_attack_time",
        "attack_total_damage_delivered": "attacker_total_damage_delivered",
    }).drop(columns=["cell_type", "damage", "dead"])

    lineage, klass = split_cell_type(attacks)
    attacks["target_lineage"] = lineage
    attacks["target_class"] = klass
    for field, value in identity.items():
        attacks[field] = value

    # Denominators, so a later comparison across sets can normalise without
    # re-reading the series: how many attackers there were and how much tumour
    # was available to attack.
    final = table[table.time == table.time.iloc[-1]]
    counts = final.cell_type.value_counts()
    is_attacker = table.cell_type == ATTACKER
    summary = {
        **identity,
        "timepoints": int(table.time.nunique()),
        "rows": int(len(table)),
        "attack_rows": int(len(attacks)),
        # Rows are cell-timepoints, so a count of them is not a count of cells:
        # one T cell alive for 337 timepoints contributes 337 rows. Both are
        # useful denominators and they differ by three orders of magnitude, so
        # they are named for what they actually are.
        "attacker_rows": int(is_attacker.sum()),
        "attackers_distinct": int(table.loc[is_attacker, "cell_id"].nunique()),
        "attackers_final": int(counts.get(ATTACKER, 0)),
        "attackers_attacking_distinct": int(attacks["attacker_id"].nunique()),
        "unresolved_targets": int(attacks.target_cell_type.isna().sum()),
    }
    for cell_type_name, n in counts.items():
        if "tumor" in str(cell_type_name):
            summary[f"final_{cell_type_name}"] = int(n)

    return attacks, pd.DataFrame([summary])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--series", help="one series file, instead of all of them")
    parser.add_argument("--limit", type=int, help="stop after N simulations")
    parser.add_argument("--force", action="store_true",
                        help="re-extract simulations that already have output")
    parser.add_argument("--no-combine", action="store_true",
                        help="skip writing the concatenated tables")
    parser.add_argument("--csv", action="store_true",
                        help="also write the combined tables as CSV, for an "
                             "environment with no parquet engine")
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    paths = ([Path(args.series)] if args.series
             else sorted(DERIVED.glob("*-series.h5ad")))
    if args.limit:
        paths = paths[:args.limit]
    if not paths:
        print(f"no series files under {DERIVED}", file=sys.stderr)
        return 1

    print(f"{len(paths)} series to read")
    started, done, skipped, total_attacks = time.time(), 0, 0, 0
    for path in paths:
        stem = path.name[: -len("-series.h5ad")]
        events_out = OUT / f"{stem}-attacks.parquet"
        summary_out = OUT / f"{stem}-summary.parquet"
        if events_out.exists() and summary_out.exists() and not args.force:
            skipped += 1
            continue

        t0 = time.time()
        attacks, summary = extract(path)
        attacks.to_parquet(events_out, index=False)
        summary.to_parquet(summary_out, index=False)
        done += 1
        total_attacks += len(attacks)
        unresolved = int(summary["unresolved_targets"].iloc[0])
        note = f"  UNRESOLVED {unresolved}" if unresolved else ""
        print(f"  {stem:52s} {len(attacks):7,} attacks  {time.time()-t0:5.1f}s{note}")

    print(f"{done} extracted, {skipped} already present, "
          f"{total_attacks:,} attack rows, {time.time()-started:.0f}s")

    if not args.no_combine:
        for kind in ("attacks", "summary"):
            parts = sorted(OUT.glob(f"*-{kind}.parquet"))
            parts = [p for p in parts if not p.name.startswith("all_")]
            if not parts:
                continue
            combined = pd.concat([pd.read_parquet(p) for p in parts],
                                 ignore_index=True)
            out = OUT / f"all_{kind}.parquet"
            combined.to_parquet(out, index=False)
            print(f"  {out.name}: {len(combined):,} rows from {len(parts)} simulations")
            if args.csv:
                csv_out = out.with_suffix(".csv")
                combined.to_csv(csv_out, index=False)
                print(f"  {csv_out.name}: {csv_out.stat().st_size / 2**20:.1f} MB "
                      "(readable without a parquet engine)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
