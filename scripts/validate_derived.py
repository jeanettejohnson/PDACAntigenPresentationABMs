#!/usr/bin/env python3
"""Check the whole derived archive, not one simulation at a time.

extract_simulation_results.py verifies each simulation as it writes it: files
present, row count against the raw .mat shapes, size, identity stamped. Those
are all within-simulation checks, and they cannot see the things that only go
wrong across simulations -- a missing run, two runs claiming the same sample,
or var_names that differ between files and would break any later concatenation.

Run after the batch. Exit 0 means the archive is complete and self-consistent,
which is the claim the whole extract exists to support.
"""

import sqlite3
import sys
from collections import Counter
from pathlib import Path

import h5py
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
DERIVED = BASE / "derived"
PCMM_DB = BASE / "data" / "pcmm.db"
KINDS = ("initial", "final", "series", "microenv")


def completed_ids():
    with sqlite3.connect(PCMM_DB) as con:
        return [r[0] for r in con.execute(
            "SELECT s.simulation_id FROM simulations s "
            "JOIN status_codes c USING(status_code_id) "
            "WHERE c.status_code = 'Completed' ORDER BY s.simulation_id")]


def sim_ids():
    """simulation id -> sim_id stem, resolved the same way the extract does."""
    import sys
    sys.path.insert(0, str(BASE))
    from scripts.resolve_samples import load_cohorts, resolve, sim_type
    kind = sim_type(BASE)
    mapping, _ = resolve(BASE, load_cohorts(BASE))
    return {s: f"{kind}-{s:03d}-{n}" for s, n in mapping.items()}


_STEMS = None


def path_for(sim, kind):
    global _STEMS
    if _STEMS is None:
        _STEMS = sim_ids()
    ext = "h5" if kind == "microenv" else "h5ad"
    return DERIVED / f"{_STEMS.get(sim, f'unresolved-{sim}')}-{kind}.{ext}"


def main():
    problems, expected = [], completed_ids()
    print(f"checking {len(expected)} completed simulations against {DERIVED}")

    present, identities, var_signature, sizes = [], {}, {}, {}
    var_names = {}
    for sim in expected:
        missing = [k for k in KINDS if not path_for(sim, k).exists()]
        if missing:
            problems.append(f"sim {sim}: missing {', '.join(missing)}")
            continue
        present.append(sim)

        series = path_for(sim, "series")
        try:
            with h5py.File(series, "r") as h:
                uns = h.get("uns", {})
                sample = uns["sample_id"][()] if "sample_id" in uns else None
                if isinstance(sample, bytes):
                    sample = sample.decode()
                geometry = uns["geometry"][()] if "geometry" in uns else None
                if isinstance(geometry, bytes):
                    geometry = geometry.decode()
                if not geometry or geometry == "None":
                    geometry = None
                if not sample or sample == "None":
                    problems.append(f"sim {sim}: series carries no sample_id")
                else:
                    identities[sim] = (sample, geometry)
                # The names themselves, not just how many: anndata concatenates
                # on var_names, so two files with the same 459 columns in a
                # different order line up silently wrong. Hashed because the
                # full tuple is 459 strings per simulation and only ever gets
                # compared for equality.
                names = tuple(n.decode() if isinstance(n, bytes) else str(n)
                              for n in h["var"]["_index"][:])
                signature = (len(names), hash(names), h["X"].dtype.str,
                             h["X"].compression or "none")
                var_names.setdefault(signature, names)
                var_signature.setdefault(signature, []).append(sim)
                sizes[sim] = series.stat().st_size / 2**20
        except OSError as e:
            problems.append(f"sim {sim}: series unreadable -- {e}")

    # Every simulation should be a distinct point in the study design. That is
    # the sample alone for the three sets where one simulation is one sample,
    # but htan_geometries runs each of 73 samples through six arrangements, so
    # there the sample repeats by design and the pair is what must be unique.
    # The geometry's c-prefix is load-bearing: c2 and c6 are the same layout,
    # and keying on the descriptive part alone collapses 438 runs to 365.
    for key, count in Counter(identities.values()).items():
        if count > 1:
            claimants = [s for s, v in identities.items() if v == key]
            sample, geometry = key
            where = f"sample {sample}" + (f" geometry {geometry}" if geometry else "")
            problems.append(f"{where} claimed by simulations {claimants}")

    # one shape, dtype and compression across the archive, or later concat breaks
    if len(var_signature) > 1:
        majority = max(var_signature, key=lambda k: len(var_signature[k]))
        for sig, sims in var_signature.items():
            problems.append(
                f"{len(sims)} simulation(s) have X {sig[0]} cols, {sig[2]}, "
                f"{sig[3]} compression: {sims[:5]}")
            if sig != majority:
                here, there = set(var_names[sig]), set(var_names[majority])
                only_here, only_there = sorted(here - there), sorted(there - here)
                if only_here or only_there:
                    problems.append(
                        f"    columns only in {sims[:3]}: {only_here[:5]}; "
                        f"only in the majority: {only_there[:5]}")
                else:
                    problems.append(
                        "    same column names in a different order -- concat "
                        "would align on names, so check the writer, not the data")

    if sizes:
        total = sum(sizes.values())
        print(f"  {len(present)}/{len(expected)} converted")
        print(f"  series total {total/1024:.1f} GB, "
              f"median {sorted(sizes.values())[len(sizes)//2]:.0f} MB, "
              f"range {min(sizes.values()):.0f}-{max(sizes.values()):.0f} MB")
        samples = {k[0] for k in identities.values()}
        print(f"  {len(samples)} distinct samples, "
              f"{len(set(identities.values()))} distinct sample/geometry pairs")

    if problems:
        print(f"\nFAIL: {len(problems)} problem(s)", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        return 1

    print("\nOK: archive complete, identities unique, layout consistent")
    return 0


if __name__ == "__main__":
    sys.exit(main())
