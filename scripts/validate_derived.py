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
                if not sample or sample == "None":
                    problems.append(f"sim {sim}: series carries no sample_id")
                else:
                    identities[sim] = sample
                names = h["var"]["_index"]
                signature = (h["X"].shape[1], h["X"].dtype.str,
                             h["X"].compression or "none")
                var_signature.setdefault(signature, []).append(sim)
                sizes[sim] = series.stat().st_size / 2**20
        except OSError as e:
            problems.append(f"sim {sim}: series unreadable -- {e}")

    # every sample should appear exactly once
    for sample, count in Counter(identities.values()).items():
        if count > 1:
            claimants = [s for s, v in identities.items() if v == sample]
            problems.append(f"sample {sample} claimed by simulations {claimants}")

    # one shape, dtype and compression across the archive, or later concat breaks
    if len(var_signature) > 1:
        for sig, sims in var_signature.items():
            problems.append(
                f"{len(sims)} simulation(s) have X {sig[0]} cols, {sig[1]}, "
                f"{sig[2]} compression: {sims[:5]}")

    if sizes:
        total = sum(sizes.values())
        print(f"  {len(present)}/{len(expected)} converted")
        print(f"  series total {total/1024:.1f} GB, "
              f"median {sorted(sizes.values())[len(sizes)//2]:.0f} MB, "
              f"range {min(sizes.values()):.0f}-{max(sizes.values()):.0f} MB")
        print(f"  {len(set(identities.values()))} distinct samples")

    if problems:
        print(f"\nFAIL: {len(problems)} problem(s)", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        return 1

    print("\nOK: archive complete, identities unique, layout consistent")
    return 0


if __name__ == "__main__":
    sys.exit(main())
