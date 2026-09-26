#!/usr/bin/env python3
"""Check the whole derived archive, not one simulation at a time.

extract_simulation_results.py verifies each simulation as it writes it: files
present, row count against the raw .mat shapes, X layout, identity stamped.
Those are all within-simulation checks, and they cannot see the things that
only go wrong across simulations -- a missing run, two runs claiming the same
sample, var_names that differ between files and would break any later
concatenation, or a stray file that the combine step would pick up.

Run after the batch. Exit 0 means the archive is complete and self-consistent,
which is the claim the whole extract exists to support.
"""

import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path

import h5py

BASE = Path(__file__).resolve().parent.parent
DERIVED = BASE / "derived"
PCMM_DB = BASE / "data" / "pcmm.db"
KINDS = ("initial", "final", "series", "microenv")

#: What combine_snapshots.py writes next to the per-simulation files.
COMBINED = {"all_initial.h5ad", "all_final.h5ad"}

sys.path.insert(0, str(BASE))


def completed_ids():
    with sqlite3.connect(PCMM_DB) as con:
        return [r[0] for r in con.execute(
            "SELECT s.simulation_id FROM simulations s "
            "JOIN status_codes c USING(status_code_id) "
            "WHERE c.status_code = 'Completed' ORDER BY s.simulation_id")]


def sim_ids(ids):
    """simulation id -> sim_id stem, resolved the same way the extract does."""
    from scripts.resolve_samples import (caf_mhc2_rate, derived_stem, load_cohorts,
                                         output_dir, resolve, sim_type)
    kind = sim_type(BASE)
    mapping, _ = resolve(BASE, load_cohorts(BASE))
    return {s: derived_stem(kind, s, mapping[s], caf_mhc2_rate(output_dir(BASE, s)))
            for s in ids if s in mapping}


_STEMS = None


def path_for(sim, kind):
    ext = "h5" if kind == "microenv" else "h5ad"
    return DERIVED / f"{_STEMS.get(sim, f'unresolved-{sim}')}-{kind}.{ext}"


def _text(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def strays(expected_paths):
    """Anything in derived/ that no completed simulation accounts for.

    combine_snapshots.py globs derived/, so a file left from an earlier run --
    a renamed simulation, a removed one -- would be combined as if it belonged.
    Subfolders are downstream outputs and are left alone, except the extract's
    own scratch folders, which only a killed task leaves behind.
    """
    out = []
    for path in sorted(DERIVED.iterdir()):
        if path.is_dir():
            if path.name.startswith(".tmp_sim_"):
                out.append(f"{path.name}/ (scratch left by a killed extract task)")
        elif not path.name.startswith(".") and path.name not in COMBINED \
                and path not in expected_paths:
            out.append(path.name)
    return out


def main():
    global _STEMS
    # --only-present validates a partial archive: every cross-simulation check
    # still runs, but a simulation with no files is skipped instead of counted
    # as missing. For a subset run -- a smoke test, or a look at a batch still
    # in progress -- completeness is the one property that is legitimately not
    # true yet, and without this the other checks cannot be seen past it.
    only_present = "--only-present" in sys.argv
    from scripts.extract_simulation_results import x_layout_problems

    problems, expected = [], completed_ids()
    _STEMS = sim_ids(expected)
    if only_present:
        print(f"checking the {len(expected)} completed simulations against "
              f"{DERIVED} (subset mode: absent simulations are skipped, "
              f"completeness is NOT checked)")
    else:
        print(f"checking {len(expected)} completed simulations against {DERIVED}")

    present, identities, var_signature, sizes = [], {}, {}, {}
    var_names, layout = {}, defaultdict(list)
    for sim in expected:
        missing = [k for k in KINDS if not path_for(sim, k).exists()]
        if missing:
            # A simulation absent entirely is "not converted yet" in subset
            # mode; one with some files but not others is broken either way.
            if not (only_present and len(missing) == len(KINDS)):
                problems.append(f"sim {sim}: missing {', '.join(missing)}")
            continue
        present.append(sim)

        series = path_for(sim, "series")
        try:
            with h5py.File(series, "r") as h:
                uns = h.get("uns", {})
                values = {k: _text(uns[k][()]) if k in uns else None
                          for k in ("sample_id", "geometry", "caf_mhc2_rate")}
                sample, geometry, rate = (values[k] if values[k] not in ("", "None") else None
                                          for k in ("sample_id", "geometry", "caf_mhc2_rate"))
                if not sample:
                    problems.append(f"sim {sim}: series carries no sample_id")
                else:
                    identities[sim] = (sample, geometry, rate)
                # The names themselves, not just how many: anndata concatenates
                # on var_names, so two files with the same 459 columns in a
                # different order line up silently wrong. Hashed because the
                # full tuple is 459 strings per simulation and only ever gets
                # compared for equality.
                names = tuple(_text(n) for n in h["var"]["_index"][:])
                signature = (len(names), hash(names), h["X"].dtype.str)
                var_names.setdefault(signature, names)
                var_signature.setdefault(signature, []).append(sim)
                # Checked outright, not only for agreement: a cohort that was
                # uncompressed throughout would agree with itself perfectly.
                # Chunk messages carry each series' own shape; grouped as one.
                for p in x_layout_problems(h["X"]):
                    layout["X chunks are not X_CHUNKS capped to the series"
                           if p.startswith("X chunks") else p].append(sim)
                sizes[sim] = series.stat().st_size / 2**20
        except OSError as e:
            problems.append(f"sim {sim}: series unreadable -- {e}")

        try:
            with h5py.File(path_for(sim, "microenv"), "r") as h:
                rows = h["field"].shape[1]
                if "row_names" not in h:
                    layout["microenv has no row_names"].append(sim)
                elif len(h["row_names"]) != rows:
                    layout[f"microenv row_names do not match its {rows} rows"].append(sim)
                if len(h["time"]) != h["field"].shape[0]:
                    layout["microenv time does not match its timepoints"].append(sim)
        except (OSError, KeyError) as e:
            problems.append(f"sim {sim}: microenv unreadable -- {e}")

    for problem, sims in layout.items():
        problems.append(f"{len(sims)} simulation(s): {problem}: {sims[:5]}")

    # Every simulation should be a distinct point in the study design. That is
    # the sample alone for the three sets where one simulation is one sample,
    # but htan_geometries runs each of 73 samples through five arrangements
    # (c1-c5), so there the sample repeats by design and the pair is what must
    # be unique. Keep the c-prefix in the key: archives made before 2026-09-25
    # also ran c6, which had c2's exact layout, and without the prefix those
    # 438 runs collapse to 365. The CAF-contact MHC-II rate is part of the
    # point too: a sample run with the rule on and off is two points.
    for key, count in Counter(identities.values()).items():
        if count > 1:
            claimants = [s for s, v in identities.items() if v == key]
            sample, geometry, rate = key
            where = f"sample {sample}" + (f" geometry {geometry}" if geometry else "")
            where += f" caf_mhc2_rate {rate}" if rate else ""
            problems.append(f"{where} claimed by simulations {claimants}")

    # one shape and dtype across the archive, or later concat breaks
    if len(var_signature) > 1:
        majority = max(var_signature, key=lambda k: len(var_signature[k]))
        for sig, sims in var_signature.items():
            problems.append(
                f"{len(sims)} simulation(s) have X {sig[0]} cols, {sig[2]}: {sims[:5]}")
            if sig != majority:
                here, there = var_names[sig], var_names[majority]
                only_here = sorted(set(here) - set(there))
                only_there = sorted(set(there) - set(here))
                if only_here or only_there:
                    problems.append(
                        f"    columns only in {sims[:3]}: {only_here[:5]}; "
                        f"only in the majority: {only_there[:5]}")
                elif here != there:
                    problems.append(
                        "    same column names in a different order -- concat "
                        "would align on names, so check the writer, not the data")
                else:
                    problems.append("    same columns in the same order; the X "
                                    "dtype is what differs")

    expected_paths = {path_for(sim, k) for sim in expected for k in KINDS}
    for name in strays(expected_paths):
        problems.append(f"not from any completed simulation: {name}")

    if sizes:
        total = sum(sizes.values())
        print(f"  {len(present)}/{len(expected)} converted")
        print(f"  series total {total/1024:.1f} GB, "
              f"median {sorted(sizes.values())[len(sizes)//2]:.0f} MB, "
              f"range {min(sizes.values()):.0f}-{max(sizes.values()):.0f} MB")
        samples = {k[0] for k in identities.values()}
        print(f"  {len(samples)} distinct samples, "
              f"{len(set(identities.values()))} distinct sample/geometry/rate points")

    if problems:
        print(f"\nFAIL: {len(problems)} problem(s)", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        return 1

    if only_present:
        print(f"\nOK: the {len(present)} converted simulation(s) are "
              "self-consistent -- completeness not checked (--only-present)")
    else:
        print("\nOK: archive complete, identities unique, layout consistent")
    return 0


if __name__ == "__main__":
    sys.exit(main())
