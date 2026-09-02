#!/usr/bin/env python3
"""Stamp simulation, sample and patient into derived files that lack them.

Needed once, for files written before the extract script stamped identity: the
uns-writing step was lost when the series concatenation was switched to
anndata's own `concat_on_disk`, which does not carry uns through. Re-converting
to add three strings would cost eight minutes a simulation; this costs a second.

Safe to re-run -- it rewrites the values rather than appending, and reports what
it touched.
"""

import sys
from pathlib import Path

import h5py
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
DERIVED = BASE / "derived"
SAMPLE_MAP = BASE / "simulation_to_sample_id.csv"


def main():
    table = pd.read_csv(SAMPLE_MAP).set_index("simulation")
    stamped, already, locked, missing = 0, 0, [], []

    for path in sorted(DERIVED.glob("sim_*.h5ad")):
        sim = int(path.name.split("_")[1])
        if sim not in table.index:
            missing.append(path.name)
            continue
        row = table.loc[sim]
        values = {
            "simulation_id": str(sim),
            "sample_id": str(row["sample_id"]),
            "patient_id": str(row["patient_id"]),
        }
        try:
            handle = h5py.File(path, "a")
        except BlockingIOError:
            # Still being written by a running extract task. HDF5 takes a write
            # lock, so this is the file telling us to come back later rather than
            # a problem to work around.
            locked.append(path.name)
            continue
        with handle as h:
            uns = h.require_group("uns")
            if all(k in uns for k in values):
                already += 1
                continue
            uns.attrs["encoding-type"] = "dict"
            uns.attrs["encoding-version"] = "0.1.0"
            for key, value in values.items():
                if key in uns:
                    del uns[key]
                d = uns.create_dataset(key, data=value)
                d.attrs["encoding-type"] = "string"
                d.attrs["encoding-version"] = "0.2.0"
        stamped += 1
        print(f"  stamped {path.name}: {values['sample_id']} / {values['patient_id']}")

    print(f"\n{stamped} stamped, {already} already carried identity")
    if locked:
        print(f"{len(locked)} still being written, re-run when the batch finishes: "
              f"{', '.join(locked[:4])}{' ...' if len(locked) > 4 else ''}")
    if missing:
        print(f"no mapping for: {', '.join(missing)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
