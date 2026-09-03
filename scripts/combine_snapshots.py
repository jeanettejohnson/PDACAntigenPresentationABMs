#!/usr/bin/env python3
"""Combine the per-simulation snapshots into one initial and one final AnnData.

The two snapshots are small enough to concatenate in memory -- 149k and 294k
cells, 0.3 and 0.5 GB raw -- so they get built with plain `anndata.concat` and
can be opened normally.

The series is deliberately not combined here. At 69.2M rows it is 118 GB
uncompressed, which `concat_on_disk` would have to write in full before the
compression pass, and the result could never be read outside `backed="r"`. The
same object costs one line at read time and no duplicate archive:

    import anndata as ad, glob
    series = ad.concat([ad.read_h5ad(f, backed="r")
                        for f in sorted(glob.glob("derived/*-series.h5ad"))])

Identity moves from uns into obs here, because it varies per row once several
simulations share a file.

Run in `physicell-analysis-260901`, the same environment that wrote the
per-simulation files. That is the point rather than a convenience: string columns
written by one anndata version can be nullable-string arrays another version reads
happily but refuses to write back without an opt-in flag. Keeping one writer across
the archive avoids the question entirely.
"""

import sys
from pathlib import Path

import anndata as ad
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
DERIVED = BASE / "derived"


def combine(kind):
    parts, labels = [], []
    # Stems are <sim_type>-<nnn>-<sample>; the numeric field is the second, and
    # no field value contains a dash, so this split is unambiguous.
    for path in sorted(DERIVED.glob(f"*-{kind}.h5ad"),
                       key=lambda p: int(p.name.split("-")[1])):
        sim = int(path.name.split("-")[1])
        a = ad.read_h5ad(path)
        # The extract already writes identity into obs, so nothing to reconstruct
        # from uns here -- but older files may predate that, hence the fallback.
        for key in ("sim_type", "sim_db_id", "sample_id", "patient_id", "geometry", "sim_id"):
            if key not in a.obs.columns:
                a.obs[key] = str(a.uns.get(key, ""))
        # obs_names are <cell>_<time>: unique within a simulation, not across.
        a.obs_names = [f"{sim}-{name}" for name in a.obs_names]
        parts.append(a)
        labels.append(sim)

    if not parts:
        raise FileNotFoundError(f"no *-{kind}.h5ad under {DERIVED}")

    combined = ad.concat(parts, axis=0, join="outer", merge="unique",
                         pairwise=True)
    for key in ("sim_type", "sim_db_id", "sample_id", "patient_id", "geometry", "sim_id"):
        if key in combined.obs.columns:
            combined.obs[key] = combined.obs[key].astype("category")
    combined.uns["source"] = f"{len(parts)} per-simulation {kind} snapshots"
    combined.uns["simulations"] = str(min(labels)) + "-" + str(max(labels))
    return combined


def main():
    for kind in ("initial", "final"):
        out = DERIVED / f"all_{kind}.h5ad"
        combined = combine(kind)
        combined.write_h5ad(out, compression="gzip")
        size = out.stat().st_size / 2**20
        print(f"{out.name}: {combined.n_obs:,} cells x {combined.n_vars} vars, "
              f"{combined.obs['sample_id'].nunique()} samples, "
              f"{len(combined.obsp)} graph matrices, {size:.0f} MB")
        if combined.obs_names.duplicated().any():
            print("  WARNING: duplicate obs_names", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
