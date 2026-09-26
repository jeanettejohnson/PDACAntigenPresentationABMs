#!/usr/bin/env python3
"""Convert PhysiCell output into a compact, lossless archive under derived/.

Four files per simulation:

    derived/<sim_id>-initial.h5ad   t=0 snapshot, with contact graphs in obsp
    derived/<sim_id>-final.h5ad     terminal snapshot, with contact graphs
    derived/<sim_id>-series.h5ad    every timepoint, every cell column, no graphs
    derived/<sim_id>-microenv.h5    the substrate fields on their voxel grid

where <sim_id> is <sim_type>-<nnn>-<sample_id>, e.g.
htan_wellmixed-001-HT056P1_S1PA or imc_spatial-003-JHH317ROI3, with -cafmhc2
appended for a run that had the CAF-contact MHC-II rule on.

A small fraction of the raw output, and nothing a figure or a spatial analysis
is likely to want is dropped. The point is to read the raw output once:
everything downstream reads derived/.

The four files are built in a scratch folder and moved into derived/ only once
every check has passed, so derived/ never holds a set that failed one -- a
rerun would otherwise take it for finished and skip it.

Run in the `physicell-analysis-260901` environment, never in physicell-sim-260606
-- see conda_env_configs/physicell_analysis_260901.yaml for why they are kept
apart. That environment covers conversion, figures and spatial analysis; the
simulation environment is left alone so its pinned figure stack never re-solves.
"""

import argparse
import re
import shutil
import sqlite3
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io

BASE = Path(__file__).resolve().parent.parent
OUTPUTS = BASE / "data" / "outputs" / "simulations"
DERIVED = BASE / "derived"
PCMM_DB = BASE / "data" / "pcmm.db"

#: X chunk shape. Measured on the 2-hour-save series (221k and 702k rows): one
#: column reads in 0.020s / 0.053s against 0.063s / 0.153s for the 16384x64 it
#: replaces, and five neighbouring columns 3x faster too, with no read pattern
#: slower -- a timestep slab and the full matrix cost the same. 1 MiB per chunk,
#: so it fits h5py's default chunk cache; narrower chunks read columns faster
#: still but slow down row reads. Capped to the series by `x_chunks`.
X_CHUNKS = (16384, 16)

#: pcdl scales every numeric column by default. For an archive that is wrong:
#: a stored "cell volume" would be a unitless 0-1 quantity with nothing saying so.
SCALE = None

#: Keep every column, including the ones constant within a simulation. 397 of the
#: 459 are constant, and gzip reduces them to almost nothing -- dropping them saves
#: 13% and costs immunogenicities, attack_rates and every other vector variable.
#: values=2 would also evaluate constancy per file, giving different columns per
#: simulation and breaking concat.
VALUES = 1

#: The first four rows of a microenvironment .mat, before one row per substrate.
VOXEL_ROWS = (("x", "micron"), ("y", "micron"), ("z", "micron"),
              ("volume", "cubic micron"))


#: Identity, resolved once per run rather than per simulation. scripts/
#: resolve_samples.py works it out from the model manager and the input variation
#: databases -- see its docstring for why a CSV is the last resort rather than
#: the first.
_IDENTITY = None


def identity():
    """(sim_type, {simulation_id: sample_id}, {variation_id: geometry})."""
    global _IDENTITY
    if _IDENTITY is None:
        sys.path.insert(0, str(BASE))
        from scripts.resolve_samples import (geometries, load_cohorts, resolve,
                                             sim_type)

        kind = sim_type(BASE)
        mapping, _report = resolve(BASE, load_cohorts(BASE))
        geo = geometries(BASE) if kind == "htan_geometries" else {}
        _IDENTITY = (kind, mapping, geo)
    return _IDENTITY


#: Patient is the leading letters+digits of the sample: HT056P1_S1PA -> HT056,
#: JHH317ROI1 -> JHH317. Jeanette's build_sample_mapping_v2.py splits on "P1_"
#: instead, which agrees on HTAN but has nothing to split on for IMC; the regex
#: covers both cohorts and resolves all 242 sample ids with none left over.
PATIENT_RE = re.compile(r"^([A-Za-z]+\d+)")


def describe(simulation_id):
    """Every identity field for one simulation.

    sim_id joins its fields with "-" while the fields themselves join words with
    "_", so the stem splits unambiguously -- no field value contains a dash.
    """
    kind, mapping, geo = identity()
    # identity() put BASE on sys.path
    from scripts.resolve_samples import caf_mhc2_rate, derived_stem

    sample = mapping.get(simulation_id)
    geometry = None
    if geo:
        with sqlite3.connect(PCMM_DB) as con:
            row = con.execute(
                "SELECT ic_cell_variation_id, config_variation_id "
                "FROM simulations WHERE simulation_id = ?", (simulation_id,)
            ).fetchone()
        if row:
            geometry = f"c{row[1]}_{geo.get(int(row[0]), 'unknown')}"
    match = PATIENT_RE.match(str(sample)) if sample else None
    rate = caf_mhc2_rate(OUTPUTS / str(simulation_id) / "output")
    return {
        "sim_type": kind,
        "sim_db_id": int(simulation_id),
        "sample_id": sample,
        "patient_id": match.group(1) if match else sample,
        "geometry": geometry,
        "caf_mhc2_rate": rate,
        "sim_id": derived_stem(kind, simulation_id, sample, rate),
    }


def completed_simulations():
    """Simulation ids the manager marks Completed, in order."""
    with sqlite3.connect(PCMM_DB) as con:
        rows = con.execute(
            "SELECT s.simulation_id FROM simulations s "
            "JOIN status_codes c USING(status_code_id) "
            "WHERE c.status_code = 'Completed' ORDER BY s.simulation_id"
        ).fetchall()
    return [r[0] for r in rows]


def x_chunks(shape):
    """The X chunk shape for a series of `shape`: X_CHUNKS, capped at the data.

    h5repack does not clip a chunk that is larger than the dataset. It drops the
    chunked layout, and the gzip filter with it, and exits 0 -- which is how an
    8,888-row series came out uncompressed.
    """
    return tuple(max(1, min(c, n)) for c, n in zip(X_CHUNKS, shape))


def x_layout_problems(x):
    """What is wrong with a series X dataset (h5py), or [] if nothing."""
    problems = []
    if x.dtype != np.float32:
        problems.append(f"X is {x.dtype}, not float32")
    if x.compression != "gzip":
        problems.append(f"X compression is {x.compression or 'none'}, not gzip")
    if x.chunks != x_chunks(x.shape):
        problems.append(f"X chunks are {x.chunks}, not {x_chunks(x.shape)}")
    return problems


def substrates(xml):
    """[(name, units)] for each substrate, in the row order of the .mat fields."""
    variables = ET.parse(xml).getroot().find("microenvironment/domain/variables")
    rows = sorted((int(v.get("ID")), v.get("name"), v.get("units", ""))
                  for v in variables.findall("variable"))
    return [(name, units) for _, name, units in rows]


def cell_labels(xml):
    """{label: row} of the cells .mat, read from the output XML."""
    labels = ET.parse(xml).getroot().find(
        "cellular_information//simplified_data/labels")
    return {label.text: int(label.get("index")) for label in labels}


def write_microenvironment(steps, path, times):
    """Stack the substrate fields into one (time, row, voxel) array.

    Read straight from the .mat rather than through pcdl: the fields are a plain
    grid with no cell axis, so a DataFrame round-trip would only make them larger
    and less obviously what they are. Rows are voxel x, y, z, volume, then one
    per substrate, named in `row_names` from the run's own XML.

    `steps` are the series' output XMLs, so field i is at the series' time i by
    construction rather than by two globs happening to sort alike.
    """
    import h5py

    files = [xml.with_name(xml.stem + "_microenvironment0.mat") for xml in steps]
    missing = [f.name for f in files if not f.exists()]
    if len(missing) == len(files):
        return None
    if missing:
        raise AssertionError(
            f"{len(missing)} of {len(files)} timesteps have no microenvironment "
            f"file, e.g. {missing[:3]}")

    fields = []
    for f in files:
        mat = scipy.io.loadmat(f)
        fields.append(mat[[k for k in mat if not k.startswith("__")][0]].astype("float32"))
    stack = np.stack(fields)

    rows = list(VOXEL_ROWS) + substrates(steps[0])
    if len(rows) != stack.shape[1]:
        raise AssertionError(
            f"the .mat fields have {stack.shape[1]} rows but the XML names "
            f"{len(rows)}: {[name for name, _ in rows]}")
    with h5py.File(path, "w") as h:
        d = h.create_dataset("field", data=stack, compression="gzip", chunks=True)
        d.attrs["layout"] = "(timepoint, row, voxel); rows named in row_names"
        h.create_dataset("row_names", data=[name for name, _ in rows],
                         dtype=h5py.string_dtype())
        h.create_dataset("row_units", data=[units for _, units in rows],
                         dtype=h5py.string_dtype())
        h.create_dataset("time", data=np.asarray(times, dtype="float64"))
    return stack.shape


def verify_unscaled(h5ad_path, mat_path, xml_path, var_name="total_volume"):
    """Check a stored column against the raw .mat it came from.

    This is the assertion the whole archive rests on. pcdl scales by default, and a
    silently normalised archive looks perfectly fine until someone reads an axis
    label -- at which point the only fix is re-reading 110 GB. Cheap to check, so
    it is checked every run rather than trusted.

    Cells are matched by ID, and only the ones at the .mat's own time are
    compared, so the same check covers a snapshot and one timestep of the series.
    Row positions come from the XML's labels rather than being assumed.
    """
    import h5py

    labels = cell_labels(xml_path)
    when = int(float(ET.parse(xml_path).getroot().find("metadata/current_time").text))

    with h5py.File(h5ad_path, "r") as h:
        names = [n.decode() if isinstance(n, bytes) else str(n)
                 for n in h["var"][h["var"].attrs.get("_index", "_index")][:]]
        if var_name not in names:
            raise AssertionError(f"{var_name!r} not stored in {h5ad_path.name}")
        column = h["X"][:, names.index(var_name)]
        obs = [n.decode() if isinstance(n, bytes) else str(n)
               for n in h["obs"][h["obs"].attrs.get("_index", "_index")][:]]

    # obs_names are <cell id>_<minutes> (name_observations)
    ids, times = zip(*(n.rsplit("_", 1) for n in obs))
    at = np.asarray(times) == str(when)
    stored = pd.Series(column[at], index=np.asarray(ids, dtype=int)[at]).sort_index()

    cells = scipy.io.loadmat(mat_path)["cells"]
    raw = pd.Series(cells[labels[var_name], :],
                    index=cells[labels["ID"], :].astype(int)).sort_index()

    if not stored.index.equals(raw.index):
        raise AssertionError(
            f"{h5ad_path.name} at t={when} holds {len(stored):,} cells, the raw "
            f".mat {len(raw):,}, or not the same ones")
    if not np.allclose(stored.values, raw.values, rtol=1e-5, atol=1e-5):
        raise AssertionError(
            f"{var_name} in {h5ad_path.name} does not match the raw .mat -- "
            f"stored range [{stored.min():.4g}, {stored.max():.4g}], "
            f"raw range [{raw.min():.4g}, {raw.max():.4g}]. "
            "If stored is 0-1, scale= was not None."
        )
    return stored.min(), stored.max()


def name_observations(adata, time):
    """Make obs_names unique and meaningful: <cell id>_<simulation minutes>.

    Cell ids repeat at every timestep, so a concatenated series would otherwise
    carry one row called "1234" per timestep. Suffixing with the timepoint makes
    each row addressable, and -- because the snapshots use the same rule -- the
    same cell at the same moment has the same name in the series and in the
    snapshot, so the two join directly.

    Time comes from the timestep rather than a file index: final.xml is not
    always the last numbered output. Simulation 10's is 4,425 cells against
    output00000336's 4,419.
    """
    ids = [str(x) for x in adata.obs_names]
    adata.obs_names = [f"{i}_{int(time)}" for i in ids]
    return adata


def provenance():
    """What produced this file. The archive exists so the raw output is never
    re-read; it should be able to say what made it."""
    import datetime
    from importlib.metadata import version

    import pcdl

    return {
        "pcdl_version": pcdl.__version__,
        "anndata_version": version("anndata"),
        "extracted": datetime.datetime.now().isoformat(timespec="seconds"),
    }


def stamp_identity(path, meta):
    """Write simulation, sample and patient into a file's uns group.

    Cheap enough to apply to a finished file, which is what lets an archive
    written before this existed be corrected without re-converting it.
    """
    import h5py

    values = {k: ("" if v is None else str(v)) for k, v in meta.items()}
    values.update(provenance())
    with h5py.File(path, "a") as h:
        uns = h.require_group("uns")
        uns.attrs["encoding-type"] = "dict"
        uns.attrs["encoding-version"] = "0.1.0"
        for key, value in values.items():
            if key in uns:
                del uns[key]
            d = uns.create_dataset(key, data=value)
            d.attrs["encoding-type"] = "string"
            d.attrs["encoding-version"] = "0.2.0"


def verify_placement(simulation_id, output_dir, report):
    """Compare the requested initial composition against the t=0 output.

    Optional, and off by default, because it checks PhysiCell rather than this
    conversion -- the counts agreed exactly on every simulation spot-checked. It
    earns its place when a run may have failed to place cells, which the archive
    would otherwise record as a real starting condition.
    """
    import collections

    sys.path.insert(0, str(BASE))
    from scripts.resolve_samples import _variation_counts

    with sqlite3.connect(PCMM_DB) as con:
        row = con.execute("SELECT ic_cell_variation_id FROM simulations "
                          "WHERE simulation_id = ?", (simulation_id,)).fetchone()
    variations = _variation_counts(BASE)
    if not row or int(row[0]) not in variations:
        report["placement"] = "no variation to compare"
        return
    requested = variations[int(row[0])]

    mat = output_dir / "initial_cells.mat"
    xml = output_dir / "output00000000.xml"
    if not (mat.exists() and xml.exists()):
        report["placement"] = "no t=0 output to compare"
        return
    codes = scipy.io.loadmat(mat)["cells"][cell_labels(xml)["cell_type"], :].astype(int)
    names = {int(c.get("ID")): c.text
             for c in ET.parse(xml).getroot().find(".//cell_types")}
    actual = collections.Counter(names.get(c, "?") for c in codes)

    differ = {t: (n, actual.get(t, 0)) for t, n in requested.items()
              if n != actual.get(t, 0)}
    if differ:
        raise AssertionError(
            f"simulation {simulation_id}: PhysiCell placed different counts than "
            f"the initial condition requested -- {differ}"
        )
    report["placement"] = f"{sum(requested.values()):,} cells as requested"


def verify_complete(paths, steps, report):
    """Refuse to call a simulation done unless its files say so.

    This exists because the failure mode here is not a crash. Three separate
    times during development -- an out-of-memory kill, a missing dask, and a bug
    in a hand-written writer -- the script exited 0 having written some files but
    not others. Every time the output looked plausible. For a job whose whole
    purpose is never re-reading the raw output, a silent partial success is worse
    than an error, because nothing downstream will question it.

    Checks that every file exists and is non-trivial, that the series holds
    exactly as many rows as the timesteps have cells between them, and that X is
    float32 and gzip-compressed in the chunks it was meant to have.
    """
    import h5py

    for name, path in paths.items():
        if not path.exists():
            raise AssertionError(f"{path.name} was not written")
        if path.stat().st_size < 1024:
            raise AssertionError(f"{path.name} is {path.stat().st_size} bytes")

    # whosmat reads the MAT header only. loadmat here meant a second full pass
    # over 917 MB per simulation purely to count columns -- about a quarter of
    # the runtime, spent re-reading data already converted.
    expected = 0
    for xml in steps:
        mat = xml.with_name(xml.stem + "_cells.mat")
        if mat.exists():
            expected += scipy.io.whosmat(str(mat))[0][1][1]

    # The layout, because the worst bugs so far were invisible to every other
    # check. concat_on_disk writes X uncompressed and in float64, and h5repack
    # skips a chunk larger than the data without a word; either way the files
    # existed, the structure was right and the row count matched. A size limit
    # used to stand in for this and could not see a small uncompressed series.
    with h5py.File(paths["series"], "r") as h:
        problems = x_layout_problems(h["X"])
        rows = h["X"].shape[0]
        uns = {k: h["uns"][k][()] for k in ("sample_id", "pcdl_version")
               if "uns" in h and k in h["uns"]}
    if problems:
        raise AssertionError(f"{paths['series'].name}: {'; '.join(problems)}")
    if rows != expected:
        raise AssertionError(
            f"series has {rows:,} rows but the timesteps hold {expected:,} cells "
            "between them -- the concatenation dropped or duplicated data"
        )
    for key in ("sample_id", "pcdl_version"):
        if not uns.get(key):
            raise AssertionError(
                f"{paths['series'].name} carries no {key} -- uns was not stamped"
            )
    report["rows_verified"] = rows
    report["series_MB"] = round(paths["series"].stat().st_size / 2**20, 1)


def write_series(steps, meta, out, scratch):
    """Every timestep as one AnnData at `out`. Returns the timestep times.

    NOT mcdsts.get_anndata(collapse=True). That path loads every timestep into
    memory and then builds the concatenated frame: at the 30-minute saves this
    was written against, 2.55M rows x 459 float64 columns was 8.7 GB for the
    result alone, before pcdl's per-timestep objects, and it died partway
    through -- silently, having already written the snapshots, which is how it
    first showed up as "exit 0 but no series file".

    Streaming instead keeps peak memory at one timestep.
    """
    import h5py
    import pcdl

    times, parts = [], []
    for i, xml in enumerate(steps):
        mcds = pcdl.TimeStep(str(xml), microenv=False, graph=False,
                             physiboss=False, verbose=False)
        step = mcds.get_anndata(values=VALUES, scale=SCALE)
        name_observations(step, mcds.get_time())
        # Identity goes on every part, not onto the finished file. uns is
        # per-file and cannot survive a concat into a combined object, so a
        # row has to carry its own -- and adding it here means concat_on_disk
        # propagates it rather than needing h5py surgery afterwards.
        for field, value in meta.items():
            step.obs[field] = "" if value is None else str(value)
        # pcdl hands back float64. These are simulation state variables, not
        # quantities with 15 significant figures; float32 halves the file the
        # concat has to write and costs nothing anyone will measure.
        step.X = step.X.astype("float32")
        times.append(float(mcds.get_time()))
        part = scratch / f"part_{i:06d}.h5ad"
        step.write_h5ad(part)          # uncompressed: read once, then deleted
        parts.append(part)
        del step, mcds

    # anndata's own on-disk concatenation. Hand-writing this was tried and
    # abandoned: with pandas 3.0 in this environment, anndata encodes string
    # columns as nullable-string-array groups (mask + values) and categoricals
    # as codes + categories groups, so the on-disk layout is both intricate and
    # dependent on the pandas version doing the writing. That is a format to
    # let its own library own. Needs dask, which is what its dense path uses.
    from anndata.experimental import concat_on_disk

    raw_out = scratch / "series.raw.h5ad"
    concat_on_disk([str(x) for x in parts], str(raw_out), axis=0, join="outer")

    # Unlink the parts now. They and the uncompressed concat are each about
    # rows x 457 x 4 bytes, so holding both through the repack doubles peak
    # scratch per simulation, on a shared mount.
    for part in parts:
        part.unlink()

    # concat_on_disk writes X uncompressed and unchunked -- it takes no
    # compression argument. h5repack applies gzip in a streaming pass, so peak
    # memory stays flat.
    # Scoped to /X on purpose. A bare "-f GZIP=4" leaves obs, obsm and var
    # uncompressed anyway -- which is what we want -- but by accident rather
    # than instruction, and that is not a thing to depend on across h5repack
    # versions.
    with h5py.File(raw_out, "r") as h:
        chunks = x_chunks(h["X"].shape)
    result = subprocess.run(
        ["h5repack",
         "-f", "/X:GZIP=4",
         "-l", f"/X:CHUNK={chunks[0]}x{chunks[1]}",
         str(raw_out), str(out)],
        capture_output=True, text=True,
    )
    if result.returncode:
        raise RuntimeError(f"h5repack exited {result.returncode}: "
                           f"{(result.stderr or result.stdout).strip()}")
    raw_out.unlink()
    with h5py.File(out, "r") as h:
        problems = x_layout_problems(h["X"])
    if problems:
        raise AssertionError(f"h5repack left {out.name} wrong: {'; '.join(problems)}")

    # concat_on_disk does not carry uns through, so identity is stamped in
    # afterwards. Without it a derived file cannot say which sample it is,
    # which is the whole reason it is stamped rather than joined later.
    stamp_identity(out, meta)
    return times


def extract(simulation_id, force=False, verify_counts=False):
    """Write the four derived files for one simulation. Returns a timing report."""
    import pcdl

    output_dir = OUTPUTS / str(simulation_id) / "output"
    if not output_dir.is_dir():
        raise FileNotFoundError(f"no output directory for simulation {simulation_id}")
    DERIVED.mkdir(parents=True, exist_ok=True)

    meta = describe(simulation_id)
    if meta["sample_id"] is None:
        raise AssertionError(
            f"simulation {simulation_id} has no sample; refusing to write a file "
            "that cannot say what it is"
        )
    paths = {
        k: DERIVED / f"{meta['sim_id']}-{k}.{ext}"
        for k, ext in [("initial", "h5ad"), ("final", "h5ad"),
                       ("series", "h5ad"), ("microenv", "h5")]
    }
    if not force and all(p.exists() for p in paths.values()):
        return {"simulation": simulation_id, "skipped": True}

    report = {"simulation": simulation_id, "sample": meta["sample_id"], **meta}
    steps = sorted(output_dir.glob("output*.xml"))

    # Everything is written here and moved into derived/ only after the checks
    # pass. A leftover folder is from a task that was killed outright, so it is
    # cleared rather than trusted.
    scratch = DERIVED / f".tmp_sim_{simulation_id}"
    shutil.rmtree(scratch, ignore_errors=True)
    scratch.mkdir(parents=True)
    work = {k: scratch / p.name for k, p in paths.items()}
    try:
        # -- snapshots: TimeStep is not the collapsed path, so obsp keeps the graphs --
        # graph=False at t=0: no neighbours are recorded there, so obsp comes back
        # empty regardless and loading three graph files is a guaranteed discard.
        for key, want_graph in (("initial", False), ("final", True)):
            t0 = time.time()
            mcds = pcdl.TimeStep(str(output_dir / f"{key}.xml"), microenv=False,
                                 graph=want_graph, physiboss=False, verbose=False)
            adata = mcds.get_anndata(values=VALUES, scale=SCALE)
            name_observations(adata, mcds.get_time())
            # float32 here too, matching the series. Not for the 50 MB it saves
            # across the archive, but because a snapshot and the series it belongs
            # to should concatenate without a dtype promotion.
            adata.X = adata.X.astype("float32")
            adata.uns.update(provenance())
            # `field`, not `key`: the enclosing loop binds `key` to the snapshot
            # name and uses it for work[key] afterwards.
            for field, value in meta.items():
                adata.uns[field] = "" if value is None else str(value)
                adata.obs[field] = "" if value is None else str(value)
            adata.write_h5ad(work[key], compression="gzip")
            report[f"{key}_s"] = round(time.time() - t0, 1)
            report[f"{key}_cells"] = adata.n_obs
            report[f"{key}_graphs"] = sorted(adata.obsp.keys())

        # -- series: one timestep at a time, concatenated on disk --
        t0 = time.time()
        times = write_series(steps, meta, work["series"], scratch)
        report["series_s"] = round(time.time() - t0, 1)
        report["series_steps"] = len(steps)

        import anndata as ad

        series = ad.read_h5ad(work["series"], backed="r")
        report["series_rows"], report["series_cols"] = series.n_obs, series.n_vars
        report["has_spatial"] = "spatial" in series.obsm
        series.file.close()

        # -- microenvironment: its own file, no cell axis to annotate --
        t0 = time.time()
        report["microenv_shape"] = write_microenvironment(steps, work["microenv"], times)
        report["microenv_s"] = round(time.time() - t0, 1)

        # -- the assertion the archive rests on, on a snapshot and on the series --
        lo, hi = verify_unscaled(work["initial"], output_dir / "initial_cells.mat",
                                 output_dir / "initial.xml")
        report["total_volume_range"] = (round(float(lo), 1), round(float(hi), 1))
        last = steps[-1]
        verify_unscaled(work["series"], last.with_name(last.stem + "_cells.mat"), last)

        verify_complete(work, steps, report)
        if verify_counts:
            verify_placement(simulation_id, output_dir, report)

        for key, path in paths.items():
            work[key].replace(path)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    report["MB"] = round(sum(p.stat().st_size for p in paths.values()) / 2**20, 1)
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--simulation", type=int, action="append",
                    help="simulation id; repeatable. Default: all Completed.")
    ap.add_argument("--verify-counts", action="store_true",
                    help="also check that PhysiCell placed the cell counts the "
                         "initial condition asked for")
    ap.add_argument("--force", action="store_true",
                    help="re-extract simulations whose files already exist")
    args = ap.parse_args()

    ids = args.simulation or completed_simulations()
    print(f"extracting {len(ids)} simulation(s) into {DERIVED}")
    for sid in ids:
        started = time.time()
        report = extract(sid, force=args.force, verify_counts=args.verify_counts)
        if report.get("skipped"):
            print(f"  sim {sid}: already extracted, skipping")
            continue
        detail = report["sample"]
        if report.get("geometry"):
            detail += f" / {report['geometry']}"
        print(f"  {report['sim_type']} sim {sid} ({detail}): {report['MB']} MB "
              f"in {time.time() - started:.0f}s")
        for k in ("sim_id", "caf_mhc2_rate", "placement", "initial_s", "final_s",
                  "series_s", "microenv_s", "series_steps", "series_rows",
                  "series_cols", "has_spatial", "initial_graphs", "final_graphs",
                  "microenv_shape", "total_volume_range"):
            print(f"      {k:20s} {report.get(k)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
