"""
imc_tm.py
---------
Shared helpers for the JHH IMC -> TissueMosaic workflow. The notebooks and the
CLI scripts in this directory both import from here, so the conversion logic,
the channel definition and the plotting conventions exist in exactly one place.

    import sys; sys.path.insert(0, "<repo>/tissuemosaic")
    import imc_tm

IMPORTANT -- do not add an __init__.py to this directory. It is named
`tissuemosaic`, and turning it into a package makes it shadow the INSTALLED
TissueMosaic library: `from tissuemosaic.data import AnndataFolderDM` then fails
with ModuleNotFoundError. Today it is safe only because a directory without
__init__.py is a namespace-package candidate, which loses to the regular package
in site-packages. Keep it that way; import these modules flat.

Importing this module applies `tm_compat`, which scopes torch's weights_only
handling (see that module). That is a deliberate import side effect: it is the
only way a notebook picks the shim up without remembering to, and it is
reported via `imc_tm.compat_summary()` rather than being silent.

Heavy dependencies (matplotlib, seaborn, torch, tissuemosaic) are imported
lazily inside the functions that need them, so `import imc_tm` stays cheap
enough for the training entry points to use.

Coordinate frame: microns, origin at the ROI centre, y increasing UPWARD.
See the README and notebook 0 for the derivation. Plot y as-is.
"""

import os
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import tm_compat

# ---------------------------------------------------------------- paths -----


def repo_root(start=None):
    """Walk upward until the directory containing prep_imc_spatial/ is found.

    Works whether the caller is a script in tissuemosaic/ or a notebook whose
    cwd is tissuemosaic/ or the repo root.
    """
    start = Path(start) if start is not None else Path(__file__).resolve().parent
    for p in [start, *start.parents]:
        if (p / "prep_imc_spatial").is_dir():
            return p
    raise RuntimeError(
        f"repo root (a directory containing prep_imc_spatial/) not found from {start}"
    )


REPO = repo_root()
ICS_DIR = REPO / "PhysiCell/user_projects/antigen_presentation/config/ics/JHH_IMC"
H5AD_DIR = REPO / "tissuemosaic/imc_anndata"
# Small result tables live beside the code, matching prep_imc_spatial/, which commits
# assignmentsummary_JHH_IMC.csv and imc_spatial_roi_specs.csv in its own directory.
# analysis/ in this repo is 10 flat .py files with no tracked subdirectories.
RESULTS_DIR = REPO / "tissuemosaic"

# Training output lives under tissuemosaic/, not at the repo root. Two reasons: this tool
# already writes its other ~36 GB (imc_anndata*/) here, so a top-level runs/ split the same
# tool across two homes; and `runs/<label>/` already means something else in slurm/ -- the
# pdac-spatial-pipeline repo on the cluster, not this one. Keeping it here also makes
# tissuemosaic/ self-contained, the way prep_imc_spatial/ is.
RUNS_DIR = REPO / "tissuemosaic" / "runs"


QUPATH_DIR = REPO / "qupath_detections"
ROI_SPECS = REPO / "prep_imc_spatial/imc_spatial_roi_specs.csv"
CONFIG_YAML = REPO / "tissuemosaic/config_dino_ssl_imc.yaml"


def compat_summary():
    """What tm_compat changed, for reporting at the top of a notebook."""
    return "tm_compat: {} numpy types allowlisted (weights_only stays True); relaxed only for {}".format(
        tm_compat.allowlisted, ", ".join(tm_compat.scoped_files)
    )


# ------------------------------------------------------------- channels -----

# Index in this list IS the image channel index. It must match
# `categories_to_channels` in config_dino_ssl_imc.yaml -- check_channels()
# asserts that, and is called on import.
CELL_TYPES = [
    "CAF",
    "apCAF",
    "CD4_Tcell",
    "CD8_Tcell",
    "epithelial_tumor_class1",
    "epithelial_tumor_class1_class2",
    "mesenchymal_tumor_class1",
    "mesenchymal_tumor_class1_class2",
    "other_tissue",
    "duct_filler",
]

# Cell types that are not measured biology: duct_filler is hex-packed synthetic
# fill, other_tissue is unclassified. Together ~71% of all cells.
BACKGROUND_TYPES = ("other_tissue", "duct_filler")
BIOLOGY_TYPES = [t for t in CELL_TYPES if t not in BACKGROUND_TYPES]

# Variants of an ROI already in the set; the only rows with NaN coordinates and
# the only source of epithelial_normal.
SKIP_PATTERN = "withductfiller"


def channels_for(n):
    """The channel-name list matching an n-column composition/one-hot block.

    A featurized ECM run carries 11 channels, not 10, and every merge helper used to
    assume len(CELL_TYPES). Silently zipping 10 names onto 11 columns drops the ECM
    channel from the plots without erroring, which is the worst possible failure here.
    """
    if int(n) == len(CELL_TYPES):
        return list(CELL_TYPES)
    if int(n) == len(ECM_CELL_TYPES):
        return list(ECM_CELL_TYPES)
    raise ValueError("no channel list for {} channels".format(n))


def categories_to_channels():
    return {t: i for i, t in enumerate(CELL_TYPES)}


def check_channels(config_path=None, strict=True):
    """Assert CELL_TYPES matches the yaml's channel order.

    A silent mismatch here permutes the image channels and trains the model on
    scrambled inputs with no visible error, so this runs on import.
    """
    import yaml

    path = Path(config_path) if config_path else CONFIG_YAML
    if not path.is_file():
        return {"checked": False, "reason": f"{path} not found"}

    cfg = yaml.safe_load(open(path))
    c2c = cfg.get("categories_to_channels", {})
    from_yaml = [t for t, _ in sorted(c2c.items(), key=lambda kv: kv[1])]

    problems = []
    if from_yaml != CELL_TYPES:
        problems.append(
            f"channel order differs:\n  imc_tm : {CELL_TYPES}\n  yaml   : {from_yaml}"
        )
    if cfg.get("image_in_ch") != len(CELL_TYPES):
        problems.append(
            "image_in_ch={} but there are {} channels".format(
                cfg.get("image_in_ch"), len(CELL_TYPES)
            )
        )
    if sorted(c2c.values()) != list(range(len(c2c))):
        problems.append(f"channel indices are not 0..N-1: {sorted(c2c.values())}")

    if problems and strict:
        raise AssertionError(
            f"CELL_TYPES and {path.name} disagree:\n" + "\n".join(problems)
        )
    return {
        "checked": True,
        "ok": not problems,
        "problems": problems,
        "config": str(path),
    }


# ----------------------------------------------------------- conversion -----


def roi_to_patient(roi):
    """JHH317ROI3 -> JHH317, JHH417RROI4 -> JHH417R."""
    m = re.match(r"(JHH\d+R?)ROI", roi)
    if m is None:
        raise ValueError(f"cannot parse patient from ROI name {roi!r}")
    return m.group(1)


def list_ics_csvs():
    """The canonical ICS CSVs, excluding the withductfiller variants."""
    return sorted(p for p in ICS_DIR.glob("*.csv") if SKIP_PATTERN not in p.name)


def load_roi_domains():
    """{roi: [x_min, x_max, y_min, y_max]} in microns, from the ROI spec table."""
    if not ROI_SPECS.is_file():
        return {}
    spec = pd.read_csv(ROI_SPECS)
    return {
        r["roi"]: [r["x_min"], r["x_max"], r["y_min"], r["y_max"]]
        for _, r in spec.iterrows()
    }


def csv_to_anndata(csv_path, status, domain=None):
    """One ICS CSV -> AnnData in the form SparseImage.from_anndata expects."""
    csv_path = Path(csv_path)
    df = pd.read_csv(csv_path)
    n_raw = len(df)
    df = df.dropna(subset=["x", "y", "type"])
    n_dropped = n_raw - len(df)

    unknown = set(df["type"]) - set(CELL_TYPES)
    if unknown:
        raise ValueError(f"{csv_path.name}: unexpected cell types {sorted(unknown)}")

    df = df.reset_index(drop=True)
    index = pd.Index(df.index.astype(str))

    # Reindexed onto the full channel list so absent types still get a zeroed
    # column -- SparseImage keys channels off these column names.
    onehot = (
        pd.get_dummies(df["type"])
        .reindex(columns=CELL_TYPES, fill_value=0)
        .astype(np.float32)
    )
    onehot.index = index

    obs = pd.DataFrame(
        {
            "x": df["x"].astype(float).values,
            "y": df["y"].astype(float).values,
            "cell_type": pd.Categorical(df["type"].values, categories=CELL_TYPES),
            "volume": df["volume"].astype(float).values,
        },
        index=index,
    )

    roi = csv_path.stem
    adata = ad.AnnData(obs=obs)
    adata.obsm["cell_type_proportions"] = onehot
    adata.uns["status"] = int(status)
    adata.uns["patient"] = roi_to_patient(roi)
    adata.uns["roi"] = roi
    # Coordinate frame provenance -- verified by exact reconstruction from the
    # QuPath detections; see notebook 0 section 2.
    adata.uns["coord_units"] = "micron"
    adata.uns["coord_origin"] = "ROI centre (0, 0)"
    adata.uns["coord_y_axis"] = "increases upward (flipped from QuPath image y-down)"
    adata.uns["coord_source"] = (
        "QuPath 'Centroid X/Y um' -> assemble_initial_conditions.py"
    )
    if domain is not None:
        adata.uns["roi_domain_um"] = np.asarray(domain, dtype=float)
    adata.uns["n_dropped_nan_rows"] = int(n_dropped)
    return adata


def build_anndata(out_dir=None, verbose=True):
    """Rebuild every .h5ad from the canonical ICS CSVs.

    DELETES any existing .h5ad and the train/test_dataset.pt caches in out_dir
    first, so the result never mixes generations. Returns a summary DataFrame.
    """
    out_dir = Path(out_dir) if out_dir else H5AD_DIR
    csvs = list_ics_csvs()
    if not csvs:
        raise RuntimeError(f"no ICS CSVs found in {ICS_DIR}")

    patients = sorted({roi_to_patient(p.stem) for p in csvs})
    status_of = {p: i for i, p in enumerate(patients)}
    domains = load_roi_domains()

    out_dir.mkdir(parents=True, exist_ok=True)
    removed = 0
    for stale in list(out_dir.glob("*.h5ad")) + list(out_dir.glob("*.pt")):
        stale.unlink()
        removed += 1

    if verbose:
        print(f"{len(csvs)} ROIs from {len(patients)} patients")
        print(f"ROI domains loaded for {len(domains)} ROIs")
        if removed:
            print(f"removed {removed} stale file(s) from {out_dir}")

    rows = []
    for csv_path in csvs:
        roi = csv_path.stem
        patient = roi_to_patient(roi)
        adata = csv_to_anndata(
            csv_path, status=status_of[patient], domain=domains.get(roi)
        )
        adata.write_h5ad(out_dir / f"{roi}.h5ad")
        rows.append(
            {
                "roi": roi,
                "patient": patient,
                "status": status_of[patient],
                "n_cells": adata.n_obs,
            }
        )

    summary = pd.DataFrame(rows)
    if verbose:
        print(
            f"wrote {len(summary)} h5ad files, {summary.n_cells.sum():,} cells total, to {out_dir}"
        )
        print("channels ({}): {}".format(len(CELL_TYPES), ", ".join(CELL_TYPES)))
    return summary


def load_all_anndata(h5ad_dir=None):
    """{roi: AnnData} for every .h5ad in the directory."""
    from anndata import read_h5ad

    d = Path(h5ad_dir) if h5ad_dir else H5AD_DIR
    return {
        f[:-5]: read_h5ad(d / f) for f in sorted(os.listdir(d)) if f.endswith(".h5ad")
    }


def roi_metadata(adatas):
    """Per-ROI summary table from a {roi: AnnData} mapping."""
    return (
        pd.DataFrame(
            [
                {
                    "roi": r,
                    "patient": a.uns["patient"],
                    "status": a.uns["status"],
                    "n_cells": a.n_obs,
                    "x_span": float(a.obs.x.max() - a.obs.x.min()),
                    "y_span": float(a.obs.y.max() - a.obs.y.min()),
                }
                for r, a in adatas.items()
            ]
        )
        .sort_values("roi")
        .reset_index(drop=True)
    )


def composition(adatas, fraction=False):
    """ROI x cell-type counts (or fractions), columns in CELL_TYPES order."""
    comp = pd.DataFrame(
        {r: a.obs["cell_type"].value_counts() for r, a in adatas.items()}
    ).T
    comp = comp.reindex(columns=CELL_TYPES).fillna(0).astype(int)
    return comp.div(comp.sum(axis=1), axis=0) if fraction else comp


# ------------------------------------------------- config / model / data -----


def load_config(path=None):
    import yaml

    return yaml.safe_load(open(Path(path) if path else CONFIG_YAML))


def make_datamodule(config=None, **overrides):
    """AnndataFolderDM built from the IMC yaml, with optional overrides."""
    from tissuemosaic.data import AnndataFolderDM

    cfg = config if isinstance(config, dict) else load_config(config)
    kwargs = AnndataFolderDM.get_default_params()
    for k in (
        "data_folder",
        "categories_to_channels",
        "category_key",
        "pixel_size",
        "x_key",
        "y_key",
        "status_key",
        "global_size",
        "local_size",
        "n_element_min_for_crop",
        "n_neighbours_moran",
    ):
        if k in cfg:
            kwargs[k] = cfg[k]
    kwargs.update(overrides)
    return AnndataFolderDM(**kwargs)


def load_model(ckpt_path, eval_mode=True):
    """Load whichever SSL model the checkpoint declares."""
    import torch
    from tissuemosaic.models.ssl_models import Barlow, Dino, Simclr, Vae

    models = {"barlow": Barlow, "dino": Dino, "simclr": Simclr, "vae": Vae}
    raw = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    name = raw["hyper_parameters"]["ssl_model"]
    model = models[name].load_from_checkpoint(
        checkpoint_path=str(ckpt_path), strict=False
    )
    return model.eval() if eval_mode else model


def to_numpy(x):
    """Tensors may be sparse and/or on GPU; ndarrays pass through."""
    import torch

    if isinstance(x, torch.Tensor):
        if x.is_sparse:
            x = x.to_dense()
        return x.detach().cpu().numpy()
    return np.asarray(x)


def tissuemosaic_run_script(name):
    """Path to one of upstream TissueMosaic's run/ entry points.

    Override the checkout location with TISSUEMOSAIC_HOME.
    """
    home = Path(os.environ.get("TISSUEMOSAIC_HOME", Path.home() / "TissueMosaic"))
    entry = home / "run" / name
    if not entry.is_file():
        raise FileNotFoundError(
            f"cannot find {entry}\nSet TISSUEMOSAIC_HOME to the TissueMosaic checkout."
        )
    return entry


# --------------------------------------------------- window-size sweep -----
#
# Supplementary experiment S1: which window size (field of view) and which
# rasterisation resolution give the most useful micro-environment embedding?
#
# FOV = pixel_size * global_size, so the two knobs are not independent. The four
# configs form an L sharing fov192_px2.0 (the production config) as the hinge:
#
#   A vs B  -- FOV 128 vs 192 at fixed pixel_size 2.0   "how much tissue?"
#   C vs B vs D -- FOV 192 at pixel_size 3.0/2.0/1.5    "how finely rasterised?"
#   A vs C  -- FOV 128 vs 192 at IDENTICAL network input (64/48), so scale is
#              isolated from model capacity and compute.
#
# global_size / local_size are multiples of 16: the backbone downsamples 16x, so
# these give clean 4x4/6x6/8x8 pre-pool feature maps. (AdaptiveAvgPool2d means
# other sizes do not error, but they downsample raggedly.)

WINDOW_SWEEP = [
    {
        "name": "A_fov128_px2.0",
        "fov_um": 128,
        "pixel_size": 2.0,
        "global_size": 64,
        "local_size": 48,
    },
    {
        "name": "B_fov192_px2.0",
        "fov_um": 192,
        "pixel_size": 2.0,
        "global_size": 96,
        "local_size": 64,
    },
    {
        "name": "C_fov192_px3.0",
        "fov_um": 192,
        "pixel_size": 3.0,
        "global_size": 64,
        "local_size": 48,
    },
    {
        "name": "D_fov192_px1.5",
        "fov_um": 192,
        "pixel_size": 1.5,
        "global_size": 128,
        "local_size": 80,
    },
    # Round 2 (2026-08-20). Round 1 found 128 um beat every 192 um variant, and 128 was
    # the SMALLEST window tested -- so the trend pointed below it, not above. These
    # extend downward at the settled 2.0 um/px. Pre-pool feature maps stay non-degenerate
    # (gs 32 -> 2x2, gs 48 -> 3x3, local 24 -> 2x2); only an input of 16 collapses to 1x1.
    {
        "name": "E_fov64_px2.0",
        "fov_um": 64,
        "pixel_size": 2.0,
        "global_size": 32,
        "local_size": 24,
    },
    {
        "name": "F_fov96_px2.0",
        "fov_um": 96,
        "pixel_size": 2.0,
        "global_size": 48,
        "local_size": 32,
    },
]

# JHH369 mean cell density, used to scale the analysis window filter with FOV.
SWEEP_DENSITY_PER_MM2 = 6348.0


def fov_scaled_threshold(fov_um, density_per_mm2=None, fraction=0.5):
    """Analysis window filter that scales with field of view.

    A FIXED threshold is itself FOV-dependent, which is the bias it was meant to
    prevent: 50 cells is 48% of an expected 128 um window but only 21% of a 192 um
    one, so it stripped borders at 128 um and none at 192 um (retained border share
    32% vs 51%). Border windows hold ~40% of the expected count at any scale, so a
    threshold at `fraction` of expected separates them consistently.

    At fraction=0.5 this reproduces ~52 for 128 um -- i.e. the 50 used in round 1 --
    while giving 13 / 29 / 117 for 64 / 96 / 192 um.
    """
    d = SWEEP_DENSITY_PER_MM2 if density_per_mm2 is None else density_per_mm2
    return max(1, int(round(fraction * d * (fov_um / 1000.0) ** 2)))


SWEEP_PATIENT = "JHH369"
SWEEP_H5AD_DIR = REPO / "tissuemosaic/imc_anndata_jhh369"
SWEEP_RUN_DIR = RUNS_DIR / "window_sweep"

# Only max_epochs deviates from the documented recipe on purpose. warm_up/
# warm_down are rescaled to keep the schedule SHAPE: milestones are
# (0, warm_up, max_epochs - warm_down, max_epochs), so leaving them at 100/100
# with max_epochs=200 would collapse the peak-LR plateau to zero width -- and it
# would not error, because the assert is warm_up + warm_down <= max_epochs.
# ======================================================================= FINAL =====
# The confirmed final sweep: two FOVs, 8 ROIs, full-length training.
#
# Differs from the S1 screen in three ways that matter:
#   1. 1000 epochs, not 200 -- these are meant to converge, not to rank.
#   2. 8 ROIs across 5 patients, not 4 ROIs of one. The single biggest open question
#      is that a JHH369-only model scored gain 0.33-0.78x on unseen patients, i.e.
#      WORSE than clustering cell counts. This is the test of whether that is a
#      data-diversity problem.
#   3. warm_up/warm_down 100/100, which is TissueMosaic's own shape at 1000 epochs
#      (their config, not their argparse defaults). The S1 screen used 20/100, a 50%
#      decay tail, because 100/100 at 200 epochs would have collapsed the plateau.
FINAL_ROIS = [
    # the S1 set -- one patient, mid-range density, kept for continuity
    "JHH369ROI1",
    "JHH369ROI2",
    "JHH369ROI3",
    "JHH369ROI4",
    # high quality
    "JHH387ROI1",  # 10,590 cells/mm2, 48.4% biology -- dense and signal-rich
    "JHH380ROI1",  #  8,262 cells/mm2, 62.1% biology -- richest classified signal
    # low quality, on ORTHOGONAL axes so the two failure modes stay separable
    "JHH357ROI5",  #  7,396 cells/mm2 but only  2.2% biology -- dense, signal-poor
    "JHH317ROI2",  #  2,643 cells/mm2 but      43.7% biology -- sparse, signal-rich
]

FINAL_SWEEP = [
    {
        "name": "final_fov96",
        "fov_um": 96,
        "pixel_size": 2.0,
        "global_size": 48,
        "local_size": 32,
    },
    {
        "name": "final_fov128",
        "fov_um": 128,
        "pixel_size": 2.0,
        "global_size": 64,
        "local_size": 48,
    },
]

FINAL_OVERRIDES = {
    "max_epochs": 1000,
    "warm_up_epochs": 100,  # (0, 100, 900, 1000): 10% warm-up, 80% plateau,
    "warm_down_epochs": 100,  # 10% cosine decay -- upstream's shape verbatim
    "checkpoint_interval_epochs": 25,  # 40 checkpoints per run, ~334 MB each
    "n_crops_for_tissue_train": 1024,
    "n_element_min_for_crop": 10,  # SAME at train and inference, by decision
    "min_weight_decay": 0.0,
    "max_weight_decay": 0.0,
    "batch_size_per_gpu": 128,
}

FINAL_RUN_DIR = RUNS_DIR / "final_sweep"


def final_h5ad_dir(name):
    """One h5ad folder PER CONFIG, deliberately.

    AnndataFolderDM caches its raster to <data_folder>/train_dataset.pt. Both final
    configs share pixel_size 2.0, so a shared folder would have them racing to write
    the same file when run in parallel. Separate folders cost ~38 MB each.
    """
    return REPO / "tissuemosaic" / f"imc_anndata_{name}"


def build_final_h5ad_dirs(overwrite=False):
    """Populate one h5ad folder per final config with the FINAL_ROIS subset."""
    import shutil

    src = H5AD_DIR
    made = []
    for spec in FINAL_SWEEP:
        d = final_h5ad_dir(spec["name"])
        if d.is_dir() and not overwrite:
            made.append((d, "exists"))
            continue
        d.mkdir(parents=True, exist_ok=True)
        for f in d.glob("*.pt"):  # stale raster cache would be silently wrong
            f.unlink()
        for roi in FINAL_ROIS:
            shutil.copyfile(src / f"{roi}.h5ad", d / f"{roi}.h5ad")
        made.append((d, f"built {len(FINAL_ROIS)} ROIs"))
    return made


# ========================================================================= ECM =====
# ECM enters as an 11th channel. It is NOT a cell: PhysiCell models it as a BioFVM
# substrate, a static scalar field on a 20 um Cartesian voxel mesh with
# diffusion_coefficient = 0 and decay_rate = 0, read by cell rules (it drives EMT and
# modulates CAF migration) but never written. So the IC file IS the field, permanently.
#
# It is given its OWN rows at the voxel centres rather than being attached to cells.
# Attaching to cells would discard the field wherever no cell was segmented, and ECM is
# measurably HIGHER there (7.15 vs 5.92 in JHH317ROI2) -- cell-poor matrix-rich regions
# are desmoplastic stroma, the thing of interest. That loss would also be
# density-dependent: 54% of voxels are cell-free in JHH317ROI2 against 4% in JHH387ROI1.
#
# Values come from _ecm_scaled.csv (what the PhysiCell XMLs actually load) divided by
# 10. Every one of the 51 files maxes at exactly 10.0, so this gives [0, 1] with no
# clipping. NOTE the scaling is PER-ROI -- scaleimage.py clips at each ROI's own 95th
# percentile -- so ECM values are NOT comparable across ROIs.
ECM_CELL_TYPES = CELL_TYPES + ["ecm"]
ECM_DIR = REPO / "PhysiCell/user_projects/antigen_presentation/config/ics/substrates"
ECM_DIVISOR = 10.0
ECM_VOXEL_UM = 20.0


def ecm_csv_path(roi):
    return ECM_DIR / "{}_ecm_scaled.csv".format(roi)


def build_ecm_anndata(rois, out_dir, src_dir=None, verbose=True):
    """Write 11-channel h5ad: the original cells plus one row per ECM voxel centre.

    Cell rows keep their one-hot in channels 0-9 and 0 in channel 10; ECM rows are 0 in
    0-9 and carry ecm/10 in channel 10. The rasteriser drops zero values, so each ECM
    row contributes exactly one entry and never mixes with a cell type.
    """
    from anndata import read_h5ad
    src = Path(src_dir) if src_dir else H5AD_DIR
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for f in list(out.glob("*.h5ad")) + list(out.glob("*.pt")):
        f.unlink()                       # stale raster cache would be silently wrong

    rows = []
    for roi in rois:
        a = read_h5ad(src / "{}.h5ad".format(roi))
        e = pd.read_csv(ecm_csv_path(roi))
        ecm_v = e["ecm"].to_numpy(dtype=np.float32) / ECM_DIVISOR
        assert ecm_v.min() >= 0.0 and ecm_v.max() <= 1.0, "ECM outside [0,1] for {}".format(roi)

        n_c, n_e = a.n_obs, len(e)
        comp = np.zeros((n_c + n_e, len(ECM_CELL_TYPES)), dtype=np.float32)
        comp[:n_c, :len(CELL_TYPES)] = np.asarray(a.obsm["cell_type_proportions"], dtype=np.float32)
        comp[n_c:, len(CELL_TYPES)] = ecm_v

        obs = pd.DataFrame({
            "x": np.concatenate([a.obs["x"].values, e["x"].to_numpy(float)]),
            "y": np.concatenate([a.obs["y"].values, e["y"].to_numpy(float)]),
            "cell_type": pd.Categorical(
                np.concatenate([np.asarray(a.obs["cell_type"].values, dtype=object),
                                np.full(n_e, "ecm", dtype=object)]),
                categories=ECM_CELL_TYPES),
            "volume": np.concatenate([a.obs["volume"].values,
                                      np.full(n_e, ECM_VOXEL_UM ** 3)]),
            "is_ecm": np.concatenate([np.zeros(n_c, bool), np.ones(n_e, bool)]),
        }, index=pd.Index([str(i) for i in range(n_c + n_e)]))

        new = ad.AnnData(obs=obs)
        new.obsm["cell_type_proportions"] = pd.DataFrame(comp, columns=ECM_CELL_TYPES,
                                                         index=obs.index)
        for k, v in a.uns.items():
            new.uns[k] = v
        new.uns["ecm_source"] = str(ecm_csv_path(roi).relative_to(REPO))
        new.uns["ecm_divisor"] = ECM_DIVISOR
        new.uns["ecm_normalisation"] = "per-ROI 95th-percentile (see prep_imc_spatial/archive/scaleimage.py)"
        new.write_h5ad(out / "{}.h5ad".format(roi))
        rows.append({"roi": roi, "cells": n_c, "ecm_voxels": n_e,
                     "ecm_min": float(ecm_v.min()), "ecm_max": float(ecm_v.max()),
                     "ecm_mean": float(ecm_v.mean())})
        if verbose:
            print("  {:<12s} {:5d} cells + {:5d} ecm rows -> {}".format(roi, n_c, n_e, out.name))
    return pd.DataFrame(rows)


ECM_SWEEP = [
    {"name": "ecm_fov96",  "fov_um":  96, "pixel_size": 2.0, "global_size": 48, "local_size": 32},
    {"name": "ecm_fov128", "fov_um": 128, "pixel_size": 2.0, "global_size": 64, "local_size": 48},
    # Added 2026-08-22. Geometry is B_fov192_px2.0's verbatim (96/64), not a fresh guess,
    # so this arm stays comparable to the S1 round-1 192 um result.
    #
    # NOTE THE PRIOR THIS IS TESTING AGAINST: S1 round 1 found 128 um beat every 192 um
    # variant WITHOUT ECM. That was 4 ROIs from one patient on a different schedule, so it
    # is not decisive, and S4's reading -- that ECM supplies large-scale context a small
    # window lacks -- predicts ECM should help LESS at 192, not more. This run is the
    # direct test of that.
    {"name": "ecm_fov192", "fov_um": 192, "pixel_size": 2.0, "global_size": 96, "local_size": 64},
]

ECM_OVERRIDES = {
    "max_epochs": 200,
    "warm_up_epochs": 20,             # (0, 20, 180, 200) = 10/80/10, upstream's SHAPE
    "warm_down_epochs": 20,
    "checkpoint_interval_epochs": 25,
    "n_crops_for_tissue_train": 1024,
    # 20 flat, by decision.
    #
    # MEASURED AFTER THE FACT, and larger than the estimate this was chosen against:
    # ECM values are all < 1, which flips CropperSparseTensor to its one-hot branch
    # (dataset.py:365), so n_elements counts NONZERO ENTRIES, not cell mass. Matrix
    # alone contributes up to 23.5 entries per 96 um window and 41.3 per 128 um window
    # (20 um lattice, ~2.2e-3 voxels/um2). So at 128 um the threshold is satisfied by
    # matrix alone and every tile passes; at 96 um it needs only a handful of cells.
    # The ECM arm therefore trained on a strictly larger, emptier tile set than the
    # 10-channel arm at threshold 10 -- a real difference between the arms, not just
    # a channel difference. Analysis re-filters both to a shared grid to compensate.
    "n_element_min_for_crop": 20,
    "min_weight_decay": 0.0,
    "max_weight_decay": 0.0,
    "batch_size_per_gpu": 128,
    "image_in_ch": len(ECM_CELL_TYPES),
    "categories_to_channels": {t: i for i, t in enumerate(ECM_CELL_TYPES)},
}

# The 1000-epoch ECM test. Fresh runs, NOT resumed from the 200-epoch checkpoints: those
# annealed to min LR at epoch 199, and restarting an annealed model into a new schedule is
# not the same experiment as training one for 1000 epochs.
#
# warm_up/warm_down go to 100/100, which is the SAME 10/80/10 shape the 200-epoch screen
# used at 20/20 and, more importantly, exactly the final sweep's protocol -- so these runs
# are readable against final_fov96/final_fov128 at any matched epoch without the LR
# confound that S4 had to work around.
#
# n_element_min_for_crop stays flat at 20, by decision. Consequence, unchanged from S4:
# ECM values are < 1, which flips CropperSparseTensor to its one-hot branch, and matrix
# alone supplies ~23 entries per 96 um window and ~41 per 128 um one -- so this arm still
# trains on a larger, emptier tile set than the 10-channel arm does at threshold 10. A win
# here bundles "extra channel" with "different training tiles"; that is a known, accepted
# limitation of this test, not an oversight.
ECM_FINAL_RUN_DIR = RUNS_DIR / "ecm_final"
ECM_FINAL_OVERRIDES = dict(ECM_OVERRIDES,
                           max_epochs=1000, warm_up_epochs=100, warm_down_epochs=100)
ECM_FINAL_FOVS = (96, 128)          # 192 lost the FOV comparison in S4c; not carried forward


def ecm_final_specs():
    return [s for s in ECM_SWEEP if s["fov_um"] in ECM_FINAL_FOVS]


def ecm_final_config_dict(spec):
    cfg = load_config()
    cfg.update(ECM_FINAL_OVERRIDES)
    cfg["data_folder"] = str(ecm_h5ad_dir(spec["name"]))
    cfg["pixel_size"] = float(spec["pixel_size"])
    cfg["global_size"] = int(spec["global_size"])
    cfg["local_size"] = int(spec["local_size"])
    assert abs(cfg["pixel_size"] * cfg["global_size"] - spec["fov_um"]) < 1e-6
    return cfg


def ecm_final_feat_dir(name):
    """Featurised output for a 1000-epoch ECM run.

    Deliberately NOT feat_dir(name): the 200-epoch screen already owns
    imc_anndata_feat_ecm_fov96 at epoch 199, and run_featurize_screen rm -rf's its output
    folder before writing. Sharing the path would silently destroy the S4b inputs.
    """
    return REPO / "tissuemosaic" / "imc_anndata_feat_long_{}".format(name)


def ecm_final_checkpoint(name, epoch=None):
    return final_checkpoint(name, epoch=epoch, root=ECM_FINAL_RUN_DIR)


def ecm_final_status():
    rows = []
    for spec in ecm_final_specs():
        name = spec["name"]
        eps = final_periodic_epochs(name, root=ECM_FINAL_RUN_DIR)
        lc = read_loss_curve(name, root=ECM_FINAL_RUN_DIR)
        rows.append({"name": name, "fov_um": spec["fov_um"],
                     "last_epoch": int(lc["epoch"].max()) if len(lc) else 0,
                     "last_loss": float(lc["loss"].iloc[-1]) if len(lc) else float("nan"),
                     "periodic_ckpts": len(eps),
                     "latest_ckpt_epoch": eps[-1] if eps else None})
    return pd.DataFrame(rows).set_index("name")


ECM_RUN_DIR = RUNS_DIR / "ecm_screen"


def ecm_h5ad_dir(name):
    """One folder per config -- both share pixel_size 2.0 and would race on the cache."""
    return REPO / "tissuemosaic" / "imc_anndata_{}".format(name)


def ecm_config_dict(spec):
    cfg = load_config()
    cfg.update(ECM_OVERRIDES)
    cfg["data_folder"] = str(ecm_h5ad_dir(spec["name"]))
    cfg["pixel_size"] = float(spec["pixel_size"])
    cfg["global_size"] = int(spec["global_size"])
    cfg["local_size"] = int(spec["local_size"])
    assert abs(cfg["pixel_size"] * cfg["global_size"] - spec["fov_um"]) < 1e-6
    return cfg


def feat_dir(name):
    """Where main_2_featurize output for a run lives.

    The final-sweep folders predate the ECM screen and are named without the `final_`
    prefix (`imc_anndata_feat_fov96`); S3 reads those paths, so they are kept as they are
    rather than renamed.
    """
    return REPO / "tissuemosaic" / "imc_anndata_feat_{}".format(
        name[len("final_"):] if name.startswith("final_") else name)


def ecm_checkpoint(name, epoch=None):
    """Checkpoint for an ECM-screen run. Same layout as the final sweep, other root."""
    return final_checkpoint(name, epoch=epoch, root=ECM_RUN_DIR)


def ecm_periodic_epochs(name):
    return final_periodic_epochs(name, root=ECM_RUN_DIR)


def ecm_status():
    """Progress of the ECM screen, safe to call mid-run."""
    rows = []
    for spec in ECM_SWEEP:
        name = spec["name"]
        eps = ecm_periodic_epochs(name)
        lc = read_loss_curve(name, root=ECM_RUN_DIR)
        rows.append({"name": name, "fov_um": spec["fov_um"],
                     "last_epoch": int(lc["epoch"].max()) if len(lc) else 0,
                     "last_loss": float(lc["loss"].iloc[-1]) if len(lc) else float("nan"),
                     "periodic_ckpts": len(eps),
                     "latest_ckpt_epoch": eps[-1] if eps else None})
    return pd.DataFrame(rows).set_index("name")


def window_key(roi, centers_um):
    """Stable identity for a tiled window: ROI plus centre to 1 nm.

    The tiling origin is `torch.randint(low=-crop_size, high=0)` seeded per ROI and the
    stride follows from global_size alone, so two featurisations of the SAME ROI at the
    SAME fov produce the same candidate centres regardless of how many channels the
    raster has. That makes the centre a usable join key between an ECM run and its
    no-ECM counterpart -- which is the only way to compare their kappa honestly, since
    n_element_min culls a different subset in each.
    """
    c = np.asarray(centers_um, dtype=float)
    return np.array(["{}|{:.3f}|{:.3f}".format(r, x, y)
                     for r, x, y in zip(np.asarray(roi), c[:, 0], c[:, 1])])


def align_to_grid(reference, other):
    """Row indices of `other` matching every row of `reference`, by window_key.

    Raises if `reference` is not a subset of `other` -- silently dropping unmatched
    windows would compare two different grids while reporting one number.
    """
    kr, ko = window_key(reference["roi"], reference["centers_um"]), \
             window_key(other["roi"], other["centers_um"])
    pos = {k: i for i, k in enumerate(ko)}
    missing = [k for k in kr if k not in pos]
    if missing:
        raise KeyError("{} reference windows absent from the other grid, e.g. {}".format(
            len(missing), missing[:3]))
    return np.array([pos[k] for k in kr])


def raster_frame(adata, pixel_size, padding=10):
    """Cell x/y moved into the raster frame windows are reported in.

    SparseImage.from_anndata subtracts each ROI's own minimum and pads by 10 PIXELS
    (datamodule.py:640, hardcoded), so obs['x'] in ROI-centred microns and a window
    centre from read_from_patch_dictionary live in different frames -- on JHH369ROI1
    that is cells spanning x[-501, 497] against windows spanning x[-8, 1048]. Overlaying
    them without this transform silently draws the cells in the wrong place.
    """
    x = np.asarray(adata.obs["x"].values, dtype=float)
    y = np.asarray(adata.obs["y"].values, dtype=float)
    off = float(padding) * float(pixel_size)
    return (x - x.min()) + off, (y - y.min()) + off


def umap_graph(features, n_neighbors=15, pca_var=0.9, seed=0):
    """SmartPca -> SmartUmap, returning the graph so Leiden can be re-run cheaply.

    cluster_embedding_tm rebuilds PCA and UMAP on every call, which makes a resolution
    or seed sweep spend all its time recomputing an embedding that never changes.

    SEEDING IS REQUIRED, and the reason is upstream. SmartUmap does hardcode
    random_state=0, but SmartPca.fit regularises the covariance matrix with an UNSEEDED
    draw before taking its SVD (validation_util.py):

        eps = 1E-4 * torch.randn(p, dtype=cov.dtype, device=cov.device)
        cov += torch.diag(eps)

    Against 512-dimensional DINO features that jitter is negligible and results come out
    bit-identical. Against a 10-dimensional composition vector -- z-scored, so the
    covariance diagonal is ~1, with many near-duplicate rows -- it perturbs the retained
    subspace enough to flip kNN ties, and the downstream Leiden partition changes with it.
    MEASURED before this seed was added: kappa for the composition baseline at 96 um drew
    anywhere in 0.147-0.160 across rebuilds of identical data, a 4% wobble on the
    denominator of every gain in S4.
    """
    import torch
    from tissuemosaic.utils import SmartPca, SmartUmap

    torch.manual_seed(seed)
    X = np.asarray(features, dtype=np.float32)
    pca = SmartPca(preprocess_strategy="z_score")
    emb = pca.fit_transform(torch.tensor(X), n_components=pca_var)
    um = SmartUmap(n_neighbors=min(n_neighbors, len(X) - 1), preprocess_strategy="raw",
                   n_components=2, min_dist=0.5, metric="euclidean")
    u = um.fit_transform(emb)
    return {"graph": um.get_graph(), "umap": to_numpy(u), "pca": to_numpy(emb)}


def leiden_on(graph, resolution=1.0, seed=0):
    from tissuemosaic.utils import SmartLeiden
    return np.asarray(SmartLeiden(graph=graph, directed=True).cluster(
        resolution=resolution, random_state=seed, partition_type="RBC")).astype(int)


def final_run_dir(name, root=None):
    return (Path(root) if root else FINAL_RUN_DIR) / name


def final_checkpoint(name, epoch=None, root=None):
    """The final checkpoint, or a periodic one.

    `ckpt_last.pt` only appears once a run finishes all 1000 epochs. Periodic
    checkpoints land every `checkpoint_interval_epochs` under .neptune/, and are
    fully usable mid-run -- which is how a 26 h run can be evaluated at epoch 500.
    """
    d = final_run_dir(name, root)
    if epoch is None:
        last = d / "ckpt_last.pt"
        if last.is_file():
            return last
        hits = sorted(
            d.glob(".neptune/**/checkpoints/periodic_checkpoint-epoch=*.ckpt"),
            key=lambda f: int(re.search(r"epoch=(\d+)", f.name).group(1)),
        )
        if not hits:
            raise FileNotFoundError(f"no checkpoint yet for {name}")
        return hits[-1]
    hits = sorted(
        d.glob(f".neptune/**/checkpoints/periodic_checkpoint-epoch={epoch}.ckpt")
    )
    if not hits:
        raise FileNotFoundError(
            f"no periodic checkpoint at epoch {epoch} for {name}; have {final_periodic_epochs(name, root)}"
        )
    return hits[0]


def final_periodic_epochs(name, root=None):
    """Which periodic checkpoint epochs exist for a final-sweep config."""
    d = final_run_dir(name, root)
    return sorted(
        int(re.search(r"epoch=(\d+)", f.name).group(1))
        for f in d.glob(".neptune/**/checkpoints/periodic_checkpoint-epoch=*.ckpt")
    )


def heldout_rois():
    """Cohort ROIs NOT in the final training set -- the real transfer test.

    A JHH369-only model scored gain 0.33-0.78x on unseen patients, i.e. worse than
    clustering cell counts. These ROIs are what says whether five patients fixed it.
    """
    everything = sorted(f[:-5] for f in os.listdir(H5AD_DIR) if f.endswith(".h5ad"))
    return [r for r in everything if r not in set(FINAL_ROIS)]


def final_status(root=None):
    """Progress of the final sweep, safe to call mid-run."""
    rows = []
    for spec in FINAL_SWEEP:
        name = spec["name"]
        eps = final_periodic_epochs(name, root)
        d = final_run_dir(name, root)
        lc = read_loss_curve(name, root=root or FINAL_RUN_DIR)
        rows.append(
            {
                "name": name,
                "fov_um": spec["fov_um"],
                "done": (d / "ckpt_last.pt").is_file(),
                "last_epoch": int(lc["epoch"].max()) if len(lc) else 0,
                "last_loss": float(lc["loss"].iloc[-1]) if len(lc) else float("nan"),
                "periodic_ckpts": len(eps),
                "latest_ckpt_epoch": eps[-1] if eps else None,
            }
        )
    return pd.DataFrame(rows).set_index("name")


def final_config_dict(spec):
    """Full config for one final-sweep entry."""
    cfg = load_config()
    cfg.update(FINAL_OVERRIDES)
    cfg["data_folder"] = str(final_h5ad_dir(spec["name"]))
    cfg["pixel_size"] = float(spec["pixel_size"])
    cfg["global_size"] = int(spec["global_size"])
    cfg["local_size"] = int(spec["local_size"])
    fov = cfg["pixel_size"] * cfg["global_size"]
    assert abs(fov - spec["fov_um"]) < 1e-6, "pixel_size x global_size != declared fov"
    return cfg


SWEEP_OVERRIDES = {
    # === Option 1 (probe P4), selected 2026-08-20 ===================================
    # Chosen over option 2 (thr 50 / warm_down 20) and option 3 (v1 proxy) on a
    # three-arm comparison at config-A geometry with IDENTICAL analysis (tiling, no
    # overlap, seed 0, 50-cell window filter), so only training differed:
    #
    #   option        final loss  clusters  coherence_z  border clusters >=80%
    #   1 (P4)          0.865        11        23.52            0
    #   2 (P5)          0.813        10        19.31            0
    #   3 (v1 proxy)    0.756         9        14.83            1  (82.4% border)
    #
    # coherence_z is a fair comparator HERE, unlike across the FOV sweep: all three
    # arms share config-A geometry and the same filter, so N is matched and z's
    # ~sqrt(N) scaling cancels. Rank the FOV sweep on coherence_kappa instead.
    #
    # +21.8% coherence over option 2, +58.6% over option 3. Option 3 still produces a
    # border cluster even WITH the analysis filter -- v1's training regime is genuinely
    # sensitive to partial-tissue windows, which filtering alone does not fix.
    "max_epochs": 200,
    "warm_up_epochs": 20,
    # NOT upstream's shape, despite sharing the number. TissueMosaic ships
    # warm_down 100 with max_epochs 1000 -- a 10% decay tail. At max_epochs 200 the
    # same 100 is a 50% tail. The shape-matching choice was warm_down 20 (probe P5);
    # it lost on measurement, so this deviates from upstream deliberately.
    "warm_down_epochs": 100,
    "checkpoint_interval_epochs": 25,
    # Power of two: ncpt=8 -> batch_size_dataloader=16, so 16 tissue draws per batch
    # AND 32 steps/epoch (matching v1's step count). Exactly 128 crops per step.
    "n_crops_for_tissue_train": 1024,
    # Also MATCHES upstream: config_dino_ssl_testis.yaml sets 10, while the barlow,
    # simclr and vae configs set 200 -- so 10 is TissueMosaic's own DINO choice.
    # TRAINING threshold stays permissive. The border problem is handled at ANALYSIS
    # time by fov_scaled_threshold (see S1 / featurize_windows); raising it in training
    # was tested as option 2 and did not help. It also roughly doubles cropping cost
    # when ncpt is small, which is precisely this batching regime.
    # NB this value is also what cropper_test uses (datamodule.py:331), so it is the
    # floor under the analysis filter too -- harmless, since the scaled thresholds
    # (13/29/52/117 for 64/96/128/192 um) all sit above it.
    "n_element_min_for_crop": 10,
    # This MATCHES upstream, which was not obvious. The argparse defaults are
    # 0.04/0.4 (dino.py:430), but run/config_dino_ssl_testis.yaml ships 0.0/0.0 --
    # and the yaml wins (main_1_train_ssl.py resolves yaml > CLI > defaults). So
    # zeroing these restores TissueMosaic's own setting rather than departing from
    # it. Independently confirmed here: on 4 ROIs / 30k cells the argparse values
    # raised the loss floor from ~0.76 to 2.82.
    "min_weight_decay": 0.0,
    "max_weight_decay": 0.0,
}


def sweep_spec(name):
    for s in WINDOW_SWEEP:
        if s["name"] == name:
            return s
    raise KeyError(
        "unknown sweep config {!r}; have {}".format(
            name, [s["name"] for s in WINDOW_SWEEP]
        )
    )


def sweep_subset_folder(patient=None, dest=None, force=False):
    """Folder holding only one patient's h5ad, as AnndataFolderDM consumes a directory.

    Uses symlinks, so it costs nothing and cannot drift from imc_anndata/.
    Any stale train/test_dataset.pt is cleared: those caches are built for one
    specific pixel_size and MUST NOT be reused across sweep configs.
    """
    patient = patient or SWEEP_PATIENT
    dest = Path(dest) if dest else SWEEP_H5AD_DIR
    dest.mkdir(parents=True, exist_ok=True)

    for stale in list(dest.glob("*.pt")) + (
        [p for p in dest.glob("*.h5ad")] if force else []
    ):
        stale.unlink()

    src = sorted(H5AD_DIR.glob(f"{patient}*.h5ad"))
    if not src:
        raise RuntimeError(
            f"no h5ad for patient {patient} in {H5AD_DIR} -- run make_imc_anndata.py first"
        )
    linked = []
    for f in src:
        link = dest / f.name
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(f)
        linked.append(link.name)
    return dest, linked


def sweep_config_dict(spec, data_folder=None):
    """Full config for one sweep entry: the production yaml plus the deltas."""
    cfg = load_config()
    cfg.update(SWEEP_OVERRIDES)
    cfg["data_folder"] = str(data_folder or SWEEP_H5AD_DIR)
    cfg["pixel_size"] = float(spec["pixel_size"])
    cfg["global_size"] = int(spec["global_size"])
    cfg["local_size"] = int(spec["local_size"])
    fov = cfg["pixel_size"] * cfg["global_size"]
    assert abs(fov - spec["fov_um"]) < 1e-6, (
        "pixel_size x global_size = {} != declared fov {}".format(fov, spec["fov_um"])
    )
    return cfg


def sweep_run_dir(name, root=None):
    return (Path(root) if root else SWEEP_RUN_DIR) / name


def sweep_checkpoint(name, root=None):
    return sweep_run_dir(name, root) / "ckpt_last.pt"


def sweep_status(root=None):
    """Which sweep configs have finished training."""
    rows = []
    for spec in WINDOW_SWEEP:
        d = sweep_run_dir(spec["name"], root)
        ck = sweep_checkpoint(spec["name"], root)
        n_periodic = (
            len(list(d.glob(".neptune/**/checkpoints/periodic*.ckpt")))
            if d.is_dir()
            else 0
        )
        rows.append(
            {
                "name": spec["name"],
                "fov_um": spec["fov_um"],
                "pixel_size": spec["pixel_size"],
                "global_size": spec["global_size"],
                "local_size": spec["local_size"],
                "trained": ck.is_file(),
                "ckpt_MB": round(ck.stat().st_size / 1e6) if ck.is_file() else 0,
                "periodic_ckpts": n_periodic,
            }
        )
    return pd.DataFrame(rows)


def read_loss_curve(name, root=None):
    """Per-epoch training loss parsed out of a sweep run's train.log.

    Lightning writes the progress bar with carriage returns, so the log is one
    enormous line until \r is translated.
    """
    log = sweep_run_dir(name, root) / "train.log"
    if not log.is_file():
        return pd.DataFrame(columns=["epoch", "loss"])
    pat = re.compile(r"Epoch (\d+):.*?loss=([0-9.]+)")
    seen = {}
    for line in log.read_text(errors="ignore").replace("\r", "\n").splitlines():
        for ep, loss in pat.findall(line):
            seen[int(ep)] = float(loss)  # keep the last value reported for each epoch
    return pd.DataFrame(sorted(seen.items()), columns=["epoch", "loss"])


# ----------------------------------------------- window featurisation -----


def featurize_windows(
    model, dm, adatas, frac_overlap=0.0, feature_key="dino", seed=0, verbose=True
):
    """Tile every ROI into windows and embed each one.

    Tiling (not random sampling) so each window has a known position, which the
    spatial-coherence metric needs. frac_overlap defaults to 0.0: overlapping
    windows share cells, which would inflate that metric.

    SEEDING IS REQUIRED, not optional. CropperSparseTensor's tiling branch draws a
    RANDOM grid origin --

        i0 = torch.randint(low=-crop_size, high=0, size=[1]).item()   # dataset.py:325
        j0 = torch.randint(low=-crop_size, high=0, size=[1]).item()
        for i in range(i0, w_img, stride): ...

    (note this is the SPARSE cropper; CropperDenseTensor at dataset.py:180 draws
    low=0, high=stride and stops at w_img - crop_size, so it emits no partial
    tiles. The sparse one starts off-image and runs past the far edge, which is
    where border windows come from -- n_element_min then culls the emptiest.)

    -- which is sensible as training augmentation but makes analysis irreproducible:
    unseeded, three identical calls here returned 356 / 345 / 341 windows with
    different contents. `seed` is re-applied per ROI (seed + index) so that adding
    or removing an ROI does not shift the tiling of the others.

    NOTE: TissueMosaic does not batch the tiling path (`## TODO: allow batching
    with tiling strategy`), so every tile of an ROI is pushed through the model
    at once. Featurising one ROI at a time keeps that bounded.

    Returns a dict of stacked arrays, one row per window.
    """
    import torch
    from tissuemosaic.data.dataset import CropperSparseTensor
    from tissuemosaic.models.patch_analyzer import Composition, SpatialAutocorrelation

    use_gpu = torch.cuda.is_available()
    if use_gpu:
        model = model.cuda()

    feats, comps, morans, rois, centers, ncells, necms, nvox = ([] for _ in range(8))
    for i_roi, (roi, adata) in enumerate(sorted(adatas.items())):
        torch.manual_seed(seed + i_roi)  # pin the random tiling origin -- see docstring
        sp = dm.anndata_to_sparseimage(adata)
        if use_gpu:
            sp = sp.cuda()
        sp.compute_patch_features(
            feature_name=feature_key,
            datamodule=dm,
            model=model,
            strategy="tiling",
            fraction_patch_overlap=frac_overlap,
            overwrite=True,
        )
        values, xywh = sp.read_from_patch_dictionary(key=feature_key)
        xywh_np = to_numpy(xywh).astype(float)

        tensors, _, _ = CropperSparseTensor.reapply_crops(sp.data, xywh)
        comp = (
            torch.stack(Composition(return_fraction=True)(tensors), dim=0).cpu().numpy()
        )
        # absolute occupancy: `comp` is FRACTIONS and sums to 1, so it cannot supply
        # this. Needed by the analysis filter and by any density-vs-scale question.
        #
        # Split by channel, not summed over all of them: on 11-channel ECM data
        # channel 10 carries a CONTINUOUS field on a 20 um lattice, so a whole-tensor
        # sum mixes 20-40 matrix voxels into what callers read as a cell count.
        n_ch = int(sp.data.shape[-3])
        n_cell_ch = len(CELL_TYPES)
        nc, ne = np.zeros(len(tensors)), np.zeros(len(tensors))
        # Entry COUNT for the ECM channels as well as their mass. Mass alone conflates
        # "how much matrix" with "how much of this window is on-tissue": a border window
        # holds fewer 20 um voxels AND less total matrix, so the two correlate at r=+0.77
        # and a cluster split on window occupancy scores as a split on matrix. The
        # geometry-free quantity is mass / voxels -- see ecm_density().
        nv = np.zeros(len(tensors))
        for i, t in enumerate(tensors):
            idx = t.indices()[0].cpu()
            cs = torch.zeros(n_ch).index_add_(0, idx, t.values().cpu().float())
            ct = torch.zeros(n_ch).index_add_(0, idx, torch.ones(idx.numel()))
            nc[i] = float(cs[:n_cell_ch].sum())
            ne[i] = float(cs[n_cell_ch:].sum())
            nv[i] = float(ct[n_cell_ch:].sum())
        # Moran's I needs more points than neighbours: sklearn raises
        # "Expected n_neighbors < n_samples_fit" on any window with <= 6 cells, which
        # a permissive n_element_min_for_crop will happily produce. Compute it only
        # where it is defined and leave the rest NaN, so a low threshold does not
        # crash featurisation -- the model forward pass has no such limit.
        n_nb = 6
        big = [i for i, t in enumerate(tensors) if t._nnz() > n_nb]
        moran = np.full((len(tensors), n_ch), np.nan, dtype=float)
        if big:
            m = (
                torch.stack(
                    SpatialAutocorrelation(
                        modality="moran", n_neighbours=n_nb, neigh_correct=True
                    )([tensors[i] for i in big]),
                    dim=0,
                )
                .cpu()
                .numpy()
            )
            moran[big] = m

        # window centres in um, relative to this ROI's raster origin. Only used
        # for within-ROI adjacency, so the origin offset is irrelevant.
        cx = (xywh_np[:, 0] + xywh_np[:, 2] / 2.0) * dm._pixel_size
        cy = (xywh_np[:, 1] + xywh_np[:, 3] / 2.0) * dm._pixel_size

        feats.append(to_numpy(values))
        comps.append(comp)
        morans.append(moran.max(axis=-1))
        centers.append(np.c_[cx, cy])
        ncells.append(nc)
        necms.append(ne)
        nvox.append(nv)
        rois += [roi] * len(values)
        if verbose:
            print(f"  {roi:<14s} {len(values):4d} windows")
        del sp

    return {
        "features": np.concatenate(feats, axis=0),
        "composition": np.concatenate(comps, axis=0),
        "moran": np.concatenate(morans, axis=0),
        "centers_um": np.concatenate(centers, axis=0),
        "n_cells": np.concatenate(ncells, axis=0),
        "n_ecm": np.concatenate(necms, axis=0),
        "n_ecm_voxels": np.concatenate(nvox, axis=0),
        "roi": np.asarray(rois),
    }


# --------------------------------------------------------- clustering -----


def ecm_density(result):
    """Mean matrix value per ECM voxel: the geometry-free measure of "how much matrix".

    Three quantities get called "ECM content" and only this one means it:

        n_ecm         sum of channel-10 values in the window. Scales with how many
                      20 um voxels the window contains, so a half-off-tissue border
                      window reads as low-matrix regardless of its actual density
                      (r = +0.77 with the voxel count alone).
        frac_ecm      n_ecm / total window mass. A RATIO against cell mass, so it moves
                      when cell density moves and is high wherever cells are sparse.
        ecm_density   n_ecm / n_ecm_voxels. Divides the geometry out; this is the mean
                      normalised matrix value over the tissue the window covers.

    MEASURED, why it matters: clustering eta^2 on n_ecm was 0.594 at 128 um, but eta^2 on
    the VOXEL COUNT alone -- pure occupancy, carrying no matrix values at all -- was 0.670.
    On this measure the same clusters score 0.322. Separation on matrix is real but
    substantially smaller than the mass-based number implies.

    NaN where a window contains no ECM voxels at all.
    """
    m = np.asarray(result["n_ecm"], dtype=float)
    n = np.asarray(result["n_ecm_voxels"], dtype=float)
    return np.where(n > 0, m / np.where(n > 0, n, 1.0), np.nan)


def merge_featurized(folder, feature_key="dino", rois=None):
    """Merge per-ROI featurized h5ad into ONE window-level AnnData for scanpy.

    main_2_featurize writes one file per ROI, with the window data buried in
    `uns['sparse_image_state_dict']['patch_properties_dict']`. anndata.concat would
    merge the CELLS and drop that, which is the wrong axis: the objects we cluster and
    plot are windows, not cells.

    So this pivots the axis -- one observation per WINDOW -- which is what makes the
    whole scanpy plotting API apply directly:

        X                    (n_windows, 512)  the embeddings
        obs['roi','patient'] provenance
        obs['x'],obs['y']    window centre, RASTER frame (see below)
        obs['n_cells']       occupancy
        obsm['spatial']      window centres -> sc.pl.embedding(basis='spatial')
        obsm['composition']  the library's own patch_ncv, as fractions
        var                  one row per feature dimension

    COORDINATE FRAME: patch_xywh is in raster pixels from the PADDED origin, while
    obs['x'/'y'] in the source file are microns centred on the ROI. Window centres here
    are microns in the RASTER frame (pixels x pixel_size). Use merge_cells() for the
    matching cell-level object -- it applies the same transform, so the two overlay.
    """
    import anndata as _ad
    from anndata import read_h5ad

    blocks = []
    for f in sorted(Path(folder).glob("*.h5ad")):
        if rois is not None and f.stem not in set(rois):
            continue
        a = read_h5ad(f)
        sd = a.uns["sparse_image_state_dict"]
        ppd = sd["patch_properties_dict"]
        px = float(sd["pixel_size"])
        feats = np.asarray(ppd[feature_key])
        xywh = np.asarray(ppd[feature_key + "_patch_xywh"], dtype=float)
        ncv = np.asarray(ppd["patch_ncv"], dtype=float)
        ncv = ncv / np.clip(ncv.sum(axis=1, keepdims=True), 1e-9, None)
        cx = (xywh[:, 0] + xywh[:, 2] / 2) * px
        cy = (xywh[:, 1] + xywh[:, 3] / 2) * px

        w = _ad.AnnData(X=np.asarray(feats, dtype=np.float32))
        w.obs_names = ["{}_w{:04d}".format(f.stem, i) for i in range(len(feats))]
        w.var_names = ["{}{:03d}".format(feature_key, i) for i in range(feats.shape[1])]
        w.obs["roi"] = f.stem
        w.obs["patient"] = str(a.uns.get("patient", roi_to_patient(f.stem)))
        w.obs["status"] = int(a.uns.get("status", -1))
        w.obs["x"] = cx
        w.obs["y"] = cy
        w.obs["fov_um"] = float(xywh[0, 2] * px)
        w.obsm["spatial"] = np.c_[cx, cy]
        w.obsm["patch_xywh"] = xywh.astype(int)   # RASTER PIXELS -- see window_mass()
        w.obsm["composition"] = ncv
        chans = channels_for(ncv.shape[1])          # 10, or 11 for an ECM run
        w.uns["channels"] = chans
        for j, t in enumerate(chans):               # also as obs cols, for sc.pl colouring
            w.obs["frac_" + t] = ncv[:, j]
        blocks.append(w)

    if not blocks:
        raise FileNotFoundError("no featurized h5ad in {}".format(folder))
    merged = _ad.concat(blocks, axis=0, join="outer", merge="same", index_unique=None)
    for c in ("roi", "patient"):
        merged.obs[c] = merged.obs[c].astype("category")
    # concat() drops uns unless uns_merge is given, so the per-block value set above does
    # not survive; re-attach it here.
    merged.uns["channels"] = list(chans)
    merged.uns["feature_key"] = feature_key
    merged.uns["source_folder"] = str(folder)
    return merged


def window_mass(windows, cells, pixel_size=2.0, padding=10):
    """Per-window, per-channel MASS and ENTRY COUNT, assigned pixel-exactly.

    WHY THIS EXISTS. `patch_ncv`, which merge_featurized exposes as obsm['composition'],
    is not a mass fraction. sparse_image.py:1026 builds it by COUNTING NONZERO ENTRIES
    per channel:

        cells = crop.indices()[[0]]
        ct_counts[code] = torch.count_nonzero(cells == code)
        ncv = ct_counts / ct_counts.sum()

    For a one-hot cell channel a count is a mass, so for channels 0-9 the two agree to
    within the odd two-cells-in-one-pixel collision. For a CONTINUOUS channel they do not.
    MEASURED on JHH369ROI1's densest 128 um window: patch_ncv puts ecm at 0.205 while the
    actual matrix share of the window's mass is 0.017 -- a 12x difference. The count-based
    number is essentially `n_ecm_voxels / (n_cells + n_ecm_voxels)`, and since ECM voxels
    sit on a fixed 20 um lattice at constant density, that is an inverse CELL-DENSITY
    measure carrying almost no information about how much matrix is present.

    Worse, it is not confined to the ECM column: those voxels are ~20% of the entries in a
    typical window, so every cell-type fraction in patch_ncv is diluted by a pure geometry
    term. Use the mass basis for anything about the ECM screen's composition.

    Membership is tested in RASTER PIXELS against the integer patch_xywh, not in microns
    against the window centre. The cropper's own test is `x_pixel >= x0 and
    x_pixel < x0 + w` (dataset.py:361); a micron box round-trips through a floor() and
    disagrees at the boundary by up to one pixel. On dense cells that is a ~6% wobble, but
    ECM sits on a 20 um lattice where one pixel of slop flips a whole row of voxels in or
    out -- which is what a micron-space check reports as a 0.09-0.19 "misalignment".

    Returns (mass, count), both (n_windows, n_channels), unnormalised.
    """
    xywh = np.asarray(windows.obsm["patch_xywh"], dtype=int)
    X = np.asarray(cells.X, dtype=float)
    w_roi = np.asarray(windows.obs["roi"]).astype(str)
    c_roi = np.asarray(cells.obs["roi"]).astype(str)
    xp = np.floor(np.asarray(cells.obsm["spatial"])[:, 0] / float(pixel_size)).astype(int)
    yp = np.floor(np.asarray(cells.obsm["spatial"])[:, 1] / float(pixel_size)).astype(int)

    mass = np.zeros((windows.n_obs, X.shape[1]))
    count = np.zeros_like(mass)
    for r in np.unique(w_roi):
        wi = np.flatnonzero(w_roi == r)
        ci = np.flatnonzero(c_roi == r)
        rx, ry, rX = xp[ci], yp[ci], X[ci]
        for i in wi:
            x0, y0, ww, hh = xywh[i]
            m = (rx >= x0) & (rx < x0 + ww) & (ry >= y0) & (ry < y0 + hh)
            if m.any():
                sub = rX[m]
                mass[i] = sub.sum(axis=0)
                count[i] = (sub > 0).sum(axis=0)
    return mass, count


# --------------------------------------------- arrangement vs composition -----
#
# "Do the motifs carry more than cell-type composition?" needs a control stronger than
# composition alone. A fraction vector carries neither cell DENSITY nor how much of the
# window is on tissue, and the embedding demonstrably encodes both -- clustering eta^2 on
# the ECM voxel count alone (pure on-tissue area, no biology at all) reached 0.65. A model
# that only knew "half off the edge, 300 cells" would beat a composition baseline and look
# like it had found motifs. So CONTROL_FEATURES below is composition + density + occupancy
# + matrix, and the claim is only ever "beyond all of that".

def window_cell_index(windows, cells, pixel_size=2.0, exclude_ecm=True):
    """Row indices of `cells` falling inside each window, tested in RASTER PIXELS.

    Same membership rule as window_mass -- integer patch_xywh, not a micron box round the
    centre -- for the same reason: a micron box disagrees at the boundary by up to a pixel.
    """
    xywh = np.asarray(windows.obsm["patch_xywh"], dtype=int)
    w_roi = np.asarray(windows.obs["roi"]).astype(str)
    c_roi = np.asarray(cells.obs["roi"]).astype(str)
    keep = ~np.asarray(cells.obs["is_ecm"], dtype=bool) if exclude_ecm else np.ones(cells.n_obs, bool)
    sp = np.asarray(cells.obsm["spatial"])
    xp = np.floor(sp[:, 0] / float(pixel_size)).astype(int)
    yp = np.floor(sp[:, 1] / float(pixel_size)).astype(int)

    out = [np.empty(0, dtype=int)] * windows.n_obs
    for r in np.unique(w_roi):
        ci = np.flatnonzero((c_roi == r) & keep)
        rx, ry = xp[ci], yp[ci]
        for i in np.flatnonzero(w_roi == r):
            x0, y0, ww, hh = xywh[i]
            m = (rx >= x0) & (rx < x0 + ww) & (ry >= y0) & (ry < y0 + hh)
            out[i] = ci[m]
    return out


# Ordered pairs worth naming for this project: antigen presentation is a question about
# who is NEXT TO whom, and every one of these is invisible to a composition vector -- two
# windows can hold identical CD8 and tumour counts with the T cells either infiltrating or
# excluded at the margin.
ARRANGEMENT_PAIRS = [
    ("apCAF", "CD8_Tcell"),
    ("CD8_Tcell", "apCAF"),
    ("CD8_Tcell", "epithelial_tumor_class1"),
    ("CD4_Tcell", "apCAF"),
    ("CAF", "CD8_Tcell"),
]


def arrangement_targets(windows, cells, k=6, min_cells=25, min_type=4,
                        pairs=None, pixel_size=2.0):
    """Per-window statistics of HOW cells are laid out, independent of how many there are.

    Every target is composition-normalised by construction: it is an observed neighbour
    fraction MINUS the window-level fraction of that type. A type being abundant therefore
    cannot make its score high -- only its being spatially clumped (positive) or dispersed
    (negative) can. That is what makes these fair targets for asking whether an embedding
    knows something a composition vector does not.

        homo_<type>   how much type t's neighbours are type t, above chance for this window
        mixing        fraction of kNN edges joining DIFFERENT types; low = segregated
        enr_A__B      how much of A's neighbourhood is B, above B's share of the window

    NaN where undefined: fewer than `min_cells` cells in the window, or fewer than
    `min_type` cells of the type in question. NaN is propagated rather than filled -- a
    window with two CD8 cells has no measurable CD8 neighbourhood, and zero-filling it
    would invent a "not clumped" observation.
    """
    from scipy.spatial import cKDTree

    pairs = ARRANGEMENT_PAIRS if pairs is None else pairs
    idx_lists = window_cell_index(windows, cells, pixel_size=pixel_size)
    types = np.asarray(cells.obs["cell_type"]).astype(object)
    xy = np.asarray(cells.obsm["spatial"], dtype=float)

    cols = (["homo_" + t for t in CELL_TYPES] + ["mixing"]
            + ["enr_{}__{}".format(a, b) for a, b in pairs])
    out = np.full((windows.n_obs, len(cols)), np.nan)

    for i, ci in enumerate(idx_lists):
        n = len(ci)
        if n < max(min_cells, k + 1):
            continue
        t_i, p_i = types[ci], xy[ci]
        kk = min(k, n - 1)
        _, nb = cKDTree(p_i).query(p_i, k=kk + 1)
        nb = nb[:, 1:]                                   # drop self
        nb_types = t_i[nb]                               # (n, kk)
        share = {t: float((t_i == t).mean()) for t in set(t_i.tolist())}

        for j, t in enumerate(CELL_TYPES):
            m = t_i == t
            if m.sum() < min_type:
                continue
            out[i, j] = (nb_types[m] == t).mean() - share.get(t, 0.0)
        out[i, len(CELL_TYPES)] = (nb_types != t_i[:, None]).mean()
        for j, (a, b) in enumerate(pairs):
            m = t_i == a
            if m.sum() < min_type or share.get(b, 0.0) == 0.0:
                continue
            out[i, len(CELL_TYPES) + 1 + j] = (nb_types[m] == b).mean() - share[b]

    return pd.DataFrame(out, columns=cols)


def control_features(windows, mass, count, n_cell_ch=None):
    """The "just composition" model: composition + density + occupancy + matrix.

    Deliberately NOT composition alone -- see the block comment above.
    """
    n_cell_ch = len(CELL_TYPES) if n_cell_ch is None else n_cell_ch
    comp = mass[:, :n_cell_ch]
    comp = comp / np.clip(comp.sum(1, keepdims=True), 1e-9, None)
    n_cells = mass[:, :n_cell_ch].sum(1)
    cols = {t: comp[:, j] for j, t in enumerate(CELL_TYPES)}
    cols["log_n_cells"] = np.log1p(n_cells)
    if mass.shape[1] > n_cell_ch:                       # 11-channel data: matrix too
        vox = count[:, n_cell_ch:].sum(1)
        cols["ecm_voxels"] = vox
        cols["ecm_density"] = np.where(vox > 0, mass[:, n_cell_ch:].sum(1) / np.clip(vox, 1, None), 0.0)
    else:
        cols["occupancy"] = count[:, :n_cell_ch].sum(1)
    return pd.DataFrame(cols)


def spatial_blocks(centers_um, roi, n_split=2):
    """Fold labels that hold out a contiguous SPATIAL BLOCK from every ROI at once.

    Two different questions hide behind "cross-validated R^2" here, and they need different
    folds:

      hold out ROIs   -- does this information TRANSFER to an unseen section/patient?
      hold out blocks -- is the information PRESENT at all?

    Only the second is what "the motifs carry more than composition" claims. Holding out
    whole ROIs also demands cross-patient generalisation, and on 8 ROIs spanning 2.6k-10.6k
    cells/mm2 the distribution shift dominates: MEASURED, the control model scored R^2 0.58
    on `mixing` but 0.003 on apCAF-CD8 enrichment, so most targets had no headroom left for
    the embedding to add anything to.

    Blocks still defeat the leak that makes random folds useless -- adjacent windows share
    cells, so a window's neighbour in the training fold makes it trivially predictable --
    because a held-out block only touches training windows along its border.

    Splits each ROI into n_split^2 quantile tiles and pools tile k across ROIs into fold k.
    """
    centers_um = np.asarray(centers_um, dtype=float)
    roi = np.asarray(roi).astype(str)
    out = np.zeros(len(roi), dtype=int)
    for r in np.unique(roi):
        m = np.flatnonzero(roi == r)
        bx = np.searchsorted(np.quantile(centers_um[m, 0], np.linspace(0, 1, n_split + 1)[1:-1]),
                             centers_um[m, 0])
        by = np.searchsorted(np.quantile(centers_um[m, 1], np.linspace(0, 1, n_split + 1)[1:-1]),
                             centers_um[m, 1])
        out[m] = bx * n_split + by
    return out


RIDGE_ALPHAS = np.logspace(-1, 5, 13)


def _grouped_r2(X, y, groups, alphas=None):
    """Out-of-fold R^2, folds held out by GROUP (ROI), ridge penalty tuned inside each fold.

    Random folds leak badly here: adjacent windows share cells and tissue, so a window's
    neighbour in the training fold makes it trivially predictable. Whole ROIs are held out.

    THE PENALTY MUST BE TUNED, not fixed. Held-out ROIs are genuinely out of distribution --
    different patient, 2.6k to 10.6k cells/mm2 -- so an under-regularised model extrapolates
    catastrophically. MEASURED with alpha fixed at 1.0: control (13 features) reached
    R2 = 0.56 while control + 16 embedding PCs collapsed to R2 = -1.76 on the same folds,
    which says nothing about the embedding and everything about the penalty. RidgeCV picks
    alpha by efficient leave-one-out on the TRAINING fold only, so no test ROI informs it.
    """
    from sklearn.linear_model import RidgeCV
    from sklearn.model_selection import LeaveOneGroupOut
    from sklearn.metrics import r2_score
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    alphas = RIDGE_ALPHAS if alphas is None else alphas
    X, y, groups = np.asarray(X, float), np.asarray(y, float), np.asarray(groups)
    pred = np.full(len(y), np.nan)
    for tr, te in LeaveOneGroupOut().split(X, y, groups):
        m = make_pipeline(StandardScaler(), RidgeCV(alphas=alphas)).fit(X[tr], y[tr])
        pred[te] = m.predict(X[te])
    return float(r2_score(y, pred))


def nested_delta_r2(control, embed, y, groups, n_perm=200, seed=0, alphas=None):
    """How much arrangement information the embedding adds OVER the control.

    Nested, not head-to-head: the embedding obviously encodes composition, and comparing
    the two directly would just re-measure that. R2(control + embedding) - R2(control)
    asks only about the part composition cannot reach.

    THE PERMUTATION NULL IS NOT OPTIONAL. A 512-d predictor beats a 13-d one on almost
    anything through flexibility alone, so a positive delta means nothing by itself. The
    null shuffles embedding rows WITHIN each ROI -- preserving every ROI-level effect and
    destroying only the window-to-window correspondence -- and regenerates delta from pure
    capacity. Only the excess over that null is evidence. Same alpha grid and same
    dimensionality throughout, so the comparison is like-for-like.
    """
    control, embed, y = np.asarray(control, float), np.asarray(embed, float), np.asarray(y, float)
    groups = np.asarray(groups)
    ok = np.isfinite(y) & np.isfinite(control).all(1) & np.isfinite(embed).all(1)
    if ok.sum() < 50 or len(np.unique(groups[ok])) < 3:
        return dict(n=int(ok.sum()), r2_control=np.nan, r2_full=np.nan, delta=np.nan,
                    null_mean=np.nan, null_p95=np.nan, p=np.nan)
    C, E, Y, G = control[ok], embed[ok], y[ok], groups[ok]

    r2_c = _grouped_r2(C, Y, G, alphas)
    r2_f = _grouped_r2(np.c_[C, E], Y, G, alphas)
    delta = r2_f - r2_c

    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for b in range(n_perm):
        Ep = E.copy()
        for g in np.unique(G):                          # shuffle WITHIN ROI
            m = np.flatnonzero(G == g)
            Ep[m] = E[rng.permutation(m)]
        null[b] = _grouped_r2(np.c_[C, Ep], Y, G, alphas) - r2_c
    return dict(n=int(ok.sum()), r2_control=r2_c, r2_full=r2_f, delta=delta,
                null_mean=float(null.mean()), null_p95=float(np.percentile(null, 95)),
                p=float((1 + (null >= delta).sum()) / (n_perm + 1)))


def merge_cells(folder, rois=None):
    """Cell-level AnnData across ROIs, in the SAME raster frame as merge_featurized.

    The cells keep their own coordinates in the source files (microns, ROI-centred),
    which do NOT line up with window coordinates. This applies the forward transform
    (subtract the per-ROI minimum, add the padding) so cells and windows overlay.
    """
    import anndata as _ad
    from anndata import read_h5ad

    blocks = []
    for f in sorted(Path(folder).glob("*.h5ad")):
        if rois is not None and f.stem not in set(rois):
            continue
        a = read_h5ad(f)
        sd = a.uns["sparse_image_state_dict"]
        px = float(sd["pixel_size"])
        pad = float(sd.get("padding", 0))
        x = a.obs["x"].values.astype(float)
        y = a.obs["y"].values.astype(float)
        cx = (x - x.min()) + pad * px
        cy = (y - y.min()) + pad * px

        X = np.asarray(a.obsm["cell_type_proportions"], dtype=np.float32)
        chans = channels_for(X.shape[1])            # 10, or 11 for an ECM run
        c = _ad.AnnData(X=X)
        c.obs_names = ["{}_c{:06d}".format(f.stem, i) for i in range(a.n_obs)]
        c.var_names = list(chans)
        c.obs["roi"] = f.stem
        c.obs["patient"] = str(a.uns.get("patient", roi_to_patient(f.stem)))
        # categories=chans, NOT CELL_TYPES: with the 10-type list every ECM row becomes
        # NaN and the matrix rows vanish from every cell panel without a warning.
        c.obs["cell_type"] = pd.Categorical(np.asarray(a.obs["cell_type"].values),
                                            categories=chans)
        c.obs["is_ecm"] = np.asarray(a.obs["is_ecm"].values, dtype=bool) \
            if "is_ecm" in a.obs.columns else np.zeros(a.n_obs, dtype=bool)
        c.obs["x"] = cx
        c.obs["y"] = cy
        c.obsm["spatial"] = np.c_[cx, cy]
        blocks.append(c)

    merged = _ad.concat(blocks, axis=0, join="outer", merge="same", index_unique=None)
    for col in ("roi", "patient"):
        merged.obs[col] = merged.obs[col].astype("category")
    merged.uns["channels"] = list(blocks[0].var_names)
    return merged


def cluster_embedding_tm(features, n_neighbors=15, resolution=1.0, pca_var=0.9, seed=0):
    """TissueMosaic's OWN clustering path: SmartPca -> SmartUmap -> Leiden on the UMAP graph.

    This mirrors tissuemosaic.utils.embedding_util.cluster without requiring the
    AnnData / sparse_image_state_dict plumbing that function expects, so it can run on
    a plain feature matrix.

    Differs from cluster_embedding (scanpy) in three ways: PCA to a variance fraction
    rather than a fixed rank, Leiden on the UMAP graph rather than a separate kNN
    graph, and RBConfiguration partitioning. MEASURED: at matched n_neighbors and
    resolution the two agree to 1.7% on gain (3.745 vs 3.808 on final_fov128 @ep399),
    so the algorithm choice is not what moves the numbers. n_neighbors and resolution
    are.

    DO NOT COPY THE TUTORIAL'S PARAMETERS. docs/source/tutorial.rst calls
    cluster(..., n_neighbors=100, leiden_res=[0.1, 0.2, 0.3]) and the function itself
    defaults to n_neighbors=500. Those suit their Slide-seq testis dataset, which has
    far more patches; on 637 IMC windows n_neighbors=100 connects each node to 16% of
    the graph and Leiden at res 0.1-0.3 collapses everything to 1-2 clusters:

        nn=15  res=1.0  -> 16 clusters, gain 3.75x
        nn=100 res=1.0  -> 10 clusters, gain 2.44x
        nn=100 res=0.3  ->  2 clusters, gain 1.92x   <- tutorial values
        nn=100 res=0.1  ->  1 cluster,  gain undefined

    Reproducible without extra effort: SmartUmap hardcodes random_state=0 and
    SmartLeiden.cluster defaults to random_state=0.
    """
    import torch
    from tissuemosaic.utils import SmartPca, SmartUmap, SmartLeiden

    # SmartPca.fit jitters the covariance matrix with an unseeded torch.randn before its
    # SVD, which makes low-dimensional inputs cluster differently run to run. See
    # umap_graph() for the measurement.
    torch.manual_seed(seed)
    X = np.asarray(features, dtype=np.float32)
    pca = SmartPca(preprocess_strategy="z_score")
    emb = pca.fit_transform(torch.tensor(X), n_components=pca_var)
    um = SmartUmap(n_neighbors=min(n_neighbors, len(X) - 1), preprocess_strategy="raw",
                   n_components=2, min_dist=0.5, metric="euclidean")
    u = um.fit_transform(emb)
    lab = SmartLeiden(graph=um.get_graph(), directed=True).cluster(
        resolution=resolution, random_state=seed, partition_type="RBC")
    return {"labels": np.asarray(lab).astype(int),
            "umap": to_numpy(u), "pca": to_numpy(emb)}


def safe_gain(kappa, kappa_comp, floor=0.02):
    """kappa / kappa_comp, but NaN when the baseline carries no spatial signal.

    The composition baseline can legitimately land at or below zero -- JHH387ROI1 came
    in at kappa_comp = -0.028, meaning composition clusters were slightly ANTI-coherent.
    Dividing by that produced a reported gain of 4e8, which is noise wearing a headline.
    Below `floor` the ratio is meaningless and kappa should be read directly instead.
    """
    kappa_comp = float(kappa_comp)
    if not np.isfinite(kappa_comp) or kappa_comp < floor:
        return float("nan")
    return float(kappa) / kappa_comp


def cluster_embedding(features, n_neighbors=15, resolution=1.0, n_pcs=30, seed=0):
    """PCA -> kNN graph -> Leiden -> UMAP, via scanpy."""
    import anndata as _ad
    import scanpy as sc

    a = _ad.AnnData(X=np.asarray(features, dtype=np.float32))
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(
        a, n_comps=min(n_pcs, min(a.shape) - 1), svd_solver="arpack", random_state=seed
    )
    sc.pp.neighbors(
        a, n_neighbors=min(n_neighbors, a.n_obs - 1), use_rep="X_pca", random_state=seed
    )
    sc.tl.leiden(
        a, resolution=resolution, key_added="leiden", flavor="igraph", random_state=seed
    )
    sc.tl.umap(a, random_state=seed)
    return {
        "labels": a.obs["leiden"].to_numpy().astype(int),
        "umap": np.asarray(a.obsm["X_umap"]),
        "pca": np.asarray(a.obsm["X_pca"]),
    }


# ------------------------------------------------------------ metrics -----


def spatial_coherence(centers_um, roi, labels, k=4, n_perm=200, seed=0):
    """Do neighbouring windows share a cluster more than chance?

    Reported against a within-ROI label permutation, because the raw value
    depends on grid size and cluster count, which differ between configs.

    Returns two standardisations, and WHICH ONE YOU WANT DEPENDS ON THE QUESTION:

      kappa  (obs - null) / (1 - null): the share of the available headroom above
             chance that the clustering actually captures. A probability
             difference, so it does NOT grow with the number of windows. This is
             the one to compare ACROSS configs.
      z      (obs - null) / null_std: significance. null_std shrinks as
             ~1/sqrt(n_pairs), so z grows as ~sqrt(N) at a FIXED effect size.
             Across a sweep whose window counts span 118 to 1078 (a 3.0x z
             inflation from N alone) it ranks small FOVs above large ones for
             purely arithmetic reasons. Use it to ask "is this above chance at
             all", never to rank configs of different FOV.

    The round-1 A-vs-B verdict (z 23.5 vs 12.2) was read off z at N 272 vs 123,
    so ~1.5x of that gap was N. A still led on kappa, but by less than z implied.
    """
    from scipy.spatial import cKDTree

    rng = np.random.default_rng(seed)
    labels = np.asarray(labels)

    pairs = []
    for r in np.unique(roi):
        m = np.flatnonzero(roi == r)
        if len(m) <= k:
            continue
        tree = cKDTree(centers_um[m])
        _, idx = tree.query(centers_um[m], k=min(k, len(m) - 1) + 1)
        for a_, row in zip(m, idx):
            for b_ in m[row[1:]]:
                pairs.append((a_, b_))
    if not pairs:
        return {
            "observed": np.nan,
            "null_mean": np.nan,
            "null_std": np.nan,
            "n_windows": len(labels),
            "n_pairs": 0,
            "kappa": np.nan,
            "z": np.nan,
        }
    pairs = np.asarray(pairs)

    obs = float((labels[pairs[:, 0]] == labels[pairs[:, 1]]).mean())
    null = np.empty(n_perm)
    for i in range(n_perm):
        shuf = labels.copy()
        for r in np.unique(roi):  # permute WITHIN ROI, preserving composition
            m = np.flatnonzero(roi == r)
            shuf[m] = rng.permutation(shuf[m])
        null[i] = (shuf[pairs[:, 0]] == shuf[pairs[:, 1]]).mean()
    sd, mu = null.std(), float(null.mean())
    return {
        "observed": obs,
        "null_mean": mu,
        "null_std": float(sd),
        "n_windows": len(labels),
        "n_pairs": len(pairs),
        "kappa": float((obs - mu) / (1.0 - mu)) if mu < 1.0 else np.nan,
        "z": float((obs - mu) / sd) if sd > 0 else np.nan,
    }


def coherence_ceiling(centers_um, roi, domain_um, n_layout=12, k=4, n_perm=100, seed=0):
    """The coherence_kappa a PERFECT model would score on THIS window grid.

    kappa strips z's ~sqrt(N) inflation but not every geometry effect, and the
    residual is large. Lay down synthetic tissue domains of a fixed physical size
    and label each window by the domain its centre falls in -- that is ground
    truth, so its kappa is the ceiling. Measured on JHH369's real window grids:

        config          N     D250   D400   D600   D800
        E  64 um     1078    0.663  0.759  0.825  0.839
        F  96 um      499    0.515  0.642  0.743  0.761
        A 128 um      272    0.373  0.536  0.652  0.681
        B 192 um      123    0.172  0.343  0.488  0.547

    At 250 um domains a perfect model scores 3.9x higher at 64 um than at 192 um
    purely because a 192 um window IS roughly one domain, so neighbours rarely
    share one, while sixteen 64 um windows tile the same domain. RAW KAPPA IS
    THEREFORE NOT COMPARABLE ACROSS FOV EITHER -- divide by this first.

    domain_um is unknown for real tissue, so report the normalised value at
    several sizes and check the RANKING is stable. JHH369 is coarser than 250 um:
    config A's observed kappa (0.415) exceeds its own D250 ceiling (0.373).

    Domains are Voronoi cells around Poisson seeds at density 1/domain_um^2,
    averaged over n_layout draws.
    """
    from scipy.spatial import cKDTree

    centers_um, roi = np.asarray(centers_um), np.asarray(roi)
    out = []
    for rep in range(n_layout):
        rng = np.random.default_rng(1000 * rep + int(domain_um) + seed)
        lab = np.empty(len(centers_um), dtype=int)
        for r in np.unique(roi):
            m = np.flatnonzero(roi == r)
            p = centers_um[m]
            lo, hi = p.min(0) - domain_um, p.max(0) + domain_um
            nseed = max(2, int(round(np.prod(hi - lo) / (domain_um**2))))
            lab[m] = cKDTree(rng.uniform(lo, hi, size=(nseed, 2))).query(p)[1]
        out.append(
            spatial_coherence(centers_um, roi, lab, k=k, n_perm=n_perm, seed=seed)[
                "kappa"
            ]
        )
    return float(np.mean(out))


def roi_mixing(roi, labels, n_perm=200, seed=0):
    """Do clusters span ROIs, or just recover ROI identity?

    Size-weighted mean of each cluster's ROI entropy, normalised to [0, 1].
    1.0 = clusters are ROI-agnostic; 0.0 = each cluster is one ROI.

    DO NOT READ THIS AS "HIGHER IS BETTER". Read it against the value you get by
    clustering the raw composition vector with no model at all, because cell-type
    proportions recur across ROIs trivially and therefore score high by default:

        config           model embedding    composition baseline
        A 128 um              0.650                0.837
        E  64 um              0.898                0.906

    E looks like the winner and is not. Its 0.898 sits ON its own baseline, so the
    high value means "the embedding added nothing beyond composition", not "found
    shared motifs" -- and indeed E's model buys only +23% spatial coherence over
    that baseline against A's +203%. A's 0.650 is far BELOW its baseline, meaning
    A's clusters carry something composition does not.

    What A's low value cannot distinguish is real per-ROI tissue architecture from
    batch. Four ROIs of one patient are not enough; that needs the full cohort.
    """
    roi, labels = np.asarray(roi), np.asarray(labels)
    uroi = np.unique(roi)
    if len(uroi) < 2:
        return {"normalised_entropy": np.nan, "null_mean": np.nan}
    maxent = np.log(len(uroi))

    def ent(lab):
        tot, acc = 0, 0.0
        for c in np.unique(lab):
            m = lab == c
            p = np.array([(roi[m] == r).mean() for r in uroi])
            p = p[p > 0]
            acc += m.sum() * float(-(p * np.log(p)).sum() / maxent)
            tot += m.sum()
        return acc / tot

    rng = np.random.default_rng(seed)
    null = np.array([ent(rng.permutation(labels)) for _ in range(n_perm)])
    return {"normalised_entropy": ent(labels), "null_mean": float(null.mean())}


def composition_silhouette(composition, labels):
    """Do the clusters correspond to distinct cell-type compositions?"""
    from sklearn.metrics import silhouette_score

    labels = np.asarray(labels)
    if len(np.unique(labels)) < 2:
        return np.nan
    return float(silhouette_score(np.asarray(composition), labels, metric="euclidean"))


def border_mask(centers_um, roi):
    """True for windows in the outermost ring of their ROI's tiling grid.

    Border windows sit partly outside the tissue -- in v1 they held a median 42
    cells against 106 interior, and formed 2 of 9 clusters (one 100% border, one
    86%). A border ring is contiguous by construction, so it inflates spatial
    coherence.

    Its share depends on the analysis filter, which is why that filter must scale
    with FOV. Under the old fixed 50-cell threshold the share ran 32% at 128 um
    against 51% at 192 -- the filter culled sparse border tiles at small FOV and
    spared them at large, biasing exactly the comparison it was meant to protect.
    Under fov_scaled_threshold the six configs sit at 13.6 / 26.8 / 15.3 / 28.2 /
    16.9 / 19.4 % (A-F). Better, and no longer FOV-monotone, but a ~2x spread
    remains at 192 um, so read cross-FOV coherence with that in mind.
    """
    centers_um = np.asarray(centers_um)
    roi = np.asarray(roi)
    edge = np.zeros(len(roi), dtype=bool)
    for r in np.unique(roi):
        m = np.flatnonzero(roi == r)
        c = np.round(centers_um[m], 3)
        xs, ys = np.unique(c[:, 0]), np.unique(c[:, 1])
        edge[m] = np.isin(c[:, 0], [xs[0], xs[-1]]) | np.isin(c[:, 1], [ys[0], ys[-1]])
    return edge


def border_cluster_report(centers_um, roi, labels):
    """Per-cluster share of border windows -- the v1 border-cluster diagnostic."""
    edge = border_mask(centers_um, roi)
    labels = np.asarray(labels)
    rows = []
    for c in np.unique(labels):
        m = labels == c
        rows.append(
            {
                "cluster": int(c),
                "n_windows": int(m.sum()),
                "n_border": int(edge[m].sum()),
                "pct_border": round(100.0 * edge[m].mean(), 1),
            }
        )
    df = pd.DataFrame(rows).set_index("cluster")
    return df, {
        "overall_border_pct": round(100.0 * edge.mean(), 1),
        "n_border_clusters_ge80pct": int((df.pct_border >= 80).sum()),
        "max_cluster_border_pct": float(df.pct_border.max()),
    }


def window_metrics(result, labels, k=4, n_perm=200, seed=0):
    """All comparison metrics for one sweep config.

    Rank configs on `coherence_kappa`, not `spatial_coherence_z` -- the sweep's
    window counts span 118 to 1078 and z grows as ~sqrt(N). See spatial_coherence.
    """
    sc_ = spatial_coherence(
        result["centers_um"], result["roi"], labels, k=k, n_perm=n_perm, seed=seed
    )
    rm = roi_mixing(result["roi"], labels, n_perm=n_perm, seed=seed)
    return {
        "n_windows": len(labels),
        "n_clusters": len(np.unique(labels)),
        "spatial_coherence": sc_["observed"],
        "coherence_kappa": sc_["kappa"],  # N-independent -- rank on this
        "spatial_coherence_z": sc_["z"],  # scales ~sqrt(N); significance only
        "roi_mixing": rm["normalised_entropy"],
        "roi_mixing_null": rm["null_mean"],
        "composition_silhouette": composition_silhouette(result["composition"], labels),
    }


# ------------------------------------------------------------- plotting -----


def _srgb_to_lab(rgb):
    """sRGB in [0,1] -> CIELAB (D65). Written out because colorspacious and
    scikit-image are both absent from this environment."""
    c = np.asarray(rgb, dtype=float)
    lin = np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)
    m = np.array(
        [
            [0.4124564, 0.3575761, 0.1804375],
            [0.2126729, 0.7151522, 0.0721750],
            [0.0193339, 0.1191920, 0.9503041],
        ]
    )
    xyz = lin @ m.T / np.array([0.95047, 1.0, 1.08883])
    e, k = 216 / 24389, 24389 / 27
    f = np.where(xyz > e, np.cbrt(xyz), (k * xyz + 16) / 116)
    return np.stack(
        [
            116 * f[..., 1] - 16,
            500 * (f[..., 0] - f[..., 1]),
            200 * (f[..., 1] - f[..., 2]),
        ],
        axis=-1,
    )


def _glasbey_fallback(
    n, l_range=(18.0, 72.0), chroma_min=18.0, grid=12, seed_rgb=(0.10, 0.35, 0.85)
):
    """Glasbey construction: greedily take the candidate colour whose minimum
    distance to everything already chosen is largest.

    This is the method of Glasbey et al. 2007, which used CIELAB exactly as here.
    colorcet's glasbey_dark instead optimises in CAM02-UCS, so the two are not
    identical -- `cluster_palette` uses the real one whenever colorcet imports.

    l_range restricts to darker colours, which is what makes glasbey_dark legible
    against the white figure background; chroma_min drops washed-out greys that
    are hard to tell apart at scatter-point size.
    """
    g = np.linspace(0, 1, grid)
    cand = np.stack(np.meshgrid(g, g, g, indexing="ij"), axis=-1).reshape(-1, 3)
    lab = _srgb_to_lab(cand)
    chroma = np.hypot(lab[:, 1], lab[:, 2])
    keep = (
        (lab[:, 0] >= l_range[0]) & (lab[:, 0] <= l_range[1]) & (chroma >= chroma_min)
    )
    cand, lab = cand[keep], lab[keep]

    chosen = [np.asarray(seed_rgb, dtype=float)]
    dmin = np.linalg.norm(lab - _srgb_to_lab(chosen[0]), axis=-1)
    while len(chosen) < n:
        i = int(np.argmax(dmin))
        chosen.append(cand[i])
        dmin = np.minimum(dmin, np.linalg.norm(lab - lab[i], axis=-1))
    return [tuple(c) for c in chosen[:n]]


# Every categorical variable in these notebooks is coloured from glasbey_dark, but
# each family takes a DIFFERENT slice so that cluster 3, a cell type and an ROI never
# share a colour -- several figures put those legends side by side. The offsets are
# not arbitrary: measured minimum separation in CIELAB shows each slice beats the
# matplotlib palette it replaced.
#
#   family      slice      min dLab    was
#   clusters    [ 0:n ]      23.7      tab20  16.6   (n=15)
#   cell types  [16:26]      37.5      tab10  27.7
#   ROIs        [40:44]      49.3      tab10  40.3
#
# (The cluster row is measured in CIELAB, which understates colorcet's palette --
# it optimises in CAM02-UCS. Both beat tab20 on the same yardstick, which is the
# comparison that matters.)
GLASBEY_OFFSET = {"cluster": 0, "cell_type": 16, "roi": 40}


def _glasbey(n, offset=0, stride=1):
    """n colours from colorcet.glasbey_dark, from `offset`, taking every `stride`-th.

    stride=2 at offset 0 vs offset 1 gives two DISJOINT palettes -- how two cluster
    sets shown side by side avoid sharing colours.
    """
    try:
        import colorcet
        from matplotlib.colors import to_rgb

        base = colorcet.glasbey_dark
        idx = list(range(offset, len(base), stride))[:n]
        cols = [to_rgb(base[i]) for i in idx]
        while len(cols) < n:  # ran off the end of a 256-colour palette
            cols.append(to_rgb(base[len(cols) % len(base)]))
        return [tuple(c) for c in cols]
    except ImportError:
        return _glasbey_fallback(offset + n * stride)[offset::stride][:n]


def roi_palette(rois):
    """Stable ROI -> colour mapping, glasbey_dark offset away from the other families."""
    rois = (
        list(dict.fromkeys(np.asarray(rois).tolist()))
        if not isinstance(rois, int)
        else None
    )
    cols = _glasbey(len(rois), offset=GLASBEY_OFFSET["roi"])
    return dict(zip(rois, cols))


def cluster_palette(n=None, labels=None, offset=None, stride=1):
    """Maximally distinct colours for cluster labels -- glasbey_dark.

    Leiden gives 6-15 clusters here, and matplotlib's tab20 wraps into
    near-duplicate pairs well before that, which made neighbouring clusters
    indistinguishable in the tissue maps. Glasbey palettes are built for exactly
    this: each successive colour is the one furthest from every colour already in
    the palette, so distinctness degrades gracefully as the count grows.

    Uses colorcet.glasbey_dark when colorcet is installed. It is NOT installed in
    tissue-mosaic-260819, so a local CIELAB Glasbey construction stands in --
    same algorithm, different colour space (see _glasbey_fallback). Installing
    colorcet switches to the real palette with no code change.

    Call with `n` for a list of colours, or with `labels` for a {label: colour}
    mapping covering exactly those labels. One of the two is required.
    """
    if n is None and labels is None:
        raise TypeError("cluster_palette() needs either n or labels")
    if labels is not None:
        uniq = list(dict.fromkeys(np.asarray(labels).tolist()))
        cols = cluster_palette(len(uniq), offset=offset, stride=stride)
        return dict(zip(uniq, cols))

    return _glasbey(n, offset=GLASBEY_OFFSET["cluster"] if offset is None else offset,
                    stride=stride)


def cluster_palette_source():
    """Which palette cluster_palette() is actually returning, for reporting."""
    try:
        import colorcet  # noqa: F401

        return "colorcet.glasbey_dark"
    except ImportError:
        return "local CIELAB Glasbey fallback (colorcet not installed)"


def plot_window_map(ax, centers_um, colors, fov_um, edge=None, lw=0.0):
    """Draw each window as a true square of its physical size, in data coordinates.

    The obvious approach -- ax.scatter with marker="s" -- sizes markers in POINTS,
    so the squares only tile the tissue if you happen to guess a size that matches
    the figure geometry. Guess low and the map becomes a sparse dot grid with white
    gutters between windows; guess high and neighbours overlap. Either way the one
    thing this panel exists to show, whether same-cluster windows form CONTIGUOUS
    blocks, is exactly what gets obscured.

    Rectangles in data coordinates are immune: a 128 um window is drawn 128 um wide
    at any figure size, zoom or dpi, so adjacent windows always share an edge and
    contiguity reads honestly.

    centers_um is (N, 2), colors is a sequence of N matplotlib colours, fov_um is
    the window side in microns (pixel_size x global_size).
    """
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle

    centers_um = np.asarray(centers_um, dtype=float)
    half = float(fov_um) / 2.0
    rects = [Rectangle((x - half, y - half), fov_um, fov_um) for x, y in centers_um]
    pc = PatchCollection(
        rects,
        facecolors=list(colors),
        edgecolors=(edge or "none"),
        linewidths=lw,
        match_original=False,
    )
    ax.add_collection(pc)
    ax.set_xlim(centers_um[:, 0].min() - half, centers_um[:, 0].max() + half)
    ax.set_ylim(centers_um[:, 1].min() - half, centers_um[:, 1].max() + half)
    ax.set_aspect("equal")
    return pc


# other_tissue and duct_filler are not measured biology -- they are unclassified fill
# and hex-packed synthetic duct. Together ~68% of every window, so giving them saturated
# hues made them dominate every tissue panel and swamp the stacked bars. Greys push them
# back without hiding them: light for the diffuse other_tissue, dark for the structured
# duct_filler, which reads as tissue architecture rather than noise.
BACKGROUND_GREYS = {"other_tissue": "#D8D8D8", "duct_filler": "#4A4A4A"}


def palette(types=None, background_grey=True):
    """Stable cell-type -> colour mapping, shared across all notebooks.

    The eight BIOLOGY_TYPES take glasbey_dark from GLASBEY_OFFSET["cell_type"], chosen
    so cell-type colours never collide with cluster colours -- several figures show both
    legends at once.

    background_grey=True (default) overrides the two BACKGROUND_TYPES with greys; see
    BACKGROUND_GREYS. Pass False for the raw 10-colour glasbey mapping.
    """
    from matplotlib.colors import to_rgb

    types = list(types) if types is not None else CELL_TYPES
    if not background_grey:
        return dict(zip(types, _glasbey(len(types), offset=GLASBEY_OFFSET["cell_type"])))

    bio = [t for t in types if t not in BACKGROUND_GREYS]
    cols = _glasbey(len(bio), offset=GLASBEY_OFFSET["cell_type"])
    out = dict(zip(bio, cols))
    for t in types:
        if t in BACKGROUND_GREYS:
            out[t] = to_rgb(BACKGROUND_GREYS[t])
    return {t: out[t] for t in types}          # keep the caller's ordering


def plot_tissue(
    adata,
    ax=None,
    biology_only=False,
    domain=None,
    s=2.0,
    background_color="0.88",
    colors=None,
    legend=False,
):
    """Scatter one ROI in its native frame.

    y is plotted as-is: the ICS frame already has y increasing upward, so the
    matplotlib default reproduces the original anatomical orientation. Never
    invert or negate it.
    """
    import matplotlib.pyplot as plt

    colors = colors or palette()
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))

    if biology_only:
        bg = adata.obs["cell_type"].isin(BACKGROUND_TYPES).values
        ax.scatter(
            adata.obs.x[bg],
            adata.obs.y[bg],
            s=s * 0.5,
            c=background_color,
            linewidths=0,
        )
        draw = BIOLOGY_TYPES
        s = s * 2
    else:
        draw = CELL_TYPES

    for t in draw:
        m = (adata.obs["cell_type"] == t).values
        if m.sum():
            ax.scatter(
                adata.obs.x[m],
                adata.obs.y[m],
                s=s,
                c=[colors[t]],
                label=t,
                linewidths=0,
            )

    if domain is not None:
        import matplotlib.pyplot as _plt

        x0, x1, y0, y1 = domain
        ax.add_patch(
            _plt.Rectangle(
                (x0, y0), x1 - x0, y1 - y0, fill=False, ec="0.4", lw=1, ls="--"
            )
        )

    ax.set_aspect("equal")
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    if legend:
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", markerscale=5, fontsize=8)
    return ax


def plot_channels(sparse_image, pixel_size=None, ncols=5, figsize=None):
    """Montage of the rasterised SparseImage, one panel per channel.

    origin='lower' so the raster keeps the same orientation as plot_tissue.
    """
    import matplotlib.pyplot as plt

    dense = to_numpy(sparse_image.data)
    n = len(CELL_TYPES)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=figsize or (4 * ncols, 4.8 * nrows),
        constrained_layout=True,
    )
    for i, (ax, t) in enumerate(zip(np.asarray(axes).ravel(), CELL_TYPES)):
        ax.imshow(dense[i], origin="lower", cmap="viridis")
        ax.set_title(f"ch {i}: {t}\n({int(dense[i].sum()):,} cells)", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in np.asarray(axes).ravel()[n:]:
        ax.axis("off")
    return fig, axes


# Fail loudly at import if the channel definition has drifted from the yaml.
CHANNEL_CHECK = check_channels()

__all__ = [
    "BACKGROUND_TYPES",
    "BIOLOGY_TYPES",
    "CELL_TYPES",
    "CONFIG_YAML",
    "H5AD_DIR",
    "ICS_DIR",
    "QUPATH_DIR",
    "REPO",
    "ROI_SPECS",
    "SKIP_PATTERN",
    "SWEEP_DENSITY_PER_MM2",
    "SWEEP_H5AD_DIR",
    "SWEEP_OVERRIDES",
    "SWEEP_PATIENT",
    "SWEEP_RUN_DIR",
    "WINDOW_SWEEP",
    "border_cluster_report",
    "border_mask",
    "build_anndata",
    "categories_to_channels",
    "check_channels",
    "cluster_embedding",
    "compat_summary",
    "composition",
    "composition_silhouette",
    "csv_to_anndata",
    "featurize_windows",
    "fov_scaled_threshold",
    "list_ics_csvs",
    "load_all_anndata",
    "load_config",
    "load_model",
    "load_roi_domains",
    "make_datamodule",
    "palette",
    "plot_channels",
    "plot_tissue",
    "read_loss_curve",
    "repo_root",
    "roi_metadata",
    "roi_mixing",
    "roi_to_patient",
    "spatial_coherence",
    "sweep_checkpoint",
    "sweep_config_dict",
    "sweep_run_dir",
    "sweep_spec",
    "sweep_status",
    "sweep_subset_folder",
    "tissuemosaic_run_script",
    "to_numpy",
    "window_metrics",
]
