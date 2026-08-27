# PDACAntigenPresentationABMs

Agent-based models of antigen presentation in pancreatic cancer, built on
[PhysiCell](https://github.com/jeanettejohnson/PhysiCell) and driven by
[PhysiCellModelManager.jl](https://github.com/drbergman-lab/PhysiCellModelManager.jl) (PCMM).

## Repository layout

| directory | what lives there |
|-----------|------------------|
| `slurm/` | job submission -- one runner per simulation, plus shared SLURM settings |
| `data/` | PCMM project: `inputs/` (configs, rulesets, custom code, ICs) and `outputs/` (results) |
| `PhysiCell/`, `PhysiCell-Studio/` | git submodules: the engine and the viewer |
| `prep_imc_spatial/` | builds the IMC spatial initial conditions (see below) |
| `analysis/` | consume results, produce processed data |
| `plotting/` | visualization |
| `sanity_checks/` | verify the validity of data and code |
| `tissuemosaic/` | self-supervised spatial modelling of the IMC cohort with [TissueMosaic](https://github.com/broadinstitute/TissueMosaic) (see its own README) |
| `figures/`, `movies/` | generated output |

The three phases below map onto that layout: **get running** (steps 1-5),
**prepare IMC spatial inputs** if you need them, then **work with results**.

---

## Getting started

### 1. Clone the repository

```
git clone git@github.com:jeanettejohnson/PDACAntigenPresentationABMs.git
cd PDACAntigenPresentationABMs
```

### 2. Populate the submodules

This repo uses git submodules for the PhysiCell engine and PhysiCell Studio. Without this step, `PhysiCell/` and `PhysiCell-Studio/` are empty directories:

```
git submodule update --init --recursive
```

### 3. Create the conda environment

```
conda env create -f ./conda_env_configs/physicell_sim_260606.yaml
```

### 4. Install PCMM in the conda environment

```
conda activate physicell-sim-260606
julia
```

In the Julia REPL, enter Pkg mode with `]` and run:

```
pkg> registry add General
pkg> registry add https://github.com/drbergman-lab/BergmanLabRegistry
pkg> add PhysiCellModelManager CSV DataFrames
```

Back at the `julia>` prompt (press backspace to leave Pkg mode), verify the install:

```julia
using PhysiCellModelManager
initializeModelManager()
```

This should return `true`. Do **not** run `createProject()` here: this repo is already an existing PCMM project (`data/`, submodules, and configs are already checked in), and PCMM's `createDefaultGitIgnore` step isn't idempotent -- running it again will duplicate `data/.gitignore`. `initializeModelManager()` is the one part of `createProject()` that actually matters at runtime -- it's already called automatically by every script in step 5; this just verifies the install works now.

See the [PhysiCellModelManager.jl guide](https://github.com/drbergman-lab/PhysiCellModelManager.jl) for background on these steps.

### 5. Run simulations

Simulations are submitted to SLURM via scripts in [`slurm/`](slurm/). The entry point is `slurm/submit_driver.sh`, which submits a driver job for one or more of:

| # | Simulation       | Script                          |
|---|------------------|----------------------------------|
| 1 | htan_wellmixed   | `slurm/run_htan_wellmixed.jl`   |
| 2 | htan_geometries  | `slurm/run_htan_geometries.jl`  |
| 3 | imc_wellmixed    | `slurm/run_imc_wellmixed.jl`    |
| 4 | imc_spatial      | `slurm/run_imc_spatial.jl`      |

```
cd slurm
./submit_driver.sh 1 3      # submit htan_wellmixed and imc_wellmixed
./submit_driver.sh          # no arguments: interactive prompt instead
```

Each driver job activates the conda environment and runs the matching Julia script, which itself submits the individual simulation jobs -- a job that submits more jobs, not one flat batch. Job resources (CPUs, memory, time limit) and the SLURM account resolve automatically; see `slurm/slurm_common.jl` to adjust them.

Driver-job logs land in `slurm/logs/` as `<job-name>_<job-id>.out`/`.err`. Individual simulations log to `data/outputs/simulations/<id>/`.

All four simulations run straight from a fresh clone -- every PCMM input they
need (configs, rulesets, custom code, and the 48 per-ROI IC folders for
imc_spatial) is tracked in `data/inputs/`. No generation step is required.

To run a subset of imc_spatial ROIs, set `IMC_SPATIAL_ROIS`:

```
IMC_SPATIAL_ROIS=JHH368 ./submit_driver.sh 4          # the 4 JHH368 ROIs
IMC_SPATIAL_ROIS=JHH368,JHH369 ./submit_driver.sh 4   # two donors
```

That is enough to get results. Everything below is optional.

---

## Preparing IMC spatial initial conditions

> **Work in progress.** The pipeline below produced the initial conditions
> currently in the repo, but it is not yet reproducible end to end -- see
> *Known gaps*. You only need this if the ICs must be rebuilt; otherwise the
> committed ICs are ready to use and step 5 above is sufficient.

Everything that turns raw IMC imaging into the PCMM inputs `run_imc_spatial.jl`
consumes lives in [`prep_imc_spatial/`](prep_imc_spatial/). Run from the repo
root. The funnel narrows at each stage: **56 detections → 51 substrates → 50
cell ICs → 48 ROI configs → 48 PCMM input pairs.**

| # | script | reads | writes |
|---|--------|-------|--------|
| 1 | `scale_biwt.py` | `ics/substrates/*_ecm.csv` | `*_ecm_scaled.csv` (95th-percentile clip) |
| 2 | `batch_assemble_ics.py` → `assemble_initial_conditions.py` | `qupath_detections/*.txt`, `ics/ductcoordinates/*.geojson`, `*_ecm_scaled.csv` | `ics/JHH_IMC/<ROI>.csv` |
| 3 | `generate_roi_configs.py` | cell ICs + substrates | `PhysiCell_settings_<ROI>.xml` ×48 |
| 4 | `make_imc_count_summary.py` | cell ICs | `prep_imc_spatial/assignmentsummary_JHH_IMC.csv` |
| 5 | `setup_imc_spatial_pcmm.py` | the 48 ROI configs | `data/inputs/` + `prep_imc_spatial/imc_spatial_roi_specs.csv` |

Stage 2 is where antigen presentation is encoded: `TYPE_MAP` maps QuPath classes
onto PhysiCell cell types (`tumor_epithelial: HLA-DR` → `epithelial_tumor_class1_class2`,
`CD4 T cell: FOXP3` → `Treg`, ...), and duct polygons are hex-packed with `duct_filler` cells.

Stage 5 is the PCMM bridge. The custom code under
`data/inputs/custom_codes/antigen_presentation/` is **not** generated -- it is
tracked in git and edited directly, matching `antigen_presentation_htan_singlecell`.
It used to be re-copied from the user project on every run, which silently
destroyed hand edits (the SVG palette among them).

**Ordering constraint:** if the ICs are rebuilt, stage 5 must be re-run. The PCMM
`ic_cell`/`ic_substrate` folders are copies one step downstream, and
`run_imc_spatial.jl` cannot tell that they have gone stale.

`prep_imc_spatial/archive/` holds superseded and one-off scripts from this
pipeline, kept for provenance. See [`prep_imc_spatial/README.md`](prep_imc_spatial/README.md)
for per-script detail.

### Known gaps

- **Stage 1 is not reproducible here.** The ECM channel TIFF → `_ecm.csv` step
  (`archive/scaleimage.py`) has a hardcoded personal path and a GUI file picker.
  Only its committed outputs exist.
- **Three scripts read a stale directory.** `generate_roi_configs.py` and two
  plotters read the repo-root `PhysiCell/config/ics/JHH_IMC`, which is an older
  generation with a different schema (22 files, volumes all `NaN`) rather than the
  canonical `user_projects/.../ics/JHH_IMC` (50 files). Repoint before relying on them.

---

## Working with results

> **Work in progress.** Most scripts in `analysis/` and `plotting/` predate the
> move to PCMM and still read the old `PhysiCell/outputs/<ROI>/` layout, which no
> longer exists -- PCMM writes to `data/outputs/simulations/<id>/output/`.
> **16 of the 21 scripts in these two directories need repointing before they run.**
> They are filed by purpose; that filing does not mean they currently work.

### `analysis/` -- consume results, produce processed data

| script | purpose |
|--------|---------|
| `antigen_class_comparison.py`<br>`antigen_class_comparison_imc.py`<br>`antigen_class_comparison_onedrive.py` | antigen-presentation class proportions, initial vs final. Three separate implementations (HTAN, IMC, OneDrive source) -- they differ substantially, so they are not copies to be deduplicated |
| `cell_type_counts_timeseries.py` | cell-type counts at every timepoint |
| `tcell_comparison.py` | T-cell subset composition, initial vs final |
| `load_HT_simulations.py` | loader for HT* simulation output folders |
| `analyzehtandata.py`, `loadhtanwustlscrnaseq.py` | HTAN scRNAseq: Synapse download and analysis |
| `images.py`, `image_annotations.py` | imaging preprocessing (Akoya MIF, QuPath annotations) |

Outputs go to `figures/`.

### `plotting/` -- visualization

| script | purpose |
|--------|---------|
| `plot_initial_counts.py`, `plot_initial_counts_imc.py` | initial cell counts per simulation / per ROI |
| `plot_agent_counts.py`, `plot_agent_counts_per_sample.py` | agent-assignment counts, aggregated and per sample |
| `plot_antigen_comparison.py`, `plot_antigen_presentation_timeseries.py` | antigen-presentation figures |
| `plot_apcaf_vs_tcells.py` | CAF-family vs T-cell counts per ROI |
| `plot_cell_size_by_type.py` | max-diameter distributions per cell type |
| `plot_ble_cell_counts.py`, `plot_therapy_cell_counts.py` | **byte-identical to each other** despite the names; both read a nonexistent external repo |
| `plot_invasiveness.py` | reads `../Downloads/` -- external |

All write PNGs to `figures/`, which is not created automatically by every script.

### Movies

`make_movies.jl` converts a simulation's SVG snapshots to an mp4 via ImageMagick
and ffmpeg (both in the conda environment; ffmpeg is *not* on the bare login shell,
so activate first).

```
julia make_movies.jl          # every simulation in data/outputs/simulations/
julia make_movies.jl 2 5      # specific simulation ids
```

Movies are written to `movies/<id>_<ROI>.mp4`. The ROI label comes from the
database, so the numeric simulation id is not the only handle. JPEG intermediates
are written into the simulation's output folder and removed afterwards -- an
interrupted run can leave stray `.jpg` files there.

### PhysiCell Studio

Studio visualizes a single simulation interactively. Launch it through PCMM's
`runStudio` -- see the [PCMM Studio guide](https://drbergman-lab.github.io/PhysiCellModelManager.jl/stable/man/physicell_studio/).

```
conda activate physicell-sim-260606
julia
```

```julia
using PhysiCellModelManager
runStudio(3, python_path="python", studio_path="PhysiCell-Studio")
```

Replace `3` with the simulation id you want to view; it must already have run.
`python_path`/`studio_path` are remembered for the session, so later calls can be
just `runStudio(<id>)`.

**Do not use the "Run" tab inside Studio** -- it can delete simulation data.
Studio opens with temporary config files; if you make edits you want to keep, use
"File > Save As", since they are otherwise lost when Studio closes.

---

## Maintenance

### Sanity checks

[`sanity_checks/`](sanity_checks/) holds scripts that verify the validity of data
and code, as opposed to producing results.

```
python3 sanity_checks/check_cell_colors.py
```

Run it after changing any cell-type palette, or after bumping the PhysiCell or
PhysiCell-Studio submodule.

### Cell-type colors

The SVG palette is hardcoded in `my_coloring_function` in each model's
`custom.cpp`, keyed by cell type name, and mirrored in
`PhysiCell-Studio/bin/cmaps.py` (`paint_clist`) for Studio's scalar-coloring view.
Studio's default cell view reads the SVG fills directly, so the C++ side covers it.

It deliberately does **not** live in `PhysiCell/modules/PhysiCell_pathology.cpp`:
it did once, and a submodule bump that rewrote `paint_by_number_cell_coloring()`
silently discarded it -- movies kept rendering, just in the wrong colors.
`user_projects/` is ours alone, so upstream merges cannot reach it.

Because the palette is duplicated across four `custom.cpp` copies plus Studio,
`check_cell_colors.py` is what keeps them honest; it fails loudly if any copy
drifts or falls out of step with `<cell_definitions>`.
