# PDACAntigenPresentationABMs

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

Simulations are submitted to SLURM via scripts in the [`slurm/`](slurm/) directory. The entry point is `slurm/submit_driver.sh`, which submits a driver job for one or more of the following simulations:

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

Each driver job activates the conda environment and runs the matching Julia script, which itself submits the individual simulation (or per-ROI, for imc_spatial) jobs to SLURM -- so this is a job that submits more jobs, not one flat batch. Job resources (CPUs, memory, time limit) and the SLURM account are resolved automatically; see `slurm/slurm_common.jl` to adjust them.

Driver-job logs land in `slurm/logs/` as separate `<job-name>_<job-id>.out`/`.err` files (SLURM's own `%x_%j` naming). The individual simulation jobs log elsewhere: PCMM's own `data/outputs/simulations/<id>/` for simulations 1-3, and `PhysiCell/outputs/<ROI>/run.out`/`run.err` for imc_spatial.

### 6. Visualize results with PhysiCell Studio

PhysiCell Studio (the `PhysiCell-Studio/` submodule from step 2) can visualize a simulation's output. It's launched through PCMM's `runStudio` -- see the [PCMM Studio guide](https://drbergman-lab.github.io/PhysiCellModelManager.jl/stable/man/physicell_studio/) for the full reference.

From the top level of the repo:

```
conda activate physicell-sim-260606
julia
```

```julia
using PhysiCellModelManager
runStudio(3, python_path="python", studio_path="PhysiCell-Studio")
```

Replace `3` with the ID of the simulation whose output you want to visualize -- that simulation must have already run and produced output (see step 5). `python_path`/`studio_path` are remembered for the rest of the Julia session, so later calls can just be `runStudio(<id>)`.

**Do not use the "Run" tab inside Studio** -- it can delete simulation data. Studio opens with temporary config files; if you make edits you want to keep, use "File > Save As", since they're otherwise lost when Studio closes.
