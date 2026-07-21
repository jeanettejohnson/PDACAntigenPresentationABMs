# PDACAntigenPresentationABMs

## Getting started

### 1. Clone the repository (with submodules)

This repo uses git submodules for the PhysiCell engine and PhysiCell Studio. Clone with `--recurse-submodules`, otherwise `PhysiCell/` and `PhysiCell-Studio/` come down as empty directories:

```
git clone --recurse-submodules git@github.com:jeanettejohnson/PDACAntigenPresentationABMs.git
cd PDACAntigenPresentationABMs
```

### 2. Create the conda environment

```
conda env create -f ./conda_env_configs/physicell_sim_260606.yaml
```

### 3. Install PCMM in the conda environment

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

This should return `true`. Do **not** run `createProject()` here: this repo is already an existing PCMM project (`data/`, submodules, and configs are already checked in), and PCMM's `createDefaultGitIgnore` step isn't idempotent -- running it again will duplicate `data/.gitignore`. `initializeModelManager()` is the one part of `createProject()` that actually matters at runtime -- it's already called automatically by every script in step 4; this just verifies the install works now.

See the [PhysiCellModelManager.jl guide](https://github.com/drbergman-lab/PhysiCellModelManager.jl) for background on these steps.

### 4. Run scripts

todo (will be added later)