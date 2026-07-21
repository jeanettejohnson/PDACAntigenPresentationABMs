#!/bin/bash
# wrap_driver.sh -- generic job script; this is what actually runs on the
# compute node. Takes the simulation id as $1 (e.g. htan_wellmixed) and runs
# slurm/run_<id>.jl. submit_driver.sh passes this argument when it submits.
set -euo pipefail

module load slurm

source /usr/local/packages/miniconda3/etc/profile.d/conda.sh
conda activate physicell-sim-260606

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIM_ID="$1"
julia "$SCRIPT_DIR/run_${SIM_ID}.jl"
