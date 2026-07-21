#!/bin/bash
# wrap_driver.sh -- generic job script; this is what actually runs on the
# compute node. Takes the simulation id as $1 (e.g. htan_wellmixed) and the
# slurm/ directory's absolute path as $2, and runs $2/run_<id>.jl.
# submit_driver.sh passes both when it submits.
#
# $2 must be passed in rather than self-located via ${BASH_SOURCE[0]}: sbatch
# copies the submitted script into its own spool directory
# (/var/spool/slurm/...) and executes that copy, so self-location would
# resolve to the spool dir, not the real slurm/ directory.
# No -u (nounset): the conda-provided julia_activate.sh hook references
# $JULIA_DEPOT_PATH without a safe default on one line, which trips nounset
# even though the variable is legitimately unset in a fresh job environment.
set -eo pipefail

module load slurm

source /usr/local/packages/miniconda3/etc/profile.d/conda.sh
conda activate physicell-sim-260606

SIM_ID="$1"
SCRIPTS_DIR="$2"
julia "$SCRIPTS_DIR/run_${SIM_ID}.jl"
