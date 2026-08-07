#!/bin/bash
# wrap_spatial_build.sh -- generic job script; this is what actually runs on
# the compute node for the "build" stage. Takes the pipeline directory
# (contains summarize_simulations.py and .venv/) as $1.
#
# $1 must be passed in rather than self-located via ${BASH_SOURCE[0]}: sbatch
# copies the submitted script into its own spool directory, so self-location
# would resolve to the spool dir, not the real pipeline directory (same
# reasoning as wrap_driver.sh).
#
# --force is intentional here, not a default you should keep long-term: an
# earlier disk-full crash may have left truncated/corrupt per-simulation CSVs
# in results/, and build_combined() silently reuses any existing
# {sim}_spatial_features.csv without checking it's complete. --force makes
# this first post-fix run recompute everything cleanly. Drop --force on
# later incremental runs (or switch to --only-missing) once you trust the
# cache again.
set -eo pipefail

module load slurm

PIPELINE_DIR="$1"
DATA_DIR="/local/projects-t3/fertig_pdacagmodel/htan_wellmixed/PDACAntigenPresentationABMs/data/outputs/simulations"

cd "$PIPELINE_DIR"

echo "Host: $(hostname)"
echo "--- memory available on this node ---"
free -h
echo "--- disk space at pipeline dir ---"
df -h "$PIPELINE_DIR"
echo "---"

"$PIPELINE_DIR/.venv/bin/python3" -u summarize_simulations.py \
    --data-dir "$DATA_DIR" \
    --output-dir results \
    --force \
    --build-only
