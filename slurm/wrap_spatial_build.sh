#!/bin/bash
# wrap_spatial_build.sh -- generic job script; this is what actually runs on
# the compute node for the "build" stage. Takes the pipeline directory
# (contains summarize_simulations.py) as $1.
#
# $1 must be passed in rather than self-located via ${BASH_SOURCE[0]}: sbatch
# copies the submitted script into its own spool directory, so self-location
# would resolve to the spool dir, not the real pipeline directory (same
# reasoning as wrap_driver.sh).
#
# No .venv here -- there isn't one. The working interpreter is whatever
# `module load python/3.11.4` puts on PATH (confirmed interactively: this is
# what fixed the "from __future__ import annotations" SyntaxError from the
# stock system python3). Dependencies (pandas, spatialtissuepy, sklearn,
# umap-learn, etc.) must be importable under that module -- if this job fails
# with ModuleNotFoundError instead of the old SyntaxError, they're installed
# to user-site (~/.local/lib/python3.11/site-packages) or some other
# mechanism that needs to be reproduced here too.
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
module load python/3.11.4

PIPELINE_DIR="$1"
DATA_DIR="/local/projects-t3/fertig_pdacagmodel/htan_wellmixed/PDACAntigenPresentationABMs/data/outputs/simulations"

cd "$PIPELINE_DIR"

echo "Host: $(hostname)"
echo "python3: $(which python3)"
echo "--- memory available on this node ---"
free -h
echo "--- disk space at pipeline dir ---"
df -h "$PIPELINE_DIR"
echo "---"

python3 -u summarize_simulations.py \
    --data-dir "$DATA_DIR" \
    --output-dir results \
    --force \
    --build-only
