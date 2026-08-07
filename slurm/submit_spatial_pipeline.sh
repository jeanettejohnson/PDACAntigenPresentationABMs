#!/bin/bash
# submit_spatial_pipeline.sh -- run this on the login node to submit the
# spatial-statistics pipeline to SLURM. Split into two independent jobs
# (build vs sweep) rather than one long job: the build step is the expensive,
# disk-heavy part (per-frame feature extraction + spatial-LDA fit across every
# discovered simulation), and the sweep step is the memory-heavy part (UMAP
# grid). Keeping them separate means a sweep failure never forces rebuilding
# the spatial-stats cache from scratch, and each stage can get its own
# resource request.
#
#   ./submit_spatial_pipeline.sh build     # build/refresh results/*.csv
#   ./submit_spatial_pipeline.sh sweep     # full UMAP grid -> runs/run1/
#                                           # (requires build to have finished,
#                                           #  and runs/run1/config.yaml to exist
#                                           #  -- e.g. via "Save config.yaml" in
#                                           #  the Streamlit Configure page)
#
# Mirrors the pattern in PDACAntigenPresentationABMs/slurm/submit_driver.sh:
# this script computes the SLURM account once and hands off to a wrap_*.sh
# script that does the real work on the compute node.
set -euo pipefail

module load slurm

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"
mkdir -p "$SCRIPT_DIR/logs"

STAGE="${1:-}"
if [[ "$STAGE" != "build" && "$STAGE" != "sweep" ]]; then
    echo "Usage: $0 {build|sweep}" >&2
    exit 1
fi

# Same lookup as submit_driver.sh -- computed once here so the job doesn't
# need to re-query sacctmgr itself.
export PCMM_SLURM_ACCOUNT=$(sacctmgr -n -P show assoc user="$(whoami)" format=account | head -1)

# --mem / --cpus-per-task are starting points, not measured limits -- check
# what your account can actually request before relying on these:
#   sinfo -o "%P %m %c"          # memory/cpus per partition
#   sacctmgr show assoc user=$(whoami) format=account,partition,maxsubmit
if [[ "$STAGE" == "build" ]]; then
    sbatch \
        --account="$PCMM_SLURM_ACCOUNT" \
        --time=12:00:00 \
        --cpus-per-task=4 \
        --mem=64G \
        --job-name="pdac-spatial-build" \
        --output="$SCRIPT_DIR/logs/%x_%j.out" \
        --error="$SCRIPT_DIR/logs/%x_%j.err" \
        "$SCRIPT_DIR/wrap_spatial_build.sh" "$PIPELINE_DIR"
else
    sbatch \
        --account="$PCMM_SLURM_ACCOUNT" \
        --time=24:00:00 \
        --cpus-per-task=4 \
        --mem=64G \
        --job-name="pdac-spatial-sweep" \
        --output="$SCRIPT_DIR/logs/%x_%j.out" \
        --error="$SCRIPT_DIR/logs/%x_%j.err" \
        "$SCRIPT_DIR/wrap_spatial_sweep.sh" "$PIPELINE_DIR"
fi