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
#   ./submit_spatial_pipeline.sh build <label> [data-dir-override]
#   ./submit_spatial_pipeline.sh sweep <label>
#
# <label> is a short tag for the dataset (e.g. htan_wellmixed, htan_geometries,
# imc_wellmixed, imc_spatial -- the same labels submit_driver.sh's SIM_NAMES
# uses) and drives three things:
#   - the data dir, by convention:
#       /local/projects-t3/fertig_pdacagmodel/<label>/PDACAntigenPresentationABMs/data/outputs/simulations
#     pass [data-dir-override] as a 3rd arg to the "build" call if your
#     dataset doesn't live at that conventional path
#   - the build stage's --output-dir: results_<label>/
#   - the sweep stage's config path: runs/<label>/config.yaml (requires the
#     build to have finished, and that config.yaml to already exist -- e.g.
#     via "Save config.yaml" in the Streamlit Configure page, with Directory
#     set to this dataset's data dir and Run output directory set to
#     runs/<label>)
#
# Keeping build/sweep as two separate jobs (rather than one long job) means a
# sweep failure never forces rebuilding the spatial-stats cache from scratch,
# and per-dataset results_<label>/runs/<label> naming means different
# datasets' caches never collide with each other.
#
# NOTE: this script operates on the pdac-spatial-pipeline repo (a separate
# Python project living at /autofs/projects-t3/fertig_pdacagmodel/pdac-spatial-pipeline
# on thanos, not part of this Julia/PhysiCell repo). It's kept here in
# slurm/ alongside submit_driver.sh/wrap_driver.sh for convenience and so the
# submission pattern stays consistent, but wrap_spatial_build.sh and
# wrap_spatial_sweep.sh cd into that other repo to do their actual work.
#
# Mirrors the pattern in submit_driver.sh: this script computes the SLURM
# account once and hands off to a wrap_*.sh script that does the real work
# on the compute node.
set -euo pipefail

module load slurm

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="/autofs/projects-t3/fertig_pdacagmodel/pdac-spatial-pipeline"
mkdir -p "$SCRIPT_DIR/logs"

STAGE="${1:-}"
LABEL="${2:-}"
DATA_DIR_OVERRIDE="${3:-}"
if [[ "$STAGE" != "build" && "$STAGE" != "sweep" ]] || [[ -z "$LABEL" ]]; then
    echo "Usage: $0 {build|sweep} <label> [data-dir-override]" >&2
    echo "  e.g. $0 build htan_geometries" >&2
    echo "       $0 build my_dataset /custom/path/to/simulations" >&2
    echo "       $0 sweep htan_geometries" >&2
    exit 1
fi

if [[ -n "$DATA_DIR_OVERRIDE" ]]; then
    DATA_DIR="$DATA_DIR_OVERRIDE"
else
    DATA_DIR="/local/projects-t3/fertig_pdacagmodel/${LABEL}/PDACAntigenPresentationABMs/data/outputs/simulations"
fi
OUTPUT_DIR="results_${LABEL}"
CONFIG="runs/${LABEL}/config.yaml"

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
        --job-name="pdac-spatial-build-${LABEL}" \
        --output="$SCRIPT_DIR/logs/%x_%j.out" \
        --error="$SCRIPT_DIR/logs/%x_%j.err" \
        "$SCRIPT_DIR/wrap_spatial_build.sh" "$PIPELINE_DIR" "$DATA_DIR" "$OUTPUT_DIR"
else
    sbatch \
        --account="$PCMM_SLURM_ACCOUNT" \
        --time=24:00:00 \
        --cpus-per-task=4 \
        --mem=64G \
        --job-name="pdac-spatial-sweep-${LABEL}" \
        --output="$SCRIPT_DIR/logs/%x_%j.out" \
        --error="$SCRIPT_DIR/logs/%x_%j.err" \
        "$SCRIPT_DIR/wrap_spatial_sweep.sh" "$PIPELINE_DIR" "$CONFIG"
fi