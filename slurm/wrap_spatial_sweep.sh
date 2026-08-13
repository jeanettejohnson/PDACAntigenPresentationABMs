#!/bin/bash
# wrap_spatial_sweep.sh -- generic job script; this is what actually runs on
# the compute node for the "sweep" stage.
#   $1 = pipeline directory (contains sweep.py)
#   $2 = path to config.yaml, relative to pipeline dir or absolute
#        (e.g. runs/htan_geometries/config.yaml)
#
# Only submit this stage after the build job has completed successfully --
# it reuses the results_<label>/ cache rather than recomputing it. Requires
# that config.yaml to already exist (e.g. from clicking "Save config.yaml" in
# the Streamlit Configure page, with Directory set to this dataset's data dir
# and Run output directory set to runs/<label>) since sweep.py takes the
# config path as its only positional argument, not a --data-dir flag.
#
# No .venv here -- there isn't one. Same interpreter resolution as
# wrap_spatial_build.sh: `module load python/3.11.4` is what actually works.
set -eo pipefail

module load slurm
module load python/3.11.4

PIPELINE_DIR="$1"
CONFIG="$2"

cd "$PIPELINE_DIR"

if [[ ! -f "$CONFIG" ]]; then
    echo "Missing $CONFIG -- save a config from the Streamlit Configure page first." >&2
    exit 1
fi

echo "Host: $(hostname)"
echo "python3: $(which python3)"
echo "Config: $CONFIG"
echo "--- memory available on this node ---"
free -h
echo "---"

python3 -u sweep.py "$CONFIG"
