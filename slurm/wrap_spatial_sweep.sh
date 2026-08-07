#!/bin/bash
# wrap_spatial_sweep.sh -- generic job script; this is what actually runs on
# the compute node for the "sweep" stage. Takes the pipeline directory
# (contains sweep.py, runs/run1/config.yaml) as $1.
#
# Only submit this stage after the build job has completed successfully --
# it reuses the results/ cache rather than recomputing it. Requires
# runs/run1/config.yaml to already exist (e.g. from clicking "Save
# config.yaml" in the Streamlit Configure page) since sweep.py takes the
# config path as its only positional argument, not a --data-dir flag.
#
# No .venv here -- there isn't one. Same interpreter resolution as
# wrap_spatial_build.sh: `module load python/3.11.4` is what actually works.
set -eo pipefail

module load slurm
module load python/3.11.4

PIPELINE_DIR="$1"
CONFIG="$PIPELINE_DIR/runs/run1/config.yaml"

cd "$PIPELINE_DIR"

if [[ ! -f "$CONFIG" ]]; then
    echo "Missing $CONFIG -- save a config from the Streamlit Configure page first." >&2
    exit 1
fi

echo "Host: $(hostname)"
echo "python3: $(which python3)"
echo "--- memory available on this node ---"
free -h
echo "---"

python3 -u sweep.py "$CONFIG"
