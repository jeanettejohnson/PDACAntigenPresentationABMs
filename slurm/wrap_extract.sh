#!/bin/bash
# wrap_extract.sh -- one array task converts one simulation into derived/.
#
# Takes the repository root as $1. Reads its own simulation id from the array
# task id via the id list written by submit_extract.sh, because simulation ids
# run 1-77 with gaps (four runs were retried and took new ids), so the array
# index is not the simulation id.
#
# $1 must be passed in rather than self-located via ${BASH_SOURCE[0]}: sbatch
# copies the submitted script into its own spool directory and runs that copy,
# so self-location resolves to the spool dir. Same reason wrap_driver.sh takes
# its path as an argument.
set -eo pipefail

module load slurm

source /usr/local/packages/miniconda3/etc/profile.d/conda.sh
conda activate pcdl-260901

REPO="$1"
ID_LIST="$REPO/derived/.extract_ids"

# sed is 1-indexed and so is SLURM_ARRAY_TASK_ID, so they line up directly.
SIM="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ID_LIST")"
if [ -z "$SIM" ]; then
    echo "no simulation id at line $SLURM_ARRAY_TASK_ID of $ID_LIST" >&2
    exit 1
fi

echo "task $SLURM_ARRAY_TASK_ID -> simulation $SIM on $(hostname)"
cd "$REPO"
exec python -u scripts/extract_simulation_results.py --simulation "$SIM"
