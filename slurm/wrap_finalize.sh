#!/bin/bash
# wrap_finalize.sh -- runs once, after every array task has succeeded.
#
# Builds the combined snapshots and validates the whole archive. Submitted by
# submit_extract.sh with --dependency=afterok on the array, so it does not run
# at all if any simulation failed to convert. That is deliberate: combining a
# partial archive would produce a file that looks complete and is not, which is
# the failure mode this pipeline has produced more often than any other.
#
# Takes the repository root as $1, for the same reason wrap_extract.sh does:
# sbatch runs its own copy of this script out of the spool directory, so
# self-location would resolve there rather than to slurm/.
set -eo pipefail

module load slurm

source /usr/local/packages/miniconda3/etc/profile.d/conda.sh

REPO="$1"
VALIDATE_ARGS="${2:-}"   # e.g. --only-present, for a --limit run
cd "$REPO"

# Combine in the pcdl environment. The per-simulation files were written by the
# anndata version living there, and keeping one writer across the archive avoids
# the nullable-string encoding mismatch documented in combine_snapshots.py.
conda activate pcdl-260901
echo "=== combining snapshots ==="
python -u scripts/combine_snapshots.py

# Validate in the analysis environment instead -- checking the archive where it
# will actually be read is the only version of this check that means anything.
conda activate physicell-sim-260606
echo
echo "=== validating archive ==="
python -u scripts/validate_derived.py $VALIDATE_ARGS
