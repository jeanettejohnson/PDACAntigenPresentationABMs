#!/bin/bash
#SBATCH --job-name=pdac-extract-submit
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=slurm/logs/submit_%j.out
#SBATCH --error=slurm/logs/submit_%j.err
#
# submit_extract.sh -- convert every completed simulation into derived/, then
# combine and validate. One command for the whole pipeline.
#
#   ./slurm/submit_extract.sh          # run here; returns as soon as it has queued
#   sbatch slurm/submit_extract.sh     # or queue the submitter itself
#   ./slurm/submit_extract.sh --dry    # print what would be submitted
#
#   CONCURRENCY=4 ./slurm/submit_extract.sh    # fewer tasks at once
#
# To change the cap on an array that is already queued or running:
#
#   scontrol update jobid=<array id> arraytaskthrottle=4
#
# The #SBATCH block above is a comment to bash, so the same file works either
# way. Running it directly is usually what you want: it exits in seconds and
# prints the job ids. Submitting it is for when the login node is not somewhere
# you want to leave a shell, or when this is itself called from something else.
#
# Two stages:
#
#   1. an array, one task per simulation, capped at twelve concurrent
#   2. a single finalize task, --dependency=afterok on the array
#
# afterok means the finalize does not run if any conversion failed. A partial
# archive that has been combined and reported as valid is worse than an obvious
# failure, and this pipeline has produced silent partial success six times.
set -eo pipefail

module load slurm

# How many array tasks run at once. Twelve by default: the work is ~85% NFS read
# wait, so beyond that tasks contend on a shared mount rather than going faster,
# and each needs about 5 GB of scratch for its intermediate parts. Lower it when
# the filesystem is busy or somebody else needs the nodes; raising it is unlikely
# to help and will be felt by everyone else on the mount.
CONCURRENCY="${CONCURRENCY:-12}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/.." && pwd)"
ID_LIST="$REPO/derived/.extract_ids"

mkdir -p "$REPO/derived" "$SCRIPT_DIR/logs"

# Ids come from the model manager's status table, not from a range: they run
# 1-77 with gaps, and only Completed runs have output worth converting.
source /usr/local/packages/miniconda3/etc/profile.d/conda.sh
conda activate pcdl-260901
python - "$REPO" > "$ID_LIST.new" <<'PY'
import sqlite3, sys
from pathlib import Path
db = Path(sys.argv[1]) / "data" / "pcmm.db"
with sqlite3.connect(db) as con:
    for (sid,) in con.execute(
        "SELECT s.simulation_id FROM simulations s "
        "JOIN status_codes c USING(status_code_id) "
        "WHERE c.status_code = 'Completed' ORDER BY s.simulation_id"
    ):
        print(sid)
PY

N="$(wc -l < "$ID_LIST.new")"
echo "$N completed simulations -> $ID_LIST ($CONCURRENCY at a time)"

# Job names carry the simulation type so the guard below, and squeue, can tell
# one clone's array from another's. Four clones now run this concurrently; a
# bare "pdac-extract" made the first submission block the other three. The type
# comes from the resolver, not from the directory name -- same rule as the
# metadata itself.
TAG="$(cd "$REPO" && python -c "
import sys; sys.path.insert(0, '.')
from scripts.resolve_samples import sim_type
from pathlib import Path
print(sim_type(Path('.')))
" 2>/dev/null)"
if [ -z "$TAG" ]; then
    echo "could not determine the simulation type for $REPO" >&2
    exit 1
fi
echo "simulation type: $TAG"

# This cluster refuses submissions without an account. submit_driver.sh resolves
# it the same way -- first association for this user -- rather than hardcoding it.
ACCOUNT="$(sacctmgr -n -P show assoc user="$(whoami)" format=account | head -1)"
if [ -z "$ACCOUNT" ]; then
    echo "no slurm account found for $(whoami)" >&2
    exit 1
fi
echo "billing to account: $ACCOUNT"

if [ "${1:-}" = "--dry" ]; then
    rm -f "$ID_LIST.new"
    echo "would submit:"
    echo "  stage 1  sbatch --job-name=pdac-extract-$TAG --array=1-$N%$CONCURRENCY wrap_extract.sh"
    echo "  stage 2  sbatch --dependency=afterok:<stage1> wrap_finalize.sh"
    exit 0
fi

# Publish the id list only now. Array tasks read it by line number for the whole
# life of the job, so rewriting it under a running array would repoint tasks at
# different simulations -- silently, and only when the id set had changed. Refuse
# instead of racing.
if squeue -h -n "pdac-extract-$TAG" -u "$(whoami)" 2>/dev/null | grep -q .; then
    rm -f "$ID_LIST.new"
    echo "a pdac-extract-$TAG array is already queued or running; not submitting another." >&2
    echo "wait for it to finish, or scancel it first." >&2
    exit 1
fi
mv "$ID_LIST.new" "$ID_LIST"

# --- stage 1: one task per simulation -------------------------------------
ARRAY_ID="$(sbatch --parsable \
    --account="$ACCOUNT" \
    --job-name="pdac-extract-$TAG" \
    --array="1-${N}%${CONCURRENCY}" \
    --cpus-per-task=2 \
    --mem=16G \
    --time=02:00:00 \
    --output="$SCRIPT_DIR/logs/extract_%A_%a.out" \
    --error="$SCRIPT_DIR/logs/extract_%A_%a.err" \
    "$SCRIPT_DIR/wrap_extract.sh" "$REPO")"
echo "stage 1  array   $ARRAY_ID  ($N tasks, $CONCURRENCY concurrent)"

# --- stage 2: combine and validate, only if every task succeeded ----------
FINAL_ID="$(sbatch --parsable \
    --account="$ACCOUNT" \
    --job-name="pdac-finalize-$TAG" \
    --dependency="afterok:${ARRAY_ID}" \
    --cpus-per-task=2 \
    --mem=32G \
    --time=01:00:00 \
    --output="$SCRIPT_DIR/logs/finalize_%j.out" \
    --error="$SCRIPT_DIR/logs/finalize_%j.err" \
    "$SCRIPT_DIR/wrap_finalize.sh" "$REPO")"
echo "stage 2  finalize $FINAL_ID  (after array succeeds)"
echo
echo "watch:   squeue -j $ARRAY_ID,$FINAL_ID"
echo "results: $SCRIPT_DIR/logs/finalize_${FINAL_ID}.out"
