#!/bin/bash
# submit_driver.sh -- run this on the login node to submit one or more of the
# numbered simulation drivers to SLURM. Pass simulation numbers as positional
# arguments, or run with no arguments to be prompted interactively.
#
#   ./submit_driver.sh 1 3      # submit htan_wellmixed and imc_wellmixed
#   ./submit_driver.sh          # interactive prompt
set -euo pipefail

module load slurm

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

declare -A SIM_NAMES=(
    [1]="htan_wellmixed"
    [2]="htan_geometries"
    [3]="imc_wellmixed"
    [4]="imc_spatial"
)

if [ "$#" -eq 0 ]; then
    echo "Available simulations:"
    for n in "${!SIM_NAMES[@]}"; do echo "  $n) ${SIM_NAMES[$n]}"; done | sort -n
    read -rp "Enter simulation number(s) to submit (space-separated): " -a REPLY_ARR
    set -- "${REPLY_ARR[@]}"
fi

# Exported so the driver job (and the julia process it launches) sees the
# same account without re-querying sacctmgr -- slurm_common.jl checks this
# env var first and only falls back to its own sacctmgr lookup when a script
# is run standalone (bypassing this submit script).
export PCMM_SLURM_ACCOUNT=$(sacctmgr -n -P show assoc user="$(whoami)" format=account | head -1)

for n in "$@"; do
    if [[ -z "${SIM_NAMES[$n]:-}" ]]; then
        echo "Unknown simulation number: $n" >&2
        exit 1
    fi
    echo "Submitting simulation $n (${SIM_NAMES[$n]})..."
    sbatch \
        --account="$PCMM_SLURM_ACCOUNT" \
        --time=7-00:00:00 \
        --cpus-per-task=1 \
        --mem=5G \
        --job-name="pcmm-driver-${SIM_NAMES[$n]}" \
        --output="$SCRIPT_DIR/logs/%x_%j.out" \
        --error="$SCRIPT_DIR/logs/%x_%j.err" \
        "$SCRIPT_DIR/wrap_driver.sh" "${SIM_NAMES[$n]}" "$SCRIPT_DIR"
done
