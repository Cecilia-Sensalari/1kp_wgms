#!/bin/bash
# Submits a generated command file (init_commands.sh / paralogs_ks_commands.sh /
# orthologs_ks_commands.sh, as produced by prepare_branch_ksrates_configs.py) as a SLURM job
# array, one array task per real command line - the way to actually run these at the real
# ~2354-branch scale, where running everything sequentially in one job would take far too long.
#
# Usage:
#   ./submit_command_array.sh <command_file.sh> [max_concurrent] [sbatch overrides...]
#
# Examples:
#   ./submit_command_array.sh /path/to/paralogs_ks_commands.sh
#   ./submit_command_array.sh /path/to/orthologs_ks_commands.sh 100
#   ./submit_command_array.sh /path/to/init_commands.sh 200 --cpus-per-task=1 --mem=2G --time=00:30:00
#
# max_concurrent (default 25) caps how many array tasks run simultaneously (SLURM's
# --array=1-N%C syntax) - important here since thousands of tasks all hitting the cluster (and
# the shared paralog Ks sqld server) at once would be its own problem, separate from the
# per-species/per-pair deduplication that already makes parallelizing this safe at all. This
# throttle can be changed on an already-submitted job without cancelling/resubmitting it:
#   scontrol update JobId=<jobid> ArrayTaskThrottle=<new_limit>
# Any extra arguments are passed straight through to sbatch, overriding
# run_command_array_task.sbatch's own #SBATCH defaults (e.g. to give paralogs-ks more memory
# than init needs).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMMAND_FILE=${1:?"Usage: $0 <command_file.sh> [max_concurrent] [sbatch overrides...]"}
COMMAND_FILE=$(realpath "$COMMAND_FILE")
MAX_CONCURRENT=${2:-25}
shift $(( $# >= 2 ? 2 : 1 ))
EXTRA_SBATCH_ARGS=("$@")

if [ ! -f "$COMMAND_FILE" ]; then
    echo "ERROR: command file not found: $COMMAND_FILE" >&2
    exit 1
fi

N=$(grep -v '^#' "$COMMAND_FILE" | grep -v '^set -euo pipefail$' | grep -cv '^[[:space:]]*$')
if [ "$N" -eq 0 ]; then
    echo "ERROR: no command lines found in $COMMAND_FILE - nothing to submit." >&2
    exit 1
fi

JOB_NAME=$(basename "$COMMAND_FILE" .sh)
LOG_DIR="$(dirname "$COMMAND_FILE")/array_logs"
mkdir -p "$LOG_DIR"

echo "Submitting $N tasks from $COMMAND_FILE as a job array (max $MAX_CONCURRENT concurrent)..."
sbatch \
    --job-name="$JOB_NAME" \
    --output="$LOG_DIR/%x_%A_%a.log" \
    --array="1-${N}%${MAX_CONCURRENT}" \
    "${EXTRA_SBATCH_ARGS[@]}" \
    "$SCRIPT_DIR/run_command_array_task.sbatch" \
    "$COMMAND_FILE"
