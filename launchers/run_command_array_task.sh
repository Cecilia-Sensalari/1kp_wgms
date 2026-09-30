#!/bin/bash
# SLURM array worker: runs exactly one line of a generated command file
# (init_commands.sh / paralogs_ks_commands.sh / orthologs_ks_commands.sh) per array task,
# selected by SLURM_ARRAY_TASK_ID. Not meant to be submitted directly - use
# submit_command_array.sh, which computes the array size and invokes this via sbatch.
#
# Each line in these files is a self-contained "( cd <dir> && cmd1 && cmd2 ... )" subshell
# already produced by prepare_branch_ksrates_configs.py's write_command_file(), so running it
# via a plain `eval` is safe - no extra parsing/splitting needed here.

#SBATCH --job-name=cmd_array
#SBATCH --output=array_logs/%x_%A_%a.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00

set -euo pipefail

COMMAND_FILE=$1
if [ -z "${COMMAND_FILE:-}" ]; then
    echo "ERROR: no command file given. Submit via submit_command_array.sh, not directly." >&2
    exit 1
fi
if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID not set - this script must be run as a SLURM array job." >&2
    exit 1
fi

# Skip the shebang, #SBATCH directives, comments, "set -euo pipefail", and blank lines - only
# the actual "( cd ... && ... )" command lines remain, one per real unit of work.
LINE=$(grep -v '^#' "$COMMAND_FILE" | grep -v '^set -euo pipefail$' | grep -v '^[[:space:]]*$' | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [ -z "$LINE" ]; then
    echo "ERROR: no command found at line index $SLURM_ARRAY_TASK_ID in $COMMAND_FILE" >&2
    echo "(the file may have fewer real command lines than the array size submitted)" >&2
    exit 1
fi

echo "[array task $SLURM_ARRAY_TASK_ID] $LINE"
eval "$LINE"
