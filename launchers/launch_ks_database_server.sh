#!/bin/bash
# Launches the dedicated Ks database server (paralog and ortholog Ks data) of the 1kp analysis, via the
# "ksrates launch-ks-database-server" command (sqld is baked into the ksrates .sif image, see ksrates/Dockerfile).
# The container image, the ksrates code directory (if one is set) and the database location are read from the
# analysis' nextflow.config, so they are the same as those used by the wrapper:
#   ksrates_sif, ksrates_code_dir, and the folder of ks_database_server_address_path.
#
# Submit this ONCE (sbatch launch_ks_database_server.sh); it then runs indefinitely (the
# underlying command execs into sqld and never returns on success). If this cluster enforces a
# walltime limit on its default partition, add a "#SBATCH --time=..." line matching the longest
# available partition, and be prepared to resubmit (server data/keys persist across restarts,
# only the address file's host:port goes stale - see docs/paralog_ks_database.rst).

#SBATCH --job-name=server1kp
#SBATCH --output=ks_database_server.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G

ANALYSIS_DIR=/group/esb/cesen/1kp/ks_analysis/1kp

module load nextflow/x86_64/25.10.3

cd "$ANALYSIS_DIR"
CONFIG=$(nextflow config -flat)

# value of "params.<name> = <value>" in the resolved config, without quotes ('null' counts as not set)
param() {
    echo "$CONFIG" | sed -n "s/^params\.$1 = //p" | sed "s/^['\"]//; s/['\"]$//" | grep -v '^null$' || true
}

SIF=$(param ksrates_sif)
CODE_DIR=$(param ksrates_code_dir)
ADDRESS_PATH=$(param ks_database_server_address_path)

[ -n "$SIF" ] || { echo "ERROR: ksrates_sif is not set in $ANALYSIS_DIR/nextflow.config" >&2; exit 1; }
[ -n "$ADDRESS_PATH" ] || { echo "ERROR: ks_database_server_address_path is not set in $ANALYSIS_DIR/nextflow.config" >&2; exit 1; }
LOCATION=$(dirname "$ADDRESS_PATH")

PYTHONPATH_OPTION=()
[ -z "$CODE_DIR" ] || PYTHONPATH_OPTION=(--env "PYTHONPATH=$CODE_DIR")

echo "Image: $SIF"
echo "Database location: $LOCATION"
echo ""

singularity exec \
    -B /group/esb/cesen/1kp/ \
    "${PYTHONPATH_OPTION[@]}" \
    "$SIF" \
    ksrates launch-ks-database-server --location "$LOCATION"
