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
# only the address file's host:port goes stale - see docs/ks_database.rst).

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

# Optional: keep the data folder of the server on the local disk of one node (see ks_database_local_data_dir in the
# nextflow.config). sqld swaps files with renameat2(RENAME_EXCHANGE) when it compacts its log, and NFS refuses that
# (EINVAL), which crashed the server. The data folder in the location becomes a symlink to the local folder, and the
# address file and the keys stay in the location, so no job needs another path. The server then has to start on that
# node (sbatch --nodelist=<node> ...), where the local folder is.
LOCAL_DATA_DIR=$(param ks_database_local_data_dir)
SERVER_NODE=$(param ks_database_server_node)

# Optional: port sqld listens on (ks_database_server_port in the nextflow.config). Unset = ksrates' default, 8080.
# Two servers cannot share a port on one node, so give each analysis its own port when they run on the same node.
PORT=$(param ks_database_server_port)
PORT_OPTION=()
[ -z "$PORT" ] || PORT_OPTION=(--port "$PORT")
LOCAL_BIND=()
if [ -n "$LOCAL_DATA_DIR" ]; then
    [ -z "$SERVER_NODE" ] || [ "$(hostname -s)" = "$SERVER_NODE" ] || \
        { echo "ERROR: ks_database_server_node is $SERVER_NODE but this job runs on $(hostname -s): submit with sbatch --nodelist=$SERVER_NODE" >&2; exit 1; }
    DATA_LINK="$LOCATION/ks_database_sqld_data"
    if [ -L "$DATA_LINK" ]; then
        [ "$(readlink "$DATA_LINK")" = "$LOCAL_DATA_DIR" ] || \
            { echo "ERROR: $DATA_LINK points to $(readlink "$DATA_LINK"), not to ks_database_local_data_dir ($LOCAL_DATA_DIR)" >&2; exit 1; }
        [ -d "$LOCAL_DATA_DIR" ] || \
            { echo "ERROR: $LOCAL_DATA_DIR does not exist on $(hostname -s): the data of the server is not on this node" >&2; exit 1; }
    elif [ -e "$DATA_LINK" ]; then
        echo "ERROR: $DATA_LINK is a real folder: copy it to $LOCAL_DATA_DIR on $(hostname -s) and move it aside first (see docs/how_to_launch_batch_wrapper.md)" >&2
        exit 1
    else
        mkdir -p "$LOCAL_DATA_DIR" && ln -s "$LOCAL_DATA_DIR" "$DATA_LINK"
    fi
    LOCAL_BIND=(-B "$LOCAL_DATA_DIR")
    echo "Data folder: $DATA_LINK -> $LOCAL_DATA_DIR (local disk of $(hostname -s))"
fi

PYTHONPATH_OPTION=()
[ -z "$CODE_DIR" ] || PYTHONPATH_OPTION=(--env "PYTHONPATH=$CODE_DIR")

echo "Image: $SIF"
echo "Database location: $LOCATION"
echo ""

singularity exec \
    -B /group/esb/cesen/1kp/ \
    "${LOCAL_BIND[@]}" \
    "${PYTHONPATH_OPTION[@]}" \
    "$SIF" \
    ksrates launch-ks-database-server --location "$LOCATION" "${PORT_OPTION[@]}"
