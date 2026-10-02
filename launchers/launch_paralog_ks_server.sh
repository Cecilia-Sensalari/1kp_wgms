#!/bin/bash
# Launches the dedicated paralog Ks database server for the real 1kp analysis, via the new
# "ksrates launch-paralog-ks-server" CLI command against the .sif built from the
# paralog_ks_db branch (has sqld baked in - see ksrates/Dockerfile).
#
# Submit this ONCE (sbatch launch_paralog_ks_server.sh); it then runs indefinitely (the
# underlying command execs into sqld and never returns on success). If this cluster enforces a
# walltime limit on its default partition, add a "#SBATCH --time=..." line matching the longest
# available partition, and be prepared to resubmit (server data/keys persist across restarts,
# only the address file's host:port goes stale - see docs/paralog_ks_database.rst).

#SBATCH --job-name=server1kp
#SBATCH --output=paralog_ks_database_server.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G

SIF=/group/esb/cesen/1kp/software/ksrates_paralog_ks_db.sif
LOCATION=/group/esb/cesen/1kp/ks_analysis/1kp/paralog_ks_database

singularity exec \
    -B /group/esb/cesen/1kp/ \
    "$SIF" \
    ksrates launch-paralog-ks-server --location "$LOCATION"
