#!/bin/bash
# Launches the dedicated Ks database server (paralog and ortholog Ks data) for the real 1kp analysis, via the new
# "ksrates launch-ks-database-server" CLI command against the ksrates .sif image
# (has sqld baked in - see ksrates/Dockerfile), with the local ksrates checkout on PYTHONPATH so the
# unified paralog+ortholog Ks database code is used.
#
# Submit this ONCE (sbatch launch_paralog_ks_server.sh); it then runs indefinitely (the
# underlying command execs into sqld and never returns on success). If this cluster enforces a
# walltime limit on its default partition, add a "#SBATCH --time=..." line matching the longest
# available partition, and be prepared to resubmit (server data/keys persist across restarts,
# only the address file's host:port goes stale - see docs/paralog_ks_database.rst).

#SBATCH --job-name=server1kp
#SBATCH --output=ks_database_server.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G

SIF=/group/esb/cesen/1kp/software/ksrates_a6822fe_20261007.sif
LOCATION=/group/esb/cesen/1kp/ks_analysis/1kp/ks_database

singularity exec \
    -B /group/esb/cesen/1kp/ \
    --env PYTHONPATH=/group/esb/cesen/1kp/code/ksrates \
    "$SIF" \
    ksrates launch-ks-database-server --location "$LOCATION"
