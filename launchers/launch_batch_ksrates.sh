#!/bin/bash
# Runs the batch_ksrates wrapper (batch_ksrates/main.nf) for the 1kp analysis. Every setting (paths, parameters,
# resources, which stages run) lives in /group/esb/cesen/1kp/ks_analysis/1kp/nextflow.config, which Nextflow loads
# from the directory this script moves to; nothing is repeated here. The stages and any parameter can still be
# chosen for one run by passing them as arguments, which override the config, e.g.
#   sbatch launch_batch_ksrates.sh                           (generate species metadata, datasets and configs only)
#   sbatch launch_batch_ksrates.sh --run_init
#   sbatch launch_batch_ksrates.sh --execute_ortholog_jobs --execute_paralog_jobs
#
# This job only runs Nextflow itself, which submits the compute jobs to SLURM; -resume lets a rerun skip what is done.
# The Ks database server must already be running (launch_ks_database_server.sh).
# The walltime must cover the whole run: raise it, or resubmit (the run continues with -resume).

#SBATCH --job-name=batch_ksrates
#SBATCH -p all
#SBATCH -c 1
#SBATCH --mem 4G
#SBATCH -o slurm.%x.%j.out
#SBATCH -e slurm.%x.%j.err

ANALYSIS_DIR=/group/esb/cesen/1kp/ks_analysis/1kp
MAIN_NF=/group/esb/cesen/1kp/code/1kp_wgms/batch_ksrates/main.nf

module load nextflow/x86_64/25.10.3

cd "$ANALYSIS_DIR"
nextflow run "$MAIN_NF" -resume "$@"
