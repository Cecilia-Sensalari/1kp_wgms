#!/bin/bash
# Runs the batch_ksrates wrapper (batch_ksrates/main.nf) for the 1kp analysis. Every setting (paths, parameters,
# resources, which stages run) lives in /group/esb/cesen/1kp/ks_analysis/1kp/nextflow.config, which Nextflow loads
# from the directory this script moves to; nothing is repeated here. Which stages run (run_init, execute_paralog_jobs,
# execute_ortholog_jobs, execute_ksrates_pipelines) is chosen there too, so that the config is the permanent record of how
# the analysis was run. Arguments given to this script (e.g. --run_init false) override the config for one launch,
# only meant for quick tests.
#   sbatch launch_batch_ksrates.sh
#
# This job only runs Nextflow itself, which submits the compute jobs to SLURM.
# The Ks database server must already be running (launch_ks_database_server.sh).
# The walltime must cover the whole run: raise it, or resubmit with -resume.
#
# Every launch starts the run from scratch. To continue the last run and skip the tasks that are done, pass -resume
# (sbatch launch_batch_ksrates.sh -resume), e.g. to retry only the failed jobs. Careful: Nextflow decides what is
# cached from the script of each task and its inputs, not from the content of the Python scripts or of the tables
# they read, so after changing those, launch without -resume. Resuming also fails (lock on the session) while the
# run it would resume is still going.

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
# -ansi-log false: in a file, the ANSI progress display would print the whole progress block again at every update.
# Progress at a glance: others/nextflow_progress.sh (per process, from the trace file), or run in a terminal.
nextflow run "$MAIN_NF" "$@" -ansi-log false
