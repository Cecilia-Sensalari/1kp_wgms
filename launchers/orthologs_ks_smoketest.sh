#!/usr/bin/env bash
#SBATCH --job-name=orthologs_ks_test
#SBATCH --output=orthologs_ks_test.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00

set -euo pipefail

( cd /group/esb/cesen/1kp/ks_analysis/1kp/ksrates_raw/branch_setup_test/ortholog_runs && singularity exec -B /group/esb/cesen/1kp /group/esb/cesen/1kp/software/ksrates_paralog_ks_db.sif ksrates orthologs-ks /group/esb/cesen/1kp/ks_analysis/1kp/ksrates_raw/branch_setup_test/configs/000002_internal/config_000002_internal_APTP.txt --expert /group/esb/cesen/1kp/ks_analysis/1kp/ksrates_raw/config_expert.txt APTP JKHA --n-threads 4 )
