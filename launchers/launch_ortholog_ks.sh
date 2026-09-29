#!/bin/bash
#
#SBATCH -p all # partition (queue)
#SBATCH -c 1 # number of cores
#SBATCH --mem 7G # memory pool for all cores
#SBATCH -t 0-48:00 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

module load python/x86_64/3.11.4
source /group/esb/cesen/.venv/venv_python_11/bin/activate

SIF=/group/esb/cesen/1kp/software/ksrates_paralog_ks_db.sif
KSRATES_RAW=/group/esb/cesen/1kp/ks_analysis/1kp/ksrates_raw
PARALOG_DB=$KSRATES_RAW/paralog_ks_database/paralog_ks_server_address.txt

python /group/esb/cesen/1kp/code/1kp_wgms/run_ksrates_raw/prepare_branch_ksrates_configs.py \
  --subtrees "$KSRATES_RAW/branch_subtrees.tsv" \
  --species-metadata "$KSRATES_RAW/species_metadata.tsv" \
  --out-dir "$KSRATES_RAW/branch_setup" \
  --expert-config "$KSRATES_RAW/config_expert.txt" \
  --ksrates-command "singularity exec -B /group/esb/cesen/1kp $SIF ksrates" \
  --paranome yes \
  --paralog-database "$PARALOG_DB" \
  --max-outgroups 4
  # --run-init intentionally omitted here: running "ksrates init" for all ~2354 branches in one
  # job would be far too slow/long-running for a single SBATCH job. Once configs look right,
  # run init_commands.sh split across many parallel jobs instead (same reasoning as
  # paralogs_ks_commands.sh/orthologs_ks_commands.sh - internally independent, safe to
  # parallelize, no duplicate targets by construction).