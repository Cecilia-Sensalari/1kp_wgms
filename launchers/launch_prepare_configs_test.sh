#!/bin/bash
# Small-subset validation of prepare_branch_ksrates_configs.py before trusting it at the full
# ~2354-branch scale (see plan: "Validate the existing ortholog-pair logic on a small subset").
# Runs against 3 branches only (branch_subtrees_subset_test.tsv), with --run-init actually
# executing "ksrates init" for each generated config, and the new --ks-database-server-address-path wiring
# writing paralogs_ks_commands.sh alongside the existing orthologs_ks_commands.sh.

#SBATCH --job-name=prep_configs_test
#SBATCH --output=prepare_configs_test.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --partition=short

module load python/x86_64/3.11.4
source /group/esb/cesen/.venv/venv_python_11/bin/activate

SIF=/group/esb/cesen/1kp/software/ksrates_a6822fe_20261007.sif
KSRATES_RAW=/group/esb/cesen/1kp/ks_analysis/1kp/ksrates_raw
KS_DB=$KSRATES_RAW/ks_database/ks_database_server_address.txt

python /group/esb/cesen/1kp/code/1kp_wgms/run_ksrates_raw/prepare_branch_ksrates_configs.py \
    --datasets "$KSRATES_RAW/branch_subtrees_subset_test.tsv" \
    --species-metadata "$KSRATES_RAW/species_metadata.tsv" \
    --out-dir "$KSRATES_RAW/branch_setup_test" \
    --expert-config "$KSRATES_RAW/config_expert.txt" \
    --ksrates-command "singularity exec -B /group/esb/cesen/1kp --env PYTHONPATH=/group/esb/cesen/1kp/code/ksrates $SIF ksrates" \
    --paranome yes \
    --ks-database-server-address-path "$KS_DB" \
    --max-outgroups 4 \
    --run-init
