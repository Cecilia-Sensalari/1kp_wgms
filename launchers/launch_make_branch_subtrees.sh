#!/bin/bash
#
#SBATCH -p all # partition (queue)
#SBATCH -c 1 # number of cores
#SBATCH --mem 7G # memory pool for all cores
#SBATCH -t 0-01:00 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

cd /group/esb/cesen/1kp/code/1kp_wgms/run_ksrates_raw

module load python/x86_64/3.11.4
source /group/esb/cesen/.venv/venv_python_11/bin/activate

output_dir="/group/esb/cesen/1kp/ks_analysis/1kp/ksrates_raw"

python3 make_branch_subtrees.py \
  --max-a 4 \
  --max-older 2 \
  --older-clades 3 \
  --out $output_dir/branch_subtrees.tsv \
  &> $output_dir/make_branch_subtrees.log