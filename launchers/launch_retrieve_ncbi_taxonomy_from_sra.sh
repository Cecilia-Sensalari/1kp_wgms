#!/bin/bash
#SBATCH -N 1
#SBATCH --mem 4G
#SBATCH -t 00:30:00
#SBATCH -J 1kptaxo
#SBATCH -o 1kptaxo.o%j
#SBATCH -e 1kptaxo.e%j

cd /group/esb/cesen/1kp

python code/1kp_wgms/retrieve_ncbi_taxonomy_from_sra.py \
  --input source_data/1.species_dataset/1kp_paper_2019_suptab1_species.xlsx \
  --output source_data/1.species_dataset/1kp_paper_2019_suptab1_species_ncbi_taxonomy.csv \
  --email ceci.sensa93@gmail.com \
  &> source_data/1.species_dataset/1kp_paper_2019_suptab1_species_ncbi_taxonomy.log
