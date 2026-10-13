#!/bin/bash
#
# BUSCO completeness of the proteins of each genome species, measured as the 1KP paper did for its
# transcriptomes, so that the genome species can get a quality score comparable with the 1KP table column
# "% BUSCOs (complete+fragmented)".
#
# What the 1KP paper did (Leebens-Mack et al. 2019, Nature 574:679, Methods and Extended Data Fig. 3; and the
# 1KP data paper, GigaScience 8:giz126): BUSCO (an early version; "v1" is only reported second-hand and could not be confirmed) with default settings and the eukaryote (eukaryota_odb9) and
# embryophyte (embryophyta_odb9) conserved gene sets, run on the proteins predicted from each assembly. The
# percentage reported is "complete + fragmented" = 100% - missing. Samples with more than 42.5% missing
# eukaryota BUSCOs were removed. The single table column used here matches the eukaryota_odb9 set (BUSCO
# values of 30-98% for all kinds of algae and land plants, which the much larger embryophyta set would not give).
#
# Here: BUSCO 3.0.2b (cluster module), a version that still reads the odb9 gene sets (BUSCO v4+ cannot), in
# protein mode with default settings on the primary-transcript CDS of each genome translated to proteins.
# The original version is not available, so small differences are possible.
#
# One array task per species of the genome manifest (rows 2..N, so --array=1-31 for 31 species).
# Before running, put the eukaryota_odb9 lineage dataset (unpacked from the archived odb9 files, no longer on
# busco.ezlab.org) in OUT_DIR/eukaryota_odb9. The BUSCO and Python modules are loaded below.
#
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem 8G
#SBATCH -t 0-12:00
#SBATCH --array=1-31
#SBATCH -o /group/esb/cesen/1kp/source_data/2.genomes/1kp/genome_busco/logs/slurm.%x.%A_%a.out
#SBATCH -e /group/esb/cesen/1kp/source_data/2.genomes/1kp/genome_busco/logs/slurm.%x.%A_%a.err

set -euo pipefail

MANIFEST=/group/esb/cesen/1kp/source_data/2.genomes/1kp/genome_species_manifest.tsv
OUT_DIR=/group/esb/cesen/1kp/source_data/2.genomes/1kp/genome_busco
LINEAGE=${LINEAGE:-/group/esb/cesen/1kp/software/busco_lineages/eukaryota_odb9}
SCRIPT_DIR=/group/esb/cesen/1kp/code/1kp_wgms/others
PYTHON_MODULE=python/x86_64/3.11.4
BUSCO_MODULE=busco/x86_64/3.0.2b
THREADS=${SLURM_CPUS_PER_TASK:-4}

row=$(( ${SLURM_ARRAY_TASK_ID:?run it as an array job} + 1 ))
species=$(awk -F'\t' -v r="$row" 'NR == r {print $1}' "$MANIFEST")
cds=$(awk -F'\t' -v r="$row" 'NR == r {print $3}' "$MANIFEST")
[ -n "$species" ] && [ -r "$cds" ] || { echo "ERROR: no species or unreadable CDS file in manifest row $row" >&2; exit 1; }
[ -d "$LINEAGE" ] || { echo "ERROR: lineage dataset [$LINEAGE] not found" >&2; exit 1; }

# Proteins: OUT_DIR/translated_fasta/<CDS file name>.translated; translated here only if not there yet.
# The translation runs in a subshell so that its Python module does not stay loaded: BUSCO needs the Python
# that comes with its own module (it imports pipebricks, which a different Python does not have).
cds_name=$(basename "$cds")
proteins="$OUT_DIR/translated_fasta/$cds_name.translated"
if [ ! -s "$proteins" ]; then
    ( module load "$PYTHON_MODULE" && python3 "$SCRIPT_DIR/translate_cds_to_protein.py" "$cds" "$proteins" )
fi

module purge
module load "$BUSCO_MODULE"
python -c "import pipebricks" || { echo "ERROR: Python [$(command -v python)] cannot import pipebricks after loading $BUSCO_MODULE" >&2; exit 1; }

# BUSCO v3 writes run_<name>/ in the current directory and refuses to overwrite it
work="$OUT_DIR/output_busco/$species"
mkdir -p "$work"
cd "$work"
[ ! -e "run_$species" ] || { echo "ERROR: [$work/run_$species] already exists, remove it to rerun" >&2; exit 1; }
run_BUSCO.py -i "$proteins" -o "$species" -l "$LINEAGE" -m proteins -c "$THREADS"
