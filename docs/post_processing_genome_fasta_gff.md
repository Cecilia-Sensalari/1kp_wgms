# Post-processing of genome FASTA and GFF files for the collinearity analysis

Processing done on 2026-10-12. Species ids are those of `genome_species_manifest.tsv`. For each processed species
the original files and the logs are in the subdirectory with its name, next to this file; the processed files are in
`genome_species/`. Species that are not described in a section below were not processed.

## Summary

| Species | Processing | Feature / attribute for the collinearity analysis |
|---|---|---|
| `Chabr_v0.1` | CDS renamed to locus tag (longest CDS per locus tag kept); 103 split `gene` lines merged | `gene` / `locus_tag` |
| `Cyame_v1.0` | CDS renamed to locus tag (no duplicates); GFF unchanged (no split genes) | `mRNA` / `locus_tag` |
| `Solly_iTAGv2.3`, `Micpu_v3.0` | FASTA headers reduced to the first word | `mRNA` / `Name` |
| `Pinta_v2.0` | GTF converted to GFF3 | `gene` / `ID` |
| `Chlva_v1.0` | none | `CDS` / `ID` (set manually, see below) |
| `Klefl_v1.0` | none (v1.0 CDS) | none: no GFF available, so no anchor-pair (collinearity) analysis |
| `Betvu_v1.1`, `Elagu_v2.0`, `Nelnu_v1.0` | none | `mRNA` / `ID` |
| `Ambtr_v1.0.27`, `Aquco_v1.1`, `Arath_TAIR10`, `Carpa_ASGPBv0.4`, `Chlre_v5.5`, `Eucgr_v1.1`, `Manes_v4.1`, `Mimgu_v2.0`, `Musac_v1.0`, `Orysa_v7.0`, `Ostlu_v2.0`, `Phavu_v1.0`, `Phypa_v3.0`, `Poptr_v3.0`, `Prupe_v1.0`, `Selmo_v1.0`, `Sorbi_v2.1`, `Spipo_v2`, `Theca_v1.1`, `Vitvi_Genoscope.12X`, `Volca_v2.0` | none | `mRNA` / `Name` |

The pairs are the best candidate of `scan_gff_feature_attribute.py` for each species (result:
`source_data/2.genomes/1kp/gff_feature_attribute_scan.tsv`), except `Chlva_v1.0`. For `Cyame_v1.0`, `Pinta_v2.0`,
`Elagu_v2.0`, `Nelnu_v1.0`, `Selmo_v1.0` and `Spipo_v2` the scan reports `TIE`: a second pair is equally valid (one GFF
line per gene, all FASTA names found), and the first candidate (`mRNA` before `gene`) was taken.
`Chlva_v1.0` has only `CDS` lines in its GFF, one per gene (9791, the number of FASTA names), so the scan, which
only considers features that span the whole transcript, reports no candidate; `CDS` / `ID` was chosen manually. Its
`CDS` lines cover the coding sequence only, not the UTRs.

## Why

The collinearity analysis (ksrates, via wgd's `gff_parser`) makes one gene-list entry for every GFF line of the
chosen feature (`gff_feature`), with the whole value of the chosen attribute (`gff_attribute`) as the gene ID and
the coordinates of that same line. The attribute values must be the sequence names of the CDS FASTA file, and
each name must be on exactly one GFF line of that feature (otherwise a gene enters the list several times).
The original files of the species processed below did not satisfy this or could not be read by the pipeline (the
NCBI files of `Chabr_v0.1` and `Cyame_v1.0`, the FASTA headers of `Solly_iTAGv2.3` and `Micpu_v3.0`, the annotation
format of `Pinta_v2.0`), so they were processed as described. The original files are kept, unchanged, in the species
subdirectory (`*_ORIGINAL.*`).

## Chabr_v0.1 (*Chara braunii*)

NCBI assembly GCA_003427395.1 Cbr_1.0, strain S276.

### Input files (originals, in `Chabr_v0.1/`)

| File | Content | md5 |
|---|---|---|
| `GCA_003427395.1_Cbr_1.0_cds_from_genomic_ORIGINAL.fna` | 35667 CDS records, names like `lcl\|BFEA01000001.1_cds_GBG58599.1_1 [locus_tag=CBR_g1] ...` | `694771c6776af12bf4b727651c516e64` |
| `GCA_003427395.1_Cbr_1.0_genomic_ORIGINAL.gff` | GFF3 written by NCBI annotwriter (source column DDBJ), 35701 `gene` lines | `12acf2ab8742c6a99e2f3715a525a023` |

### Steps

**1. Rename the CDS sequences to their `locus_tag`** (`>lcl|..._GBG58599.1_1 [locus_tag=CBR_g1] ...` becomes
`>CBR_g1`) so that the names are the values the GFF uses to identify the genes.
- 842 of the 34721 locus tags have several CDS records (alternative isoforms with different `protein_id`
  and length; the sequences differ in 835 of the 842 cases, the other 7 are identical). Sequence names must be
  unique, so one record was kept per locus tag: the **longest** CDS (the first one if equally long), i.e. one
  representative transcript per gene. 946 records were dropped.
- Result: 34721 sequences, one per locus tag.
- Script: `code/1kp_wgms/others/post_process_genome_fasta_gff/rename_fasta_to_locus_tag.sh`
- Logs: `GCA_003427395.1_Cbr_1.0_cds_from_genomic.locus_tag.fna.mapping.tsv` (new name, original name, length of
  the kept records) and `...locus_tag.fna.dropped.tsv` (dropped records and the record kept instead).

**2. Merge the split `gene` lines of the GFF.** 103 genes are split into several `gene` lines (382 lines in
total) with the same ID, Name and `locus_tag` but different coordinates (NCBI `is_ordered=true` fragments).
Each was replaced by one `gene` line from the lowest start to the highest end of its fragments (the attributes
of the first fragment are kept and `;merged_fragments=N` is appended). Fragments are merged only if they are
on the same contig and strand, which was the case for all 103. For all 103, the merged span is identical to the
span of the `mRNA` line of the same locus tag; the largest merged gene spans 417819 bases (CBR_g30033). All other
lines are unchanged. `gene` lines: 35701 before, 35422 after.
- Script: `code/1kp_wgms/others/post_process_genome_fasta_gff/merge_gff_gene_fragments.sh`
- Log: `GCA_003427395.1_Cbr_1.0_genomic.merged_genes.gff.merge_log.tsv` (per gene: fragments, merged span,
  `mRNA` span, status).

**3. Choice of the feature and attribute** with `code/1kp_wgms/others/post_process_genome_fasta_gff/scan_gff_feature_attribute.py`, which
tests every (feature, attribute) pair of the GFF against the FASTA names and keeps only pairs that identify each
gene exactly once and whose feature spans the whole transcript (`mRNA`, `gene` or `transcript`).
- `gene` / `locus_tag`: 34719 of the 34721 FASTA names found (99.99%), one line per gene. **Chosen.**
  The two names without a `gene` line are `CBR_g146333` and `CBR_g46342` because they are labeled as `pseudogene`.
- `gene` / `Name`: 99.7%. For 99 genes `Name` is a gene symbol (e.g. `CNGC2`) and not the locus tag.
- `mRNA` / `locus_tag`: not valid. Step 1 reduced only the FASTA to one sequence per gene; the original GFF still has
  one `mRNA` line per isoform, also for the isoforms whose CDS was dropped. So after the renaming the same FASTA name is
  found on several `mRNA` lines (e.g. `CBR_g1047` on `rna-CBR_g1047` and `rna-CBR_g1047-2`, both `149067..177435`), and
  the parser would list the gene once per line. The 842 multi-isoform genes have 944 `mRNA` lines beyond one per gene.
  (`gene` lines do not have this problem: one per gene after step 2.)

### Resulting files and settings

| File (in `genome_species/`) | md5 |
|---|---|
| `GCA_003427395.1_Cbr_1.0_cds_from_genomic.locus_tag.fna` | `22ad023a5de054a63cc092229a99bda8` |
| `GCA_003427395.1_Cbr_1.0_genomic.merged_genes.gff` | `d180732ae62bffce69adf0ae5e2fbd06` |

`cds_filepath` and `gff_filepath` of `Chabr_v0.1` in `genome_species_manifest.tsv` point to these two files.
For the collinearity analysis use `gff_feature = gene` and `gff_attribute = locus_tag` (ksrates passes an
attribute other than ID, Name and Parent unchanged, with a warning that it is not a common one).

### Reproduction

```bash
cd /path/to/1kp/source_data/2.genomes/1kp/genome_species
O=post_processing_fasta_gff_logs/Chabr_v0.1
bash /path/to/1kp/code/1kp_wgms/others/post_process_genome_fasta_gff/rename_fasta_to_locus_tag.sh \
    $O/GCA_003427395.1_Cbr_1.0_cds_from_genomic_ORIGINAL.fna \
    GCA_003427395.1_Cbr_1.0_cds_from_genomic.locus_tag.fna
bash /path/to/1kp/code/1kp_wgms/others/post_process_genome_fasta_gff/merge_gff_gene_fragments.sh \
    $O/GCA_003427395.1_Cbr_1.0_genomic_ORIGINAL.gff \
    GCA_003427395.1_Cbr_1.0_genomic.merged_genes.gff
```

## Cyame_v1.0 (*Cyanidioschyzon merolae*)

NCBI RefSeq assembly ASM9120v1 (GCF_000091205.1), per the manifest.

### Input files (in `Cyame_v1.0/` and `genome_species/`)

| File | Content | md5 |
|---|---|---|
| `Cyame_v1.0/GCF_000091205.1_ASM9120v1_cds_from_genomic.fna` | 4803 CDS records, names like `lcl\|NC_010127.1_cds_XP_005535014.1_1 [locus_tag=CYME_CMA001C] [db_xref=GeneID:16992110] ...` | `79aa8baa1518adaa72c1e13d547317a4` |
| `GCF_000091205.1_ASM9120v1_genomic.gff` (not modified, so no `_ORIGINAL` copy) | GFF3 written by NCBI annotwriter (source column RefSeq), 5373 `gene` and 4803 `mRNA` lines | `c7e798e2c041af9f9d5086b6035aa234` |

### Steps

**1. Rename the CDS sequences to their `locus_tag`** (`>lcl|NC_010127.1_cds_XP_005535014.1_1 [locus_tag=CYME_CMA001C]
...` becomes `>CYME_CMA001C`), as for Chabr_v0.1.
- Every locus tag has exactly one CDS record, so all 4803 sequences were kept and none dropped (the
  `...dropped.tsv` log is empty).
- Result: 4803 sequences, one per locus tag.
- Script: `code/1kp_wgms/others/post_process_genome_fasta_gff/rename_fasta_to_locus_tag.sh`
- Logs in `Cyame_v1.0/`: `GCF_000091205.1_ASM9120v1_cds_from_genomic.locus_tag.fna.mapping.tsv` and
  `...locus_tag.fna.dropped.tsv` (empty).

**2. Check of the GFF for split `gene` lines**, with `code/1kp_wgms/others/post_process_genome_fasta_gff/merge_gff_gene_fragments.sh`, the
script used for Chabr_v0.1. No locus tag has more than one `gene` line, so nothing was merged and the output was
identical to the input (same md5): the GFF was **not changed** and the original file is used. Only the log (header
only, no merged gene) is kept, in `Cyame_v1.0/`: `GCF_000091205.1_ASM9120v1_genomic.merged_genes.gff.merge_log.tsv`.
The output copy, identical to the input, is not kept.

**3. Choice of the feature and attribute** with `code/1kp_wgms/others/post_process_genome_fasta_gff/scan_gff_feature_attribute.py`, against the
renamed FASTA. Three pairs find all 4803 FASTA names (coverage 1.000) with exactly one GFF line per gene, and
are therefore all valid (the scan reports `TIE`):
- `mRNA` / `locus_tag` (the transcript line; one `mRNA` per gene). **Chosen** (first candidate of the scan),
- `gene` / `Name` (the `Name` is the locus tag for every gene),
- `gene` / `locus_tag`.

The 5373 `gene` lines also include 570 non-coding genes (521 `misc_RNA`, 37 `tRNA`, 12 `rRNA`), which have no
CDS and are ignored.

### Resulting files and settings

| File (in `genome_species/`) | md5 |
|---|---|
| `GCF_000091205.1_ASM9120v1_cds_from_genomic.locus_tag.fna` | `0c9b51c0a37a4dfdee28601fd63310d6` |
| `GCF_000091205.1_ASM9120v1_genomic.gff` (unchanged) | `c7e798e2c041af9f9d5086b6035aa234` |

`cds_filepath` of `Cyame_v1.0` in `genome_species_manifest.tsv` points to the renamed FASTA; `gff_filepath` still
points to the original GFF. For the collinearity analysis use `gff_feature = mRNA` and `gff_attribute = locus_tag`.

### Reproduction

```bash
cd /path/to/1kp/source_data/2.genomes/1kp/genome_species
O=post_processing_fasta_gff_logs/Cyame_v1.0
bash /path/to/1kp/code/1kp_wgms/others/post_process_genome_fasta_gff/rename_fasta_to_locus_tag.sh \
    $O/GCF_000091205.1_ASM9120v1_cds_from_genomic.fna \
    GCF_000091205.1_ASM9120v1_cds_from_genomic.locus_tag.fna
# GFF check: expect 5373 gene lines in input and output and no status line listed. The output
# ($O/...merged_genes.gff) is identical to the input and can be deleted; only the log is kept.
bash /path/to/1kp/code/1kp_wgms/others/post_process_genome_fasta_gff/merge_gff_gene_fragments.sh \
    GCF_000091205.1_ASM9120v1_genomic.gff \
    $O/GCF_000091205.1_ASM9120v1_genomic.merged_genes.gff
```

## Solly_iTAGv2.3 (*Solanum lycopersicum*) and Micpu_v3.0 (*Micromonas pusilla*)

The CDS FASTA headers carry a long description after the gene name, e.g.
`>Solyc00g005000.2.1 genomic_reference:SL2.40ch00 gene_region:16437-18189 ...` and
`>19233 pacid=27338642 polypeptide=19233 locus=MicpuC2.e_gw1.8.246.1 ID=19233.3.0.228 annot-version=v3.0`.
The FASTA reader of the wgd code in the ksrates package (`read_fasta` in `wgd_ksrates/utils.py`) raises an error if a
header contains `:`, `;` or `?` anywhere, also in the description, and all 34727 headers of `Solly_iTAGv2.3` and 896 of
the 10660 headers of `Micpu_v3.0` did. Each header was therefore reduced to its first word (`>Solyc00g005000.2.1`,
`>19233`), which is also the name that the GFF is matched against, so the choice of the feature and attribute is not
affected. The sequence lines are unchanged.
- Result: `Solly_iTAGv2.3` 34727 sequences, `Micpu_v3.0` 10660 sequences; in both the sequence lines are
  byte-identical to the originals, all names are unique, and no header contains `:`, `;`, `?` or a space.
- Scripts (awk one-liners): `Solly_iTAGv2.3/rename_sequence_headers.sh` and `Micpu_v3.0/rename_sequence_headers.sh`.
  They refuse to overwrite an existing output, and refer to the original files under their names before they were
  moved; the commands below use the moved files.
- Logs: `ITAG2.3_cds.short_ids.fasta.mapping.tsv` and `MpusillaCCMP1545_228_v3.0.cds_primaryTranscriptOnly.short_ids.fa.mapping.tsv`
  (new name, then the full original header).
- The GFFs are not changed (`ITAG2.3_gene_models.gff3`, `MpusillaCCMP1545_228_v3.0.gene.gff3`). Pair for both:
  `mRNA` / `Name`.

| File | md5 |
|---|---|
| `Solly_iTAGv2.3/ITAG2.3_cds_ORIGINAL.fasta` | `a31e96e08d31de13671fbedf0fd5ad9a` |
| `ITAG2.3_cds.short_ids.fasta` (in `genome_species/`) | `ff46a9374754aa2780b8345665100373` |
| `Micpu_v3.0/MpusillaCCMP1545_228_v3.0.cds_primaryTranscriptOnly_ORIGINAL.fa` | `c40a72eb575fc96eb4a3b6080734de0d` |
| `MpusillaCCMP1545_228_v3.0.cds_primaryTranscriptOnly.short_ids.fa` (in `genome_species/`) | `bcb34ce17c6354f2bfc1bb8cd6e5ee6f` |

`cds_filepath` of both species in the manifest points to the `.short_ids` files.

### Reproduction

```bash
cd /path/to/1kp/source_data/2.genomes/1kp/genome_species
O=post_processing_fasta_gff_logs/Solly_iTAGv2.3
awk -v map="$O/ITAG2.3_cds.short_ids.fasta.mapping.tsv" '/^>/ { id = substr($1, 2); print ">" id; print id "\t" substr($0, 2) > map; next } { print }' \
    $O/ITAG2.3_cds_ORIGINAL.fasta > ITAG2.3_cds.short_ids.fasta
O=post_processing_fasta_gff_logs/Micpu_v3.0
awk -v map="$O/MpusillaCCMP1545_228_v3.0.cds_primaryTranscriptOnly.short_ids.fa.mapping.tsv" '/^>/ { id = substr($1, 2); print ">" id; print id "\t" substr($0, 2) > map; next } { print }' \
    $O/MpusillaCCMP1545_228_v3.0.cds_primaryTranscriptOnly_ORIGINAL.fa > MpusillaCCMP1545_228_v3.0.cds_primaryTranscriptOnly.short_ids.fa
```

## Pinta_v2.0 (*Pinus taeda*)

The manifest notes that the assembly on disk is v2.01 (v2.0 had no annotation files). Its annotation (`Pita.2_01.gtf`,
source column GFACS) is a GTF whose column 9 holds only the bare gene ID (e.g. `PITA_00001`), without `key=value`
attributes, so the GFF parser cannot use it. It was converted to GFF3 with an awk command: the header
`##gff-version 3` is added; in column 9, `gene` lines get `ID=PITA_00001;Name=PITA_00001` and all other lines
(`CDS`, `intron`, `start_codon`, `stop_codon`) get `Parent=PITA_00001`. Columns 1 to 8 are unchanged.
- Result: 519618 lines in, 519619 out (the header). 51751 `gene` lines, identical to the 51751 FASTA names
  (`Pita.2_01.cds.fa`, not changed), each name on exactly one `gene` line. Pair: `gene` / `ID` (equally valid:
  `gene` / `Name`).
- Limits: this is a minimal GFF3 that the parser can read, not a strictly valid one (there are no `mRNA` lines and the
  `Parent` of the other lines is the gene). The parser uses only the feature, the coordinates, the strand and the
  chosen attribute.
- Script: `Pinta_v2.0/convert_gft_to_gff3.sh` (the file name has a typo, "gft"). It refuses to overwrite an existing
  output, and reads `Pita.2_01.gtf` from `genome_species/`, the name before the original was moved; the command below
  uses the moved file.

| File | md5 |
|---|---|
| `Pinta_v2.0/Pita.2_01_ORIGINAL.gtf` | `f244b011891ebe7b071fbccc163b0c10` |
| `Pita.2_01.gff3` (in `genome_species/`) | `7940fa07610b4dd297d07ac634062320` |

`gff_filepath` of `Pinta_v2.0` in the manifest points to `Pita.2_01.gff3`.

### Reproduction

```bash
cd /path/to/1kp/source_data/2.genomes/1kp/genome_species
awk 'BEGIN { FS = OFS = "\t"; print "##gff-version 3" } /^#/ { print; next } NF == 9 { if ($3 == "gene") $9 = "ID=" $9 ";Name=" $9; else $9 = "Parent=" $9 } { print }' \
    post_processing_fasta_gff_logs/Pinta_v2.0/Pita.2_01_ORIGINAL.gtf > Pita.2_01.gff3
```

## Common to all species

`rename_fasta_to_locus_tag.sh` and `merge_gff_gene_fragments.sh` refuse to overwrite an existing output and never
change their input; each writes its logs next to its output. The scan is run for all genome species with
`python3 code/1kp_wgms/others/post_process_genome_fasta_gff/scan_gff_feature_attribute.py` (result: `source_data/2.genomes/1kp/gff_feature_attribute_scan.tsv`,
overwritten at every run).

Checks on 2026-10-13: every processed file listed above was regenerated from its original with the commands given,
into a scratch directory, and was identical to the file kept here, logs included. For every species with a GFF
(`Klefl_v1.0` has none), a separate count confirmed that the pair in the summary table finds each FASTA name on exactly
one GFF line, for all names (34719 of 34721 for `Chabr_v0.1`, see above).

Scripts version: not committed yet when this was written; commit hash of `code/1kp_wgms`: _to be filled in_.
