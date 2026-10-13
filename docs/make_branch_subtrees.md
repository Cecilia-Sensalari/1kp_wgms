# `make_branch_subtrees.py`

This script creates the tables used to set up branch-centered `ksrates` analyses:
`focal_branch_to_dataset.tsv` (one row per branch) and `subtree_datasets.tsv` (one row per distinct dataset).
See `dataset_generation_summary.md` for a one-page description of how the datasets are built. It reads the rooted 1KP species tree, visits every non-root
branch, and writes a small local subtree around that branch.

## How the Local Subtree Is Defined

For each branch, the node below the branch is called `A`.
The script then walks from `A` toward the root and records progressively older
sister clades:

```text
(((A, B), C), D)
```

The retained species are:

- `A_species`: up to four descendants of the test branch.
- `B_species`: up to two species from the sister clade of `A`.
- `C_species`: up to two species from the next older sister clade.
- `D_species`: up to two species from the next older sister clade, if present.

For `A`, the rule is based on the two child subclades below the test node:

- leaf branch: keep the single leaf;
- two descendant leaves: keep both;
- larger clades: keep up to two species from each child subclade.

When `A` is an internal node, its two child subclades can be thought of as
`A'` and `A''`. If `A_species` contains three or four leaves, those leaves are
not just the best three or four species anywhere below `A`; they are sampled
from both `A'` and `A''` so that the descendant side of the test branch remains
phylogenetically balanced.

When there are more candidate leaves than slots, representatives are chosen by
phylogenetic balance first, then by transcriptome quality.

## Quality Score

The script reads the 1KP species table:

```text
/path/to/1kp/source_data/1.species_dataset/1kp_paper_2019_suptab1_species.tsv
```

and scores transcriptomes as:

```text
quality_score = 0.8 * BUSCO_percentile + 0.2 * TransRate_percentile
```

Higher scores are preferred when choosing among leaves within the same selected
subclade.

## Genome BUSCO scores (methods)

The 31 genome-based species are not in the 1KP score table, and TransRate cannot be computed for them (it needs sequencing reads). To give them a quality score comparable with the transcriptomes, BUSCO completeness was measured for them as close as possible to the 1KP protocol.

**What 1KP did** (Leebens-Mack et al. 2019, *Nature* 574:679-685, Methods and Extended Data Fig. 3; Carpenter et al. 2019, *GigaScience* 8:giz126): BUSCO with default settings, using the eukaryote (`eukaryota_odb9`) and embryophyte (`embryophyta_odb9`) gene sets, on the proteins predicted from each assembly. BUSCO completeness is reported as "complete + fragmented" (= 100% - missing), and samples with more than 42.5% missing eukaryota BUSCOs were removed. The BUSCO version is not stated unambiguously in the sources we could read. The BUSCO column used here (`% BUSCOs (complete+fragmented)`) is taken to be the `eukaryota_odb9` set, because its values (30-98%, mean 87%) are in the range of that 303-gene set and not of the 1440-gene embryophyta set; the paper does not say so explicitly.

**What was done for the genomes**

1. Input: the CDS FASTA of the primary transcripts of each genome (`genome_species_manifest.tsv`, column `cds_filepath`), after the FASTA/GFF post-processing described in the genome log.
2. Translation of each CDS to protein with the standard genetic code (trailing stop removed, internal stops and ambiguous codons written as `X`), one protein file per species. Script: `others/translate_cds_to_protein.py` (output in `genome_busco/translated_fasta/`).
3. BUSCO 3.0.2b (cluster module `busco/x86_64/3.0.2b`, protein mode `-m proteins`, default settings, `-c 4`) against `eukaryota_odb9` (303 BUSCOs). BUSCO 3 is the latest version that reads the odb9 gene sets (BUSCO 4 and later require odb10), which are needed to match 1KP. The odb9 data are no longer distributed by the BUSCO authors and were taken from the archived copy at https://zenodo.org/records/20987183. Script: `others/run_busco_genomes.sh` (one SLURM array task per species).
4. Collection: `others/collect_busco_scores.py` reads the BUSCO short summaries and writes `genome_busco_scores.tsv` with `busco_complete_fragmented = 100 - missing`, plus the complete, fragmented and missing percentages and the number of BUSCOs (303 for all species).

**How the score is used.** The genome BUSCO values are ranked together with the 1KP transcriptome values, giving each genome a BUSCO percentile. The quality score of a genome is this percentile alone, since it has no TransRate (for transcriptomes: 0.8 x BUSCO percentile + 0.2 x TransRate percentile). TransRate measures read-mapping errors of an assembly and is not defined for an annotated genome, so it cannot be used. This is a convention, not a claim that a genome would score like a transcriptome on TransRate: the BUSCO weight is rescaled from 0.8 to 1 so that the only available measure carries the whole score. A genome and a transcriptome are therefore compared on BUSCO alone for the genome and on BUSCO plus TransRate for the transcriptome, and a transcriptome with a high BUSCO but a low TransRate percentile scores lower than a genome with the same BUSCO percentile.

**Limitations to state when publishing**

- BUSCO 3.0.2b is not necessarily the version 1KP used, and BUSCO was applied to annotated gene models of genomes, not to proteins predicted from assemblies; the values are comparable in kind, not identical in method.
- Most genomes score above the highest transcriptome (98.3%), so they are preferred as representatives whenever a clade offers one. Pinta (62.7%) and Chlva (75.2%) rank low.
- The score only decides which species represent a clade in a dataset; it does not enter the Ks analysis itself.

## Output Columns: `focal_branch_to_dataset.tsv`

One row per branch.

- `dataset_id`: first 8 characters of the SHA-1 of the sorted selected species. It depends only on the species set, so different branches that select the same species get the same id.
- `traversal_order`: position of the branch in the preorder walk of the tree (unique per branch).
- `newick_tree`: pruned Newick tree containing only the selected species.
- `focal_species`: species that should be treated as focal species by downstream `ksrates` setup.
  This is `A_species + B_species + C_species`.
- `target_species`: all selected species in the subtree, including `D_species`.
- `A_species`, `B_species`, `C_species`, `D_species`: selected representatives by local clade.
- `n_selected_species`: total number of species retained in the row.

## Output Columns: `subtree_datasets.tsv`

One row per distinct `dataset_id`; this is the input of `prepare_branch_ksrates_configs.py`. Branches with the same `dataset_id` have the same tree and target species, so their ksrates runs are identical and they are generated once. For example the two sister branches below a node can select the same species set, with A and B swapped.

- `dataset_id`, `newick_tree`, `target_species`, `n_selected_species`: as above.
- `focal_species`: union of the focal species of the branches that share the dataset.
- `branches`: the `traversal_order` values of those branches.

To interpret a branch, look up its `dataset_id` in `focal_branch_to_dataset.tsv`, then use the output of the focal species of its `A_species` under `datasets/<dataset_id>/`. The script stops with an error if two branches with the same `dataset_id` have different trees or species.

The largest expected subtree has ten species: four from `A` and two each from
`B`, `C`, and `D`.

## Focal Species Logic

Only `A`, `B`, and `C` representatives are focal species by default:

- `A_species`: expected to share a putative WGM on the tested branch.
- `B_species`: closest species expected not to share that WGM.
- `C_species`: older comparison that gives more depth to the Ks plot, to see how the signal relates to the putative WGM on the tested branch.

`D_species` is kept in the subtree as older context for correction/outgroup
choice, but is not treated as focal by default.

## Run

```bash
python3 code/1kp_wgms/run_ksrates_raw/make_branch_subtrees.py \
  --out /path/to/ks_analysis/1kp/ksrates_raw/focal_branch_to_dataset.tsv \
  --datasets-out /path/to/ks_analysis/1kp/ksrates_raw/subtree_datasets.tsv
```

The script requires `ete3` in the active Python environment.
