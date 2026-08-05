# 1kp_wgms

WGM inference in One Thousand Plant (1KP) transcriptomes using rate-adjusted Ks distributions.

## Annotating 1KP WGMs on the species tree

`annotate_wgms_on_tree.py` reconstructs a practical WGM/WGD annotation for the rooted 1KP ASTRAL species tree. The 1KP/Barker data release provides the WGM summary table and the figure PDF, but not an author-supplied Newick file with all WGM labels already attached to branches. This script therefore infers branch placements from the supplementary evidence tables.

The current implementation uses ETE3 for tree parsing, MRCA calculation, and NHX tree writing. The supplementary `.xlsx` files are read directly with the Python standard library, so the only non-stdlib dependency is ETE3.

## Dependency

Install ETE3 in the Python environment where you will run the script:

```bash
pip install ete3
```

## Inputs

By default, the script uses these files in this project:

- `source_data/1.species_dataset/1kp_paper_2019_suptab2_wgm.xlsx`
  - Supplementary Table 2.
  - Defines the official WGM/WGD IDs to annotate.
- `source_data/1.species_dataset/1kp_paper_2019_suptab3_ks.xlsx`
  - Supplementary Table 3.
  - Lists the WGM/WGD history for each 1KP species code in columns `WGD 1`, `WGD 2`, and `WGD 3`.
- `source_data/3.phylogenetic_tree/1kp_trees/astral_trees_33_percent-FAA_estimated_species_tree.rooted.tree`
  - Rooted 1KP ASTRAL species tree.
  - Tip labels are expected to be 1KP sample/species codes matching Table 3.

## Placement Rule

For each WGM/WGD ID:

1. Read all 1KP species codes in Supplementary Table 3 that list that ID in `WGD 1`, `WGD 2`, or `WGD 3`.
2. Keep only IDs also present in Supplementary Table 2, unless `--include-unlisted` is used.
3. Find those taxa in the rooted species tree.
4. Use ETE3 to place the WGM/WGD label on the branch subtending their MRCA.

This is an inferred placement. It is useful for downstream tree annotation, but it is not the same thing as recovering a hidden author-provided annotation file. Some placements are necessarily approximate, especially when an event is supported by only one sampled taxon or when the supporting taxa are a subset of a larger MRCA clade.

## Outputs

The default outputs are written beside the input rooted tree:

- `astral_trees_33_percent-FAA_estimated_species_tree.rooted.wgm_suptab3_mrca.nhx.tree`
  - Newick tree with NHX features on annotated branches.
  - Example annotation for two events on one branch: `[&&NHX:WGM_label=AMBOalpha,ALINbeta:WGM_count=2:WGM1=AMBOalpha:WGM2=ALINbeta]`
- `astral_trees_33_percent-FAA_estimated_species_tree.rooted.wgm_suptab3_mrca.tsv`
  - Audit table with one row per placed WGM/WGD.
  - Includes support taxon counts, missing taxa, MRCA clade size, and a placement note.

Important TSV fields:

- `WGM`: WGM/WGD ID from Supplementary Table 2.
- `node_id`: internal script node ID for the placed MRCA branch.
- `support_taxa_count`: number of Table 3 taxa listing the event.
- `present_taxa_count`: support taxa found in the species tree.
- `missing_taxa_count`: support taxa absent from the species tree.
- `mrca_tip_count`: total number of tips under the placed MRCA.
- `exact_support_clade`: `true` if the support taxa exactly equal the MRCA clade.
- `placement_note`: short interpretation of the placement.

## Usage

Run with the project defaults:

```bash
python code/1kp_wgms/annotate_wgms_on_tree.py
```

Specify custom inputs or outputs:

```bash
python code/1kp_wgms/annotate_wgms_on_tree.py \
  --tree source_data/3.phylogenetic_tree/1kp_trees/astral_trees_33_percent-FAA_estimated_species_tree.rooted.tree \
  --wgm-table source_data/1.species_dataset/1kp_paper_2019_suptab2_wgm.xlsx \
  --ks-table source_data/1.species_dataset/1kp_paper_2019_suptab3_ks.xlsx \
  --out-tree source_data/3.phylogenetic_tree/1kp_trees/1kp_wgm_annotated.nhx.tree \
  --out-tsv source_data/3.phylogenetic_tree/1kp_trees/1kp_wgm_placements.tsv
```

Include event IDs found in Supplementary Table 3 even if they are absent from Supplementary Table 2:

```bash
python code/1kp_wgms/annotate_wgms_on_tree.py --include-unlisted
```

By default, the script reads and writes with ETE3 Newick `format=0`, which is intended for the ASTRAL-style tree with internal support values and branch lengths. If needed, change this with `--ete-format`.

## Interpreting the Annotated Tree

The WGM annotation is written as NHX features attached to the node/branch selected by the MRCA rule. Multiple WGMs on the same branch are written as separate keys (`WGM1`, `WGM2`, ...), plus `WGM_label` for a combined display string. Many tree tools preserve NHX comments, but some tools drop comments when re-saving trees. Keep the TSV as the authoritative audit trail, and prefer the iTOL text dataset for display.

Recommended checks after running:

1. Inspect rows where `present_taxa_count` is `1`; these are terminal-branch placements.
2. Inspect rows where `exact_support_clade` is `false`; these are MRCA placements where the support taxa do not cover the whole clade.
3. Compare high-interest events against `1KP_WGD_phylogeny.pdf` or Supplementary Figure 8 before using them as final curated branch placements.

## Visualization Notes

For quick inspection, the ETE3 PDF is useful. For final exploration or figures, the following approaches may work better for a 1000-tip tree:

1. Use iTOL with branch annotation datasets generated from the TSV. This gives zooming, clade collapsing, searchable labels, and cleaner export than a static all-in-one PDF.
2. Split the tree into major clades and render one PDF per clade. WGM labels become much easier to read when the tree is not forced onto one page.
3. Use a circular overview PDF without tip labels plus a separate TSV/table for exact placements. This is often clearer than trying to show every taxon name.
4. For publication-style figures, curate a simplified backbone tree and show WGM markers on that backbone, with the full annotated tree retained as supplementary data.



