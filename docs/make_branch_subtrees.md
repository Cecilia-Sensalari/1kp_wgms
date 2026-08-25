# `make_branch_subtrees.py`

This script creates the `branch_subtrees.tsv` input used to set up branch-centered
`ksrates` analyses. It reads the rooted 1KP species tree, visits every non-root
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
- three descendant leaves: keep all three, split as `1 + 2` across the two child subclades;
- larger clades: keep up to two species from each child subclade.

When there are more candidate leaves than slots, representatives are chosen by
phylogenetic balance first, then by transcriptome quality.

## Quality Score

The script reads the 1KP species table:

```text
/group/esb/cesen/1kp/source_data/1.species_dataset/1kp_paper_2019_suptab1_species.tsv
```

and scores transcriptomes as:

```text
quality_score = 0.8 * BUSCO_percentile + 0.2 * TransRate_percentile
```

Higher scores are preferred when choosing among leaves within the same selected
subclade.

## Output Columns

- `branch_id`: traversal number plus either the leaf 1KP ID or `internal`, e.g. `000651_JKAA` or `000650_internal`.
- `newick_tree`: pruned Newick tree containing only the selected species.
- `focal_species`: currently all selected species, for downstream ksrates config generation.
- `target_species`: currently all selected species.
- `A_species`, `B_species`, `C_species`, `D_species`: selected representatives by local clade.
- `n_selected_species`: total number of species retained in the row.

The largest expected subtree has ten species: four from `A` and two each from
`B`, `C`, and `D`.

## Run

```bash
python3 code/1kp_wgms/run_ksrates_raw/make_branch_subtrees.py \
  --out /group/esb/cesen/1kp/ks_analysis/1kp/ksrates_raw/branch_subtrees.tsv
```

The script requires `ete3` in the active Python environment.
