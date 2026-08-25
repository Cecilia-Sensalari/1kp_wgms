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

## Output Columns

- `branch_id`: traversal number plus either the leaf 1KP ID or `internal`, e.g. `000651_JKAA` or `000650_internal`.
- `newick_tree`: pruned Newick tree containing only the selected species.
- `focal_species`: species that should be treated as focal species by downstream `ksrates` setup.
  This is `A_species + B_species + C_species`.
- `target_species`: all selected species in the subtree, including `D_species`.
- `A_species`, `B_species`, `C_species`, `D_species`: selected representatives by local clade.
- `n_selected_species`: total number of species retained in the row.

The largest expected subtree has ten species: four from `A` and two each from
`B`, `C`, and `D`.

## Focal Species Logic

Only `A`, `B`, and `C` representatives are focal species by default:

- `A_species`: expected to share a putative WGM on the tested branch.
- `B_species`: closest species expected not to share that WGM.
- `C_species`: older negative comparison, used to confirm the signal outside `A`.

`D_species` is kept in the subtree as older context for correction/outgroup
choice, but is not treated as focal by default.

## Run

```bash
python3 code/1kp_wgms/run_ksrates_raw/make_branch_subtrees.py \
  --out /path/to/ks_analysis/1kp/ksrates_raw/branch_subtrees.tsv
```

The script requires `ete3` in the active Python environment.
