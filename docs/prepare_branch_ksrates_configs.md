The script is meant to be a **ksrates setup planner** for branch-centered subtrees.

Script:

```bash
code/1kp_wgms/run_ksrates_raw/prepare_branch_ksrates_configs.py
```

It takes subtree definitions, writes ksrates config files, optionally runs `ksrates init`, then collects the ortholog-pair lists that ksrates says are needed.

**1. Prepare Inputs**

Species metadata TSV:

```tsv
species	fasta_filename	latin_name	gff_filename
A1	/path/to/A1.fasta	A species 1	/path/to/A1.gff3
A2	/path/to/A2.fasta	A species 2	/path/to/A2.gff3
B1	/path/to/B1.fasta	B species 1	/path/to/B1.gff3
```

Required columns:

```text
species
fasta_filename
```

Optional:

```text
latin_name
gff_filename
gff_feature
gff_attribute
collinearity
```

`gff_feature` and `gff_attribute` are the GFF feature and attribute of the focal species that identify its genes (empty: the command-line defaults are used). `collinearity` (`yes` or `no`, `no` if the column is missing) is written into the config of that focal species and decides whether anchor pairs are requested from the Ks database; it can only be `yes` for a species with a GFF. `generate_species_metadata.py` fills all three columns.

Subtree TSV:

```tsv
branch_id	newick_tree	focal_species	target_species
branch_001	(((A1,A2),(B1,B2)),(C1,C2));	A1,A2,B1,B2	A1,A2,B1,B2,C1,C2
```

Required columns:

```text
branch_id
newick_tree
```

Optional:

```text
focal_species
target_species
```

If `focal_species` is omitted, the script uses **every species in the subtree as a focal species**.

**2. Basic Use**

Generate configs and command files only:

```bash
python3 code/1kp_wgms/run_ksrates_raw/prepare_branch_ksrates_configs.py \
  --datasets subtree_datasets.tsv \
  --species-metadata species_metadata.tsv \
  --out-dir ksrates_branch_setup
```

This writes configs but does not run ksrates yet.

**3. `ksrates init`**

The script does not run `ksrates init`. The Nextflow wrapper does (`run_init = true` in its `nextflow.config`, see `how_to_launch_batch_wrapper.md`): it starts one job for every config with `init_done = no` in `config_manifest.tsv` and then runs this script again, which collects the generated `ortholog_pairs_<focal>` files and deduplicates them.

**4. Main Outputs**

Inside `ksrates_branch_setup/`:

```text
configs/
  branch_001/config_branch_001_A1.txt

init_runs/
  branch_001/A1/rate_adjustment/A1/ortholog_pairs_A1.tsv

config_manifest.tsv
ortholog_pairs.tsv
ortholog_pairs_by_source.tsv
branches/branch_001/ortholog_pairs.tsv
```

Most important files:

- `ortholog_pairs.tsv`: global deduplicated pair list across all branches/focals, with a `queued` column (`no` for pairs whose peak and Ks list are already in the Ks database, which no ortholog job is started for; always `yes` when no `--ks-database-server-address-path` is given)
- `ortholog_pairs_by_source.tsv`: shows which focal config required each pair
- `paralog_runs/paralog_species.tsv`: one deduplicated species per row with the config chosen for it and a `queued` column; the Nextflow wrapper starts the Ks jobs from this table and from `ortholog_pairs.tsv`, there are no command files for them

**5. Pair Scope Modes**

Default mode:

```bash
--pair-scope init-pairs
```

This trusts ksrates’ `ortholog_pairs_<focal>` output directly. This is what I would start with.

Tighter mode:

```bash
--pair-scope target-trios
```

This reads `ortholog_trios_<focal>` and derives only pairs for the `target_species` listed in your subtree TSV. Use this if you want only the branch-flanking divergences rather than the full focal correction setup.

For your current idea, I’d start with:

```bash
--pair-scope init-pairs --max-outgroups 4
```

Then inspect `ortholog_pairs_by_source.tsv` to see whether ksrates is adding extra pairs beyond what we want.