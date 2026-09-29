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
```

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
  --subtrees branch_subtrees.tsv \
  --species-metadata species_metadata.tsv \
  --out-dir ksrates_branch_setup
```

This writes configs but does not run ksrates yet.

**3. Run `ksrates init` Automatically**

```bash
python3 code/1kp_wgms/run_ksrates_raw/prepare_branch_ksrates_configs.py \
  --subtrees branch_subtrees.tsv \
  --species-metadata species_metadata.tsv \
  --out-dir ksrates_branch_setup \
  --max-outgroups 4 \
  --run-init
```

After each `init`, the script collects the generated `ortholog_pairs_<focal>` files and deduplicates them.

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
init_commands.sh
orthologs_ks_commands.sh
branches/branch_001/ortholog_pairs.tsv
```

Most important files:

- `ortholog_pairs.tsv`: global deduplicated pair list across all branches/focals
- `ortholog_pairs_by_source.tsv`: shows which focal config required each pair
- `orthologs_ks_commands.sh`: commands to calculate the needed ortholog Ks distributions

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