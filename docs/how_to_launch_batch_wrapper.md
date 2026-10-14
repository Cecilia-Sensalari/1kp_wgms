# How to launch the batch_ksrates wrapper

## 1. Make sure the paralog Ks database server is running

The server is independent of the wrapper - launch it once, separately, and leave it running.

```bash
cd /group/esb/cesen/1kp/ks_analysis/1kp
sbatch /group/esb/cesen/1kp/code/1kp_wgms/launchers/launch_ks_database_server.sh
```

It reuses the existing `ks_database/ks_database_sqld_data`/`ks_database_sqld_keys` if
already present (no data lost across restarts), and writes a fresh `host:port` to
`ks_database/ks_database_server_address.txt`. Check that file's content to see which
node it landed on, and confirm it's actually up before running the wrapper:

```bash
cat /group/esb/cesen/1kp/ks_analysis/1kp/ks_database/ks_database_server_address.txt
```

If this job ever gets cancelled or hits a walltime limit, just resubmit the same command - the
data/keys persist, only the address file goes stale until resubmitted.

## 2. Configure and run the wrapper

All settings of the wrapper live in Nextflow config files, not in `main.nf`:

- `batch_ksrates/nextflow.config` (in the repo) holds the defaults of every parameter and of every
  process setting (executor, memory, CPUs, time, number of jobs running at once, error handling,
  `module load` lines). Do not edit it for one analysis.
- Each analysis has its own `nextflow.config` in the directory it is launched from (for example
  `ks_analysis/1kp_quicktest/nextflow.config`), containing only what differs: at least `base_dir`,
  `ks_database_server_address_path` and `expert_config`; and, for a bigger run, higher `maxForks`,
  memory, and so on. Nextflow loads it automatically and it overrides the defaults. Order of
  priority, highest first: command-line `--params`, `-c` files, `nextflow.config` in the launch
  directory, the repo's `nextflow.config`.
- The resources of the per-dataset ksrates pipelines (stage 4) are in `batch_ksrates/ksrates_run.config`.
  To change them for an analysis, copy that file, edit it and set `params.ksrates_run_config` in the
  analysis' `nextflow.config`.

Example `nextflow.config` of an analysis:

```groovy
params {
    base_dir = '/group/esb/cesen/1kp/ks_analysis/1kp'
    ks_database_server_address_path = '/group/esb/cesen/1kp/ks_analysis/1kp/ks_database/ks_database_server_address.txt'
    expert_config = '/group/esb/cesen/1kp/ks_analysis/1kp/config_expert.txt'
}
process {
    withName: 'runOrthologJob' { maxForks = 200 }
    withName: 'runParalogJob'  { maxForks = 200 }
}
```

Then launch it with the launcher, which moves to the analysis directory (so its `nextflow.config` is used) and runs Nextflow with `-resume`; it repeats none of the settings. Which stages run, like every other setting, is chosen in the analysis' `nextflow.config` and not on the command line: the config file is the permanent record of how the analysis was run, and it is the same for every rerun.

```bash
sbatch /group/esb/cesen/1kp/code/1kp_wgms/launchers/launch_batch_ksrates.sh
```

A command-line option such as `--run_init false` overrides the config for one launch (Nextflow gives it priority). Use it only for a quick test, since nothing keeps a record of it.

The Ks database server launcher reads the container image, ksrates code directory and database location from the same `nextflow.config`.

Parameters to set in the analysis' `nextflow.config` (the full list, with defaults and descriptions, is in
`batch_ksrates/nextflow.config`):

- `paranome`, `collinearity`, `reciprocal_retention` (`"yes"` or `"no"`): which analysis types
  to request in every generated config. Collinearity is only used for the focal species that have a GFF
  (the genome species); the others always get `no`.
- `run_init = true`: runs `ksrates init` as one Nextflow job per config whose init output is missing (`init_done = no`
  in `config_manifest.tsv`; configs already initialised are skipped), and when they have all ended runs the config
  step again, which collects the ortholog pairs from the init output. The ortholog jobs need these pairs, so
  leave it on until init has been run for all configs. The paralog jobs do not need init and start at once.
- `execute_paralog_jobs`, `execute_ortholog_jobs`, `execute_real_pipelines` (`true` or `false`): which
  stages to actually run as jobs. They can be combined in one launch (also with `run_init = true`): the per-dataset pipelines then wait until every ortholog and paralog job of that launch has ended (also the failed ones, which are ignored), so they find the Ks values in the database. Check the Ks jobs before relying on the pipelines of a combined run.
- `expert_config`: needed if you want non-default ksrates expert settings; `use_ks_database_server = yes`
  must be set there, or the wrapper stops.

## 3. What the wrapper actually does

Three steps, the first two in parallel:

1. Regenerates `species_metadata.tsv` from scratch (CNGB transcriptomes + the genome manifest
   for the 31 species without transcriptomes). For each genome species it also takes the GFF feature and
   attribute that identify its genes (`genome_gff_settings.tsv`, from the checked `gff_feature_attribute_scan.tsv`;
   Chlva uses CDS/ID by hand; Klefl has no GFF) into the columns `gff_feature` and `gff_attribute`. The config
   of each focal species gets its own pair, and `collinearity = no` if the focal species has no GFF.
2. Regenerates the branch table (`focal_branch_to_dataset.tsv`, one row per branch) and the dataset list (`subtree_datasets.tsv`, one row per distinct dataset: branches that select the same species set are merged) from scratch, from the backbone tree. The configs are made from `subtree_datasets.tsv`.
3. Generates one ksrates config per (dataset, focal species) under `datasets/<hash>/`, plus
   `config_manifest.tsv`, `ortholog_runs/ortholog_pairs.tsv` and (since
   `ks_database_server_address_path` is set) `paralog_runs/paralog_species.tsv`. The two last tables have a
   column `queued` (`yes`/`no`) that marks the pairs and species whose Ks data is not yet in the Ks database
   (checked with `ksrates check-ortholog-db` and `check-paralog-db`), so rerunning after adding datasets only
   queues what is new. If the database cannot be reached the wrapper stops with an error.

While the `execute_...` parameters and `run_init` are `false`, nothing expensive *runs*: this only generates configs and tables.
With `execute_ortholog_jobs` and `execute_paralog_jobs` set to `true`, one Nextflow job is started for every row marked
`yes` in `ortholog_pairs.tsv` and `paralog_species.tsv`. The ksrates command is written once, in the process
(`runOrthologJob`, `runParalogJob` in `main.nf`), from the species, the pair and the config of the row. Each job
writes its ksrates log to `ortholog_runs/ortholog_ksrates_logs/<A>_<B>.log` or
`paralog_runs/paralog_ksrates_logs/<species>.log`; the log of an earlier run is renamed with its modification
time, so the plain name is always the latest run.

`ortholog_pairs.tsv` has no rows until `ksrates init` has been run for at least some configs (`run_init = true`).
The log of each init job is `datasets/<hash>/<focal>/init.log`.

## 4. What to check afterward

- The printed summary from `generateSpeciesMetadata` should say something like
  `(X from transcriptomes, 31 from genomes)`.
- `config_manifest.tsv` should include datasets referencing the 31 genome species without
  errors.
- `paralog_runs/paralog_species.tsv` should exist now, with a `queued` column showing which
  species still need `paralogs-ks` run versus which are already fully populated in the
  database.
