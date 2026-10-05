#!/usr/bin/env nextflow

/*
 * batch_ksrates: generates species_metadata.tsv, the dataset list (branch_subtrees.tsv), and
 * per-dataset ksrates configs/command files for the centralized paralog/ortholog Ks database
 * workflow (wraps generate_species_metadata.py, make_branch_subtrees.py, and
 * prepare_branch_ksrates_configs.py from ../run_ksrates_raw).
 *
 * Does NOT launch/stop the paralog Ks sqld server - that stays its own independent,
 * long-running process; point --paralog_database at its already-running address file.
 * Does NOT execute any of the generated commands (ksrates init/orthologs-ks/paralogs-ks stay a
 * separate, later, separately-parallelized step) unless --run_init is explicitly turned on.
 */

LOG_OUTPUT = true

/*
 * Pipeline input parameters
 */

params.base_dir = null  // required: output base directory, e.g. /group/esb/cesen/1kp/ks_analysis/1kp

// generate_species_metadata.py passthrough (defaults match that script's own hardcoded constants)
params.species_table = "/group/esb/cesen/1kp/source_data/1.species_dataset/1kp_paper_2019_suptab1_species.tsv"
params.transcriptomes_dir = "/group/esb/cesen/1kp/source_data/2.transcriptomes/1kp/cngb_transcriptomes"

// make_branch_subtrees.py passthrough (defaults match that script's own argparse defaults)
params.tree = "/group/esb/cesen/1kp/source_data/3.phylogenetic_tree/1kp_trees/astral_trees_33_percent-FAA_estimated_species_tree.rooted.wgm_suptab3_mrca.nhx.manual_fixes.tree"
params.score_table = "/group/esb/cesen/1kp/source_data/1.species_dataset/1kp_paper_2019_suptab1_species.tsv"
params.max_a = 4
params.max_older = 2
params.older_clades = 3
params.min_species = 3
params.busco_weight = 0.8
params.transrate_weight = 0.2
// If set, keep only the first N rows of the generated dataset list - for quick end-to-end
// testing (e.g. with --run_init and --paralog_database also on) without paying for the full
// ~2354-dataset tree. Unset (null) by default: full real run, unaffected.
params.test_max_datasets = null

// prepare_branch_ksrates_configs.py passthrough (defaults match that script's own argparse defaults)
params.expert_config = null
params.paralog_database = null
params.paralog_threads = 4
params.paralog_check_workers = 16
params.ortholog_threads = 4
// Plain production form by default; override with the "--env PYTHONPATH=/path/to/ksrates"
// variant (see run_in_container.sh in the ksrates checkout) when testing uncommitted ksrates
// code changes against the container's baked-in (stale) package instead.
params.ksrates_command = "singularity exec -B /group/esb/cesen/1kp /group/esb/cesen/1kp/software/ksrates_paralog_ks_db.sif ksrates"
params.pair_scope = "init-pairs"
params.max_outgroups = 4
params.consensus_mode = "mean among outgroups"
params.paranome = "yes"
params.collinearity = "no"
params.reciprocal_retention = "no"
params.gff_feature = "mrna"
params.gff_attribute = "id"
params.num_bootstrap_iterations = 200
params.max_ks_paralogs = "5"
params.max_ks_orthologs = "10"
params.keep_relative_paths = false
// Actually execute "ksrates init" for every generated config inline. Off by default: running
// init for all ~2354 real datasets in one step would be far too slow (see
// launchers/launch_ortholog_ks.sh's own comment) - init_commands.sh is still generated either
// way, for running separately, split across many parallel jobs.
params.run_init = false

// Actually submit orthologs_ks_commands.sh / paralogs_ks_commands.sh as native Nextflow array
// jobs on the slurm executor, instead of leaving them for a human to run manually. Off by
// default; independently controllable since you might want to launch one without the other.
// ksrates init (run_init/init_commands.sh above) is explicitly NOT covered by these flags - it
// keeps its current inline, sequential behavior for now.
params.execute_ortholog_jobs = false
params.execute_paralog_jobs = false
params.ortholog_job_memory = '4GB'
params.paralog_job_memory = '8GB'  // paralogs-ks is the heavier of the two
params.ortholog_job_concurrency = 25
params.paralog_job_concurrency = 25

// Stage 4: actually launch the "real" per-(dataset, focal) ksrates pipeline (nextflow run
// VIB-PSB/ksrates) for every row in config_manifest.tsv, now that the shared paralog/ortholog
// Ks data has already been populated by the jobs above. Deliberately NOT auto-gated on
// execute_ortholog_jobs/execute_paralog_jobs finishing in the same invocation - this is its own,
// separately-triggered step, run once you've checked Stage 2/3's results yourself.
params.execute_real_pipelines = false
// Local checkout's main.nf, not "VIB-PSB/ksrates" pulled from GitHub - this has the DB-aware
// setParalogAnalysis/doRateAdjustment skip gates (see docs/paralog_ks_database.rst) on its
// paralog_ks_db branch, which isn't published/merged upstream yet. Using the local path also
// means no network pull is needed and whatever's actually checked out here is what runs,
// including any local-only commits ahead of origin. The default VIB-PSB/ksrates branch has no
// DB-aware gate and would schedule wgdParalogs/wgdOrthologs for every focal species regardless
// of what the shared database already has (the Python package's own internal skip check was
// removed in favor of this Nextflow-level one - see the "DONE" plan section on this). Revisit
// once paralog_ks_db is actually merged/published - then this can go back to
// "VIB-PSB/ksrates -r <tag>".
params.ksrates_main_nf = '/group/esb/cesen/1kp/code/ksrates/main.nf'
params.ksrates_container = '/group/esb/cesen/1kp/software/ksrates_paralog_ks_db.sif'
params.real_pipeline_concurrency = 25
params.real_pipeline_time = '04:00:00'

if (!params.base_dir) {
    error "base_dir is required, e.g. --base_dir /group/esb/cesen/1kp/ks_analysis/1kp"
}
base_dir = file(params.base_dir)

if (params.execute_real_pipelines && !params.expert_config) {
    error "execute_real_pipelines requires --expert_config (it needs use_paralog_ks_database " +
        "= yes to actually benefit from the already-populated database, rather than silently " +
        "recomputing everything from scratch)"
}

// Same Python environment used by the existing launcher scripts (e.g.
// 1kp_wgms/launchers/launch_prepare_configs_test.sh) - these are plain host-side Python
// scripts, not run inside the ksrates container (only the "ksrates" subcommands prepareConfigs
// shells out to, via --ksrates-command, run inside the container).
PYTHON_SETUP = "module load python/x86_64/3.11.4 && source /group/esb/cesen/.venv/venv_python_11/bin/activate"

SCRIPT_DIR = "${workflow.projectDir}/../run_ksrates_raw"

// Generated once at startup (not per-task - nothing here depends on any process output), only
// when actually needed. Mirrors ks_analysis/test_subtree/nextflow.config's shape (singularity
// profile + per-process resource withName blocks) but points at our own local container/bind
// mounts instead of the public vibpsb/ksrates:latest image, and is rewritten fresh on every run
// so a changed --paralog_database/--transcriptomes_dir is always reflected.
if (params.execute_real_pipelines) {
    bind_mounts = [base_dir.toString(), params.transcriptomes_dir]
    if (params.paralog_database) {
        bind_mounts << file(params.paralog_database).getParent().toString()
    }
    bind_options = bind_mounts.unique().collect { "-B ${it}" }.join(' ')

    file("${base_dir}/ksrates_run.nextflow.config").text = """
profiles {
    singularity {
        singularity.enabled = true
        singularity.cacheDir = '${params.base_dir}'
        singularity.autoMounts = true
        singularity.runOptions = "${bind_options}"
    }
}

executor {
    name = 'slurm'
}

process {
    container = '${params.ksrates_container}'

    withName: 'wgdParalogs' {
        cpus = ${params.paralog_threads}
        memory = '${params.paralog_job_memory}'
    }
    withName: 'wgdOrthologs' {
        cpus = ${params.ortholog_threads}
        memory = '${params.ortholog_job_memory}'
    }
    withName: 'estimatePeaks' {
        memory = '4GB'
    }
    withName: 'paralogsAnalyses' {
        memory = '8GB'
    }
    withName: 'plotOrthologDistrib' {
        memory = '2GB'
    }
}
""".stripIndent()
}

log.info """
         batch_ksrates
         =============
         base_dir:          ${params.base_dir}
         tree:               ${params.tree}
         species_table:      ${params.species_table}
         expert_config:      ${params.expert_config ?: '(none)'}
         paralog_database:   ${params.paralog_database ?: '(none)'}
         run_init:           ${params.run_init}
         test_max_datasets:  ${params.test_max_datasets ?: '(none - full run)'}
         execute_ortholog_jobs: ${params.execute_ortholog_jobs}
         execute_paralog_jobs:  ${params.execute_paralog_jobs}
         execute_real_pipelines: ${params.execute_real_pipelines}
         ksrates_main_nf:    ${params.ksrates_main_nf}
         """
         .stripIndent()


workflow {
    generateSpeciesMetadata()
    generateDatasetList()
    prepareConfigs(generateSpeciesMetadata.out.species_metadata, generateDatasetList.out.dataset_list)
    prepareConfigs.out.summary.view()

    // Gated behind prepareConfigs.out.summary (not a bare Channel.fromPath on the command file)
    // so Nextflow waits for prepareConfigs to actually finish writing it first - these files
    // aren't Nextflow outputs themselves, just plain paths on disk.
    if (params.execute_ortholog_jobs) {
        ortholog_commands = prepareConfigs.out.summary
            .flatMap { file("${base_dir}/orthologs_ks_commands.sh").readLines() }
            .filter { line -> line.trim() && !line.startsWith('#') && line.trim() != 'set -euo pipefail' }
        runOrthologJob(ortholog_commands)
    }
    if (params.execute_paralog_jobs) {
        paralog_commands = prepareConfigs.out.summary
            .flatMap { file("${base_dir}/paralogs_ks_commands.sh").readLines() }
            .filter { line -> line.trim() && !line.startsWith('#') && line.trim() != 'set -euo pipefail' }
        runParalogJob(paralog_commands)
    }
    if (params.execute_real_pipelines) {
        real_pipeline_rows = prepareConfigs.out.summary
            .flatMap { file("${base_dir}/config_manifest.tsv").readLines().drop(1) }
            .map { line -> line.split('\t') }
            .map { fields -> tuple(fields[0], fields[1], fields[2], fields[3]) }
        runRealPipeline(real_pipeline_rows)
    }
}


/*
 * Regenerates species_metadata.tsv from scratch every run, by joining the fixed 1kp species
 * table against whatever CNGB transcriptome FASTA files are actually present on disk - safe to
 * overwrite unconditionally (see generate_species_metadata.py's own docstring).
 */
process generateSpeciesMetadata {

    executor 'local'
    beforeScript PYTHON_SETUP

    output:
        env SPECIES_METADATA, emit: species_metadata

    script:
    """
    python3 ${SCRIPT_DIR}/generate_species_metadata.py \
        --species-table ${params.species_table} \
        --transcriptomes-dir ${params.transcriptomes_dir} \
        --output ${base_dir}/species_metadata.tsv
    SPECIES_METADATA=${base_dir}/species_metadata.tsv
    """
}


/*
 * Regenerates the full dataset list from scratch every run, derived entirely from the one fixed
 * backbone tree file - safe to overwrite unconditionally (see make_branch_subtrees.py's own
 * docstring). Runs independently of generateSpeciesMetadata; Nextflow parallelizes them.
 */
process generateDatasetList {

    executor 'local'
    beforeScript PYTHON_SETUP

    output:
        env DATASET_LIST, emit: dataset_list

    script:
    // Groovy-level conditional, not a bash-level one: interpolating an unset params.* value
    // directly into the script string would literally paste the text "null" into the bash
    // script below, not behave as empty/false. Passed straight to make_branch_subtrees.py so
    // it stops traversing the tree early instead of walking/pruning all ~2354 branches and
    // throwing most of the output away afterwards.
    max_datasets_arg = params.test_max_datasets ? "--max-datasets ${params.test_max_datasets}" : ""
    """
    python3 ${SCRIPT_DIR}/make_branch_subtrees.py \
        --tree ${params.tree} \
        --out ${base_dir}/branch_subtrees.tsv \
        --score-table ${params.score_table} \
        --max-a ${params.max_a} \
        --max-older ${params.max_older} \
        --older-clades ${params.older_clades} \
        --min-species ${params.min_species} \
        --busco-weight ${params.busco_weight} \
        --transrate-weight ${params.transrate_weight} \
        ${max_datasets_arg}
    DATASET_LIST=${base_dir}/branch_subtrees.tsv
    """
}


/*
 * Generates one ksrates config per (dataset, focal species) under base_dir/datasets/<hash>/,
 * plus the deduplicated command files (init_commands.sh, orthologs_ks_commands.sh, and
 * paralogs_ks_commands.sh if --paralog_database is set). Depends on both processes above;
 * Nextflow waits for both automatically via the input channels below.
 */
process prepareConfigs {

    executor 'local'
    beforeScript PYTHON_SETUP

    input:
        val species_metadata
        val dataset_list

    output:
        stdout emit: summary

    script:
    expert_arg = params.expert_config ? "--expert-config ${params.expert_config}" : ""
    paralog_db_arg = params.paralog_database ? "--paralog-database ${params.paralog_database}" : ""
    run_init_arg = params.run_init ? "--run-init" : ""
    keep_relative_arg = params.keep_relative_paths ? "--keep-relative-paths" : ""
    """
    python3 ${SCRIPT_DIR}/prepare_branch_ksrates_configs.py \
        --datasets ${dataset_list} \
        --species-metadata ${species_metadata} \
        --out-dir ${base_dir} \
        ${expert_arg} \
        ${paralog_db_arg} \
        --paralog-threads ${params.paralog_threads} \
        --paralog-check-workers ${params.paralog_check_workers} \
        --ortholog-threads ${params.ortholog_threads} \
        --ksrates-command "${params.ksrates_command}" \
        --pair-scope ${params.pair_scope} \
        --max-outgroups ${params.max_outgroups} \
        --consensus-mode "${params.consensus_mode}" \
        --paranome ${params.paranome} \
        --collinearity ${params.collinearity} \
        --reciprocal-retention ${params.reciprocal_retention} \
        --gff-feature ${params.gff_feature} \
        --gff-attribute ${params.gff_attribute} \
        --num-bootstrap-iterations ${params.num_bootstrap_iterations} \
        --max-ks-paralogs ${params.max_ks_paralogs} \
        --max-ks-orthologs ${params.max_ks_orthologs} \
        ${run_init_arg} \
        ${keep_relative_arg}
    """
}


/*
 * Runs one line of orthologs_ks_commands.sh (already a complete, self-contained
 * "( cd <cwd> && orthologs-ks ... && orthologs-analysis ... )" subshell from write_command_file())
 * as its own slurm job. One process instance per unique ortholog pair.
 */
process runOrthologJob {

    executor 'slurm'
    // A failing pair (e.g. too divergent for max_ks_orthologs, see batch_ksrates docs/notes)
    // shouldn't abort the whole batch - this job family is expected to have some unavoidable
    // failures at ~2354-dataset scale. Rerun with -resume later to retry just the failed ones
    // once addressed (or accept them as permanently uncomputable under current settings).
    errorStrategy 'ignore'
    maxForks params.ortholog_job_concurrency
    cpus params.ortholog_threads
    memory params.ortholog_job_memory

    input:
        val command_line

    script:
    """
    ${command_line}
    """
}


/*
 * Runs one line of paralogs_ks_commands.sh (already a complete, self-contained
 * "( cd <cwd> && paralogs-ks ... )" subshell from write_command_file()) as its own slurm job.
 * One process instance per unique paralog species.
 */
process runParalogJob {

    executor 'slurm'
    // Same reasoning as runOrthologJob's errorStrategy above: don't let one bad species abort
    // the whole batch.
    errorStrategy 'ignore'
    maxForks params.paralog_job_concurrency
    cpus params.paralog_threads
    memory params.paralog_job_memory

    input:
        val command_line

    script:
    """
    ${command_line}
    """
}


/*
 * Runs the "real" ksrates pipeline for one (dataset, focal) combination from
 * config_manifest.tsv, from that combination's own datasets/<id>/<focal>/ directory (the same
 * cwd "ksrates init" already used). Runs params.ksrates_main_nf directly (a local checkout path,
 * see its own comment above) instead of "nextflow run VIB-PSB/ksrates -r <branch>" - avoids
 * depending on the branch being published/reachable on GitHub at all. Reads whatever's already
 * in config_manifest.tsv and the shared databases at the time this runs - same as a human
 * manually running these commands one by one would do.
 */
process runRealPipeline {

    tag "${dataset_id}/${focal_species}"
    executor 'slurm'
    errorStrategy 'ignore'
    maxForks params.real_pipeline_concurrency
    time params.real_pipeline_time
    memory '2GB'
    beforeScript "module load nextflow"

    input:
        tuple val(dataset_id), val(focal_species), val(config_path), val(init_run_dir)

    script:
    """
    cd ${init_run_dir}
    nextflow run ${params.ksrates_main_nf} -profile singularity \
        --config ${config_path} \
        --expert ${params.expert_config} \
        -c ${base_dir}/ksrates_run.nextflow.config
    """
}
