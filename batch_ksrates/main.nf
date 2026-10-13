#!/usr/bin/env nextflow

/*
 * batch_ksrates: generates species_metadata.tsv, the branch table (focal_branch_to_dataset.tsv), the dataset list (subtree_datasets.tsv), and
 * per-dataset ksrates configs/command files for the centralized paralog/ortholog Ks database
 * workflow (wraps generate_species_metadata.py, make_branch_subtrees.py, and
 * prepare_branch_ksrates_configs.py from ../run_ksrates_raw).
 *
 * Does NOT launch/stop the paralog Ks sqld server - that stays its own independent,
 * long-running process; point --ks_database_server_address_path at its already-running address file.
 * Does NOT execute any of the generated commands (ksrates init/orthologs-ks/paralogs-ks stay a
 * separate, later, separately-parallelized step) unless --run_init is explicitly turned on.
 */

LOG_OUTPUT = true

/*
 * Configuration: all parameters (params.*) and all process settings (executor, memory, CPUs, time,
 * concurrency, error handling, module loads) live in nextflow.config, not in this script. The
 * nextflow.config next to this script holds the defaults; each analysis puts its own
 * nextflow.config in the directory the pipeline is launched from, with only what differs.
 * Values that depend on other settings are derived here, since a config file cannot re-derive them
 * after another config file has overridden a base value.
 */

if (!params.base_dir) {
    error "base_dir is required, e.g. --base_dir /group/esb/cesen/1kp/ks_analysis/1kp"
}
base_dir = file(params.base_dir)

// Full 'ksrates' command line used by prepare_branch_ksrates_configs.py (init, orthologs-ks, paralogs-ks,
// check-*-db): the container, with the local ksrates checkout first on PYTHONPATH if ksrates_code_dir is set
ksrates_command = params.ksrates_command ?: "singularity exec -B /group/esb/cesen/1kp ${params.ksrates_code_dir ? '--env PYTHONPATH=' + params.ksrates_code_dir + ' ' : ''}${params.ksrates_sif} ksrates"
ksrates_container = params.ksrates_container ?: params.ksrates_sif

if (params.execute_real_pipelines && !params.expert_config) {
    error "execute_real_pipelines requires --expert_config (it needs use_ks_database_server " +
        "= yes to actually benefit from the already-populated database, rather than silently " +
        "recomputing everything from scratch)"
}

SCRIPT_DIR = "${workflow.projectDir}/../run_ksrates_raw"

// Generated once at startup (not per-task - nothing here depends on any process output), only
// when actually needed. Mirrors ks_analysis/test_subtree/nextflow.config's shape (singularity
// profile + per-process resource withName blocks) but points at our own local container/bind
// mounts instead of the public vibpsb/ksrates:latest image, and is rewritten fresh on every run
// so a changed --ks_database_server_address_path/--transcriptomes_dir is always reflected.
if (params.execute_real_pipelines) {
    bind_mounts = [base_dir.toString(), params.transcriptomes_dir]
    if (params.ksrates_code_dir) {
        bind_mounts << params.ksrates_code_dir
    }
    if (params.ks_database_server_address_path) {
        bind_mounts << file(params.ks_database_server_address_path).getParent().toString()
    }
    bind_options = bind_mounts.unique().collect { "-B ${it}" }.join(' ')

    env_line = params.ksrates_code_dir ? "env.PYTHONPATH = '${params.ksrates_code_dir}'" : ""
    file("${base_dir}/ksrates_run.nextflow.config").text = """
${env_line}
profiles {
    singularity {
        singularity.enabled = true
        singularity.cacheDir = '${params.base_dir}'
        singularity.autoMounts = true
        singularity.runOptions = "${bind_options}"
    }
}

process {
    container = '${ksrates_container}'
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
         ks_database_server_address_path:  ${params.ks_database_server_address_path ?: '(none)'}
         run_init:           ${params.run_init}
         test_max_datasets:  ${params.test_max_datasets ?: '(none - full run)'}
         execute_ortholog_jobs: ${params.execute_ortholog_jobs}
         execute_paralog_jobs:  ${params.execute_paralog_jobs}
         execute_real_pipelines: ${params.execute_real_pipelines}
         ksrates_main_nf:    ${params.ksrates_main_nf}
         """
         .stripIndent()


workflow {
    // 1) Get metadata per species
    generateSpeciesMetadata()
    // 2) Organize species in datasets
    generateDatasetList()
    // 3) Make config per species and writes cmd to be run
    prepareConfigs(generateSpeciesMetadata.out.species_metadata, generateDatasetList.out.dataset_list)
    prepareConfigs.out.summary.view()

    // Gated behind prepareConfigs.out.summary (not a bare Channel.fromPath on the command file)
    // so Nextflow waits for prepareConfigs to actually finish writing it first - these files
    // aren't Nextflow outputs themselves, just plain paths on disk.
    // 4) Compute and store ortholog data
    if (params.execute_ortholog_jobs) {
        ortholog_commands = prepareConfigs.out.summary
            .flatMap { file("${base_dir}/orthologs_ks_commands.sh").readLines() }
            .filter { line -> line.trim() && !line.startsWith('#') && line.trim() != 'set -euo pipefail' }
        runOrthologJob(ortholog_commands)
    }
    // 5) Compute and store paralog data
    if (params.execute_paralog_jobs) {
        paralog_commands = prepareConfigs.out.summary
            .flatMap { file("${base_dir}/paralogs_ks_commands.sh").readLines() }
            .filter { line -> line.trim() && !line.startsWith('#') && line.trim() != 'set -euo pipefail' }
        runParalogJob(paralog_commands)
    }
    // 6) Execute ksrates pipelines
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
        --out ${base_dir}/focal_branch_to_dataset.tsv \
        --datasets-out ${base_dir}/subtree_datasets.tsv \
        --score-table ${params.score_table} \
        --genome-busco-table ${params.genome_busco_table} \
        --max-a ${params.max_a} \
        --max-older ${params.max_older} \
        --older-clades ${params.older_clades} \
        --min-species ${params.min_species} \
        --busco-weight ${params.busco_weight} \
        --transrate-weight ${params.transrate_weight} \
        ${max_datasets_arg}
    DATASET_LIST=${base_dir}/subtree_datasets.tsv
    """
}


/*
 * Generates one ksrates config per (dataset, focal species) under base_dir/datasets/<hash>/,
 * plus the deduplicated command files (init_commands.sh, orthologs_ks_commands.sh, and
 * paralogs_ks_commands.sh if --ks_database_server_address_path is set). Depends on both processes above;
 * Nextflow waits for both automatically via the input channels below.
 */
process prepareConfigs {


    input:
        val species_metadata
        val dataset_list

    output:
        stdout emit: summary

    script:
    expert_arg = params.expert_config ? "--expert-config ${params.expert_config}" : ""
    ks_db_arg = params.ks_database_server_address_path ? "--ks-database-server-address-path ${params.ks_database_server_address_path}" : ""
    db_check_workers_arg = params.db_check_workers ? "--db-check-workers ${params.db_check_workers}" : ""
    run_init_arg = params.run_init ? "--run-init" : ""
    keep_relative_arg = params.keep_relative_paths ? "--keep-relative-paths" : ""
    """
    python3 ${SCRIPT_DIR}/prepare_branch_ksrates_configs.py \
        --datasets ${dataset_list} \
        --species-metadata ${species_metadata} \
        --out-dir ${base_dir} \
        ${expert_arg} \
        ${ks_db_arg} \
        --paralog-threads ${params.paralog_threads} \
        ${db_check_workers_arg} \
        --ortholog-threads ${params.ortholog_threads} \
        --ksrates-command "${ksrates_command}" \
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

    input:
        tuple val(dataset_id), val(focal_species), val(config_path), val(init_run_dir)

    script:
    """
    cd ${init_run_dir}
    nextflow run ${params.ksrates_main_nf} -profile singularity \
        --config ${config_path} \
        --expert ${params.expert_config} \
        -c ${base_dir}/ksrates_run.nextflow.config \
        -c ${params.ksrates_run_config}
    """
}
