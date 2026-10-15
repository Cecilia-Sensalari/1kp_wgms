#!/usr/bin/env nextflow

/*
 * batch_ksrates: generates species_metadata.tsv, the branch table (focal_branch_to_dataset.tsv), the dataset list (subtree_datasets.tsv), and
 * per-dataset ksrates configs and the tables of jobs for the centralized paralog/ortholog Ks database
 * workflow (wraps generate_species_metadata.py, make_branch_subtrees.py, and
 * prepare_branch_ksrates_configs.py from ../run_ksrates_raw).
 *
 * Does NOT launch/stop the paralog Ks sqld server - that stays its own independent,
 * long-running process; point ks_database_server_address_path at its already-running address file.
 * Runs ksrates init (one job per config without init output) only with run_init = true, and the Ks jobs (paralogs-ks,
 * orthologs-ks and orthologs-analysis) only with execute_paralog_jobs and execute_ortholog_jobs = true: each job takes
 * its config, species or pair from a table written by prepare_branch_ksrates_configs.py, and the ksrates command is
 * written once, in its process.
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

if (!params.container_mounts) {
    error "container_mounts is required: the folders to mount in the containers, e.g. container_mounts = ['/group/esb/cesen/1kp']"
}
if (!params.base_dir) {
    error "base_dir is required: set it in the nextflow.config of the analysis, e.g. base_dir = '/group/esb/cesen/1kp/ks_analysis/1kp'"
}
base_dir = file(params.base_dir)

// Table with the A, B, C and D species of every branch: next to the dataset list (written there by generateDatasetList)
branch_table_file = "${params.dataset_list ? file(params.dataset_list).parent : base_dir}/focal_branch_to_dataset.tsv"

// Folders mounted in every container the wrapper starts: the ksrates command of the processes below and the containers of
// the per-dataset ksrates pipelines (stage 4). One list, container_mounts, for both.
mount_options = params.container_mounts.collect { "-B ${it}" }.join(' ')

// Full 'ksrates' command line used by the processes that run ksrates (init, orthologs-ks, orthologs-analysis, paralogs-ks) and by
// prepare_branch_ksrates_configs.py (check-*-db): the container, with the local ksrates checkout first on PYTHONPATH if ksrates_code_dir is set
container_exec = "singularity exec ${mount_options} ${params.ksrates_code_dir ? '--env PYTHONPATH=' + params.ksrates_code_dir + ' ' : ''}${params.ksrates_sif}"
ksrates_command = params.ksrates_command ?: "${container_exec} ksrates"
ksrates_container = params.ksrates_container ?: params.ksrates_sif

// Options of prepare_branch_ksrates_configs.py that do not depend on the process inputs. They are shared by the
// processes prepareConfigs (before ksrates init) and collectOrthologPairs (after it), which run the same script
prepare_args = [
    params.expert_config ? "--expert-config ${params.expert_config}" : "",
    params.ks_database_server_address_path ? "--ks-database-server-address-path ${params.ks_database_server_address_path}" : "",
    "--paralog-threads ${params.paralog_threads}",
    params.db_check_workers ? "--db-check-workers ${params.db_check_workers}" : "",
    "--ortholog-threads ${params.ortholog_threads}",
    "--ksrates-command \"${ksrates_command}\"",
    "--pair-scope ${params.pair_scope}",
    "--max-outgroups ${params.max_outgroups}",
    "--consensus-mode \"${params.consensus_mode}\"",
    "--paranome ${params.paranome}",
    "--reciprocal-retention ${params.reciprocal_retention}",
    "--gff-feature ${params.gff_feature}",
    "--gff-attribute ${params.gff_attribute}",
    "--num-bootstrap-iterations ${params.num_bootstrap_iterations}",
    "--max-ks-paralogs ${params.max_ks_paralogs}",
    "--max-ks-orthologs ${params.max_ks_orthologs}",
    params.test_last_datasets ? "--last-datasets ${params.test_last_datasets}" : "",
    params.keep_relative_paths ? "--keep-relative-paths" : "",
].findAll { it }.join(' ')

if (params.execute_ksrates_pipelines && !params.expert_config) {
    error "execute_ksrates_pipelines requires expert_config (it needs use_ks_database_server " +
        "= yes to actually benefit from the already-populated database, rather than silently " +
        "recomputing everything from scratch)"
}

SCRIPT_DIR = "${workflow.projectDir}/../run_ksrates_raw"

// Generated once at startup (not per-task - nothing here depends on any process output), only
// when actually needed. Mirrors ks_analysis/test_subtree/nextflow.config's shape (singularity
// profile + per-process resource withName blocks) but points at our own local container/bind
// mounts instead of the public vibpsb/ksrates:latest image, and is rewritten fresh on every run
// so a changed container_mounts, ksrates_code_dir or container is always reflected.
if (params.execute_ksrates_pipelines) {
    // Resources of the processes of the ksrates pipeline: the threads follow paralog_threads and ortholog_threads, the memory
    // comes from ksrates_pipeline_memory (process name -> memory)
    process_resources = params.ksrates_pipeline_memory.collect { process_name, process_memory ->
        def threads = process_name == 'wgdParalogs' ? "        cpus = ${params.paralog_threads}\n" : (process_name == 'wgdOrthologs' ? "        cpus = ${params.ortholog_threads}\n" : "")
        "    withName: '${process_name}' {\n${threads}        memory = '${process_memory}'\n    }"
    }.join('\n')
    env_line = params.ksrates_code_dir ? "env.PYTHONPATH = '${params.ksrates_code_dir}'" : ""
    file("${base_dir}/ksrates_run.nextflow.config").text = """
${env_line}
profiles {
    singularity {
        singularity.enabled = true
        singularity.cacheDir = '${params.base_dir}'
        singularity.autoMounts = true
        singularity.runOptions = "${mount_options}"
    }
}

executor {
    name = '${params.ksrates_pipeline_executor}'
}

process {
    container = '${ksrates_container}'
${process_resources}
}
""".stripIndent()
}



// Plain-language summary of this run, printed at the start
dataset_source = params.dataset_list ? "existing table ${params.dataset_list} (not regenerated)" : "generated in this run from the tree ${params.tree}"
if (!params.dataset_list && params.test_max_datasets) {
    dataset_source += ", keeping only the first ${params.test_max_datasets}"
}
dataset_use = params.test_last_datasets ? "only the last ${params.test_last_datasets} datasets of the table (test)" : "all datasets"
yes_no = { flag -> flag ? 'yes' : 'no' }

log.info """
         batch_ksrates
         =============
         Output directory:  ${params.base_dir}
         Species table:     ${params.species_table}
         Datasets from:     ${dataset_source}
         Datasets used:     ${dataset_use}
         Ks database:       ${params.ks_database_server_address_path ?: '(none)'}
         Expert config:     ${params.expert_config ?: '(none)'}

         Stages of this run:
           species metadata and configs    yes (always)
           ksrates init                    ${yes_no(params.run_init)}
           ortholog Ks jobs                ${yes_no(params.execute_ortholog_jobs)}
           paralog Ks jobs                 ${yes_no(params.execute_paralog_jobs)}
           branch trees                    ${yes_no(params.draw_branch_trees)}
           per-dataset ksrates pipelines   ${yes_no(params.execute_ksrates_pipelines)}${params.execute_ksrates_pipelines ? ' (script: ' + params.ksrates_main_nf + ')' : ''}
         """
         .stripIndent()


workflow {
    // 1) Get metadata per species
    generateSpeciesMetadata()
    // 2) Organize species in datasets, unless an already generated dataset list is given
    if (params.dataset_list) {
        if (!file(params.dataset_list).exists()) {
            error "dataset_list does not exist: ${params.dataset_list}"
        }
        dataset_list_ch = Channel.value(file(params.dataset_list).toString())
    } else {
        generateDatasetList()
        dataset_list_ch = generateDatasetList.out.dataset_list
    }
    // 3) Colours of the divergence lines of every run, from the branch table written together with the dataset list
    makeFocalColorings(dataset_list_ch)
    // 4) Make one config per run and the tables of the Ks jobs still to run (the jobs are started from those tables below)
    prepareConfigs(generateSpeciesMetadata.out.species_metadata, dataset_list_ch, makeFocalColorings.out.colorings)
    prepareConfigs.out.summary.view()

    // 5a) ksrates init, only for the configs whose init output is missing (init_done = no in config_manifest.tsv), and
    // then the second pass of the config script, which collects the ortholog pairs from the init output. The paralog
    // jobs below do not need init and start right after prepareConfigs.
    if (params.run_init) {

        init_jobs = prepareConfigs.out.summary
            .flatMap { file("${base_dir}/config_manifest.tsv").readLines().drop(1) }
            .map { line -> line.split('\t') }
            .filter { fields -> fields[5] == 'no' }
            .map { fields -> tuple(fields[0], fields[1], fields[2], fields[3]) }    // dataset, focal species, config, init dir
        runInit(init_jobs)

        // One item once every init job has ended (the extra item keeps the list from being empty when none ran)
        init_done = runInit.out.done.mix(Channel.of('none')).toList()
        collectOrthologPairs(generateSpeciesMetadata.out.species_metadata, dataset_list_ch, makeFocalColorings.out.colorings, init_done)
        collectOrthologPairs.out.summary.view()
        pair_tables_ready = collectOrthologPairs.out.summary
        
    } else {
        pair_tables_ready = prepareConfigs.out.summary
    }

    // 5b) The tree of every tested branch with the focal branch highlighted, drawn next to the plots of each run (after init,
    // which creates the folders). One job per dataset; it does not affect the Ks values.
    if (params.draw_branch_trees) {
        branch_tree_jobs = pair_tables_ready
            .flatMap { file("${base_dir}/config_manifest.tsv").readLines().drop(1) }
            .map { line -> line.split('\t')[0] }
            .unique()
        drawBranchTrees(branch_tree_jobs)
    }

    // The jobs of the next stages are read from tables that prepareConfigs writes. Each read starts from
    // prepareConfigs.out.summary, so that Nextflow waits for prepareConfigs to finish writing them first
    // (the tables are not Nextflow outputs, just plain files on disk).
    // 6) Compute and store ortholog data: one job per pair marked yes in the column "queued" of ortholog_pairs.tsv
    if (params.execute_ortholog_jobs) {
        ortholog_jobs = pair_tables_ready
            .flatMap { file("${base_dir}/ortholog_runs/ortholog_pairs.tsv").readLines().drop(1) }
            .map { line -> line.split('\t') }
            .filter { fields -> fields[6] == 'yes' }
            .map { fields -> tuple(fields[0], fields[1], fields[4]) }    // species_1, species_2, chosen config
        runOrthologJob(ortholog_jobs)
        // One item once every ortholog job has ended (the extra item keeps the list from being empty when no job ran)
        ortholog_done = runOrthologJob.out.done.mix(Channel.of('none')).toList()
    }
    // 7) Compute and store paralog data: one job per species marked yes in the column "queued" of paralog_species.tsv
    if (params.execute_paralog_jobs) {
        if (!params.ks_database_server_address_path) {
            error "execute_paralog_jobs requires ks_database_server_address_path: the paralog jobs are only listed when the Ks database is used"
        }
        paralog_jobs = prepareConfigs.out.summary
            .flatMap { file("${base_dir}/paralog_runs/paralog_species.tsv").readLines().drop(1) }
            .map { line -> line.split('\t') }
            .filter { fields -> fields[3] == 'yes' }
            .map { fields -> tuple(fields[0], fields[2]) }               // species, chosen config
        runParalogJob(paralog_jobs)
        paralog_done = runParalogJob.out.done.mix(Channel.of('none')).toList()
    }
    // 8) Execute ksrates pipelines. If init and/or the ortholog and paralog jobs are run in the same launch, wait until
    // all of them have ended, so that the pipelines find the init output and the Ks values in the database.
    if (params.execute_ksrates_pipelines) {
        pipelines_ready = prepareConfigs.out.summary
        if (params.run_init) {
            pipelines_ready = pipelines_ready.combine(init_done)
        }
        if (params.execute_ortholog_jobs) {
            pipelines_ready = pipelines_ready.combine(ortholog_done)
        }
        if (params.execute_paralog_jobs) {
            pipelines_ready = pipelines_ready.combine(paralog_done)
        }
        ksrates_pipeline_rows = pipelines_ready
            .flatMap { file("${base_dir}/config_manifest.tsv").readLines().drop(1) }
            .map { line -> line.split('\t') }
            .map { fields -> tuple(fields[0], fields[1], fields[2], fields[3]) }
        runKsratesPipeline(ksrates_pipeline_rows)
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
        --genome-manifest ${params.genome_manifest} \
        --genome-gff-settings ${params.genome_gff_settings} \
        --collinearity ${params.collinearity} \
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
 * Writes focal_colorings.tsv in base_dir: the colours of the divergence lines in the plot of every focal species of every
 * dataset, by the role of each divergence node in the tested branch (see make_focal_colorings.py). It reads
 * focal_branch_to_dataset.tsv, which has to be in the folder of the dataset list (generateDatasetList writes both there).
 * Takes seconds, so it runs on every launch.
 */
process makeFocalColorings {


    input:
        val dataset_list

    output:
        env COLORINGS, emit: colorings

    script:
    branch_table = "${file(dataset_list).parent}/focal_branch_to_dataset.tsv"
    """
    python3 ${SCRIPT_DIR}/make_focal_colorings.py --branches ${branch_table} --out ${base_dir}/focal_colorings.tsv --palette "${params.divergence_palette}"
    COLORINGS=${base_dir}/focal_colorings.tsv
    """
}


/*
 * First pass of prepare_branch_ksrates_configs.py: generates one ksrates config per (dataset, focal species) under
 * base_dir/datasets/<hash>/, config_manifest.tsv (with a column init_done: yes if ksrates init has already
 * written its output for that config), and the tables of the Ks jobs still to run, with a column "queued":
 * ortholog_runs/ortholog_pairs.tsv (from the init output that exists at this point) and
 * paralog_runs/paralog_species.tsv (if ks_database_server_address_path is set).
 * Depends on both processes above; Nextflow waits for both automatically via the input channels below.
 */
process prepareConfigs {


    input:
        val species_metadata
        val dataset_list
        val colorings

    output:
        stdout emit: summary

    script:
    """
    echo "=== config pass 1: the configs and the tables of the Ks jobs ==="
    python3 ${SCRIPT_DIR}/prepare_branch_ksrates_configs.py --datasets ${dataset_list} --colorings ${colorings} --species-metadata ${species_metadata} --out-dir ${base_dir} ${prepare_args}
    """
}


/*
 * Runs ksrates init for one (dataset, focal species) config, in its own datasets/<hash>/<focal>/ directory. The log
 * goes to init.log there. Only started (with run_init = true) for the configs marked init_done = no in config_manifest.tsv.
 */
process runInit {

    tag "${dataset_id}/${focal_species}"

    input:
        tuple val(dataset_id), val(focal_species), val(config_path), val(init_run_dir)

    output:
        val 'done', emit: done

    script:
    expert_arg = params.expert_config ? "--expert ${params.expert_config}" : ""
    """
    mkdir -p ${init_run_dir}
    cd ${init_run_dir}
    ${ksrates_command} init ${config_path} ${expert_arg} > init.log 2>&1
    """
}


/*
 * Second pass of prepare_branch_ksrates_configs.py, after every ksrates init has ended: the same script and options,
 * now finding the init output of the new configs, so that it collects their ortholog pairs and rewrites
 * ortholog_pairs.tsv with the pairs still to run. The input init_done only makes this process wait for the init jobs.
 */
process collectOrthologPairs {


    input:
        val species_metadata
        val dataset_list
        val colorings
        val init_done

    output:
        stdout emit: summary

    script:
    """
    echo "=== config pass 2, after ksrates init: the same, now with the ortholog pairs of the new init output ==="
    python3 ${SCRIPT_DIR}/prepare_branch_ksrates_configs.py --datasets ${dataset_list} --colorings ${colorings} --species-metadata ${species_metadata} --out-dir ${base_dir} ${prepare_args}
    """
}


/*
 * Draws the trees of the tested branches of one dataset (draw_branch_trees.py), inside the ksrates container, next to the
 * plots of each run: <run folder>/rate_adjustment/<focal>/tree_<focal>.b<branch>.pdf, with the focal branch in the colour of
 * "inside A" and the species names in the colour of their node (A, B, C; D in black). One job per dataset, since starting
 * the container costs more than drawing. Not needed for the Ks values, so a failure is ignored (see nextflow.config).
 */
process drawBranchTrees {

    tag "${dataset_id}"

    input:
        val dataset_id

    script:
    """
    ${container_exec} python3 ${SCRIPT_DIR}/draw_branch_trees.py \
        --dataset ${dataset_id} \
        --branches ${branch_table_file} \
        --manifest ${base_dir}/config_manifest.tsv \
        --species-metadata ${base_dir}/species_metadata.tsv \
        --palette "${params.divergence_palette}"
    """
}


/*
 * Runs paralogs-ks for one species, in paralog_runs/ (where ksrates writes paralog_distributions/). One job per
 * species marked queued in paralog_species.tsv. The ksrates output goes to
 * paralog_runs/paralog_ksrates_logs/<species>.log; the log of an earlier run is kept, renamed with the time of
 * its last modification.
 */
process runParalogJob {

    tag "${species}"

    input:
        tuple val(species), val(config_path)

    output:
        val 'done', emit: done

    script:
    expert_arg = params.expert_config ? "--expert ${params.expert_config}" : ""
    run_dir = "${base_dir}/paralog_runs"
    log_dir = "${run_dir}/paralog_ksrates_logs"
    """
    mkdir -p ${log_dir}
    log=${log_dir}/${species}.log
    if [ -e "\$log" ]; then mv "\$log" "\${log%.log}.\$(date -r "\$log" +%Y%m%d_%H%M%S).log"; fi
    cd ${run_dir}
    ${ksrates_command} paralogs-ks ${config_path} ${expert_arg} --n-threads ${params.paralog_threads} > "\$log" 2>&1
    """
}


/*
 * Runs the Ks estimation of one ortholog species pair, in ortholog_runs/: orthologs-ks (raw Ks values) and then
 * orthologs-analysis (bootstrap peak estimation, which is what stores the peak and the Ks list; orthologs-ks alone
 * never does). One job per pair marked queued in ortholog_pairs.tsv. The analysis is pointed at a one-pair table
 * (ortholog_runs/pairs/<A>_<B>.tsv, written here), so that it only processes this pair and not every pair of the
 * focal species. The ksrates output of both commands goes to ortholog_runs/ortholog_ksrates_logs/<A>_<B>.log; the
 * log of an earlier run is kept, renamed with the time of its last modification.
 */
process runOrthologJob {

    tag "${species_1}_${species_2}"

    input:
        tuple val(species_1), val(species_2), val(config_path)

    output:
        val 'done', emit: done

    script:
    expert_arg = params.expert_config ? "--expert ${params.expert_config}" : ""
    run_dir = "${base_dir}/ortholog_runs"
    log_dir = "${run_dir}/ortholog_ksrates_logs"
    pair_table = "${run_dir}/pairs/${species_1}_${species_2}.tsv"
    """
    mkdir -p ${log_dir} ${run_dir}/pairs
    printf 'Species1\\tSpecies2\\n${species_1}\\t${species_2}\\n' > ${pair_table}
    log=${log_dir}/${species_1}_${species_2}.log
    if [ -e "\$log" ]; then mv "\$log" "\${log%.log}.\$(date -r "\$log" +%Y%m%d_%H%M%S).log"; fi
    cd ${run_dir}
    {
        ${ksrates_command} orthologs-ks ${config_path} ${expert_arg} ${species_1} ${species_2} --n-threads ${params.ortholog_threads} &&
        ${ksrates_command} orthologs-analysis ${config_path} ${expert_arg} --ortholog-pairs ${pair_table}
    } > "\$log" 2>&1
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
process runKsratesPipeline {

    tag "${dataset_id}/${focal_species}"

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
