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

if (!params.base_dir) {
    error "base_dir is required, e.g. --base_dir /group/esb/cesen/1kp/ks_analysis/1kp"
}
base_dir = file(params.base_dir)

// Same Python environment used by the existing launcher scripts (e.g.
// 1kp_wgms/launchers/launch_prepare_configs_test.sh) - these are plain host-side Python
// scripts, not run inside the ksrates container (only the "ksrates" subcommands prepareConfigs
// shells out to, via --ksrates-command, run inside the container).
PYTHON_SETUP = "module load python/x86_64/3.11.4 && source /group/esb/cesen/.venv/venv_python_11/bin/activate"

SCRIPT_DIR = "${workflow.projectDir}/../run_ksrates_raw"

log.info """
         batch_ksrates
         =============
         base_dir:          ${params.base_dir}
         tree:               ${params.tree}
         species_table:      ${params.species_table}
         expert_config:      ${params.expert_config ?: '(none)'}
         paralog_database:   ${params.paralog_database ?: '(none)'}
         run_init:           ${params.run_init}
         """
         .stripIndent()


workflow {
    generateSpeciesMetadata()
    generateDatasetList()
    prepareConfigs(generateSpeciesMetadata.out.species_metadata, generateDatasetList.out.dataset_list)
    prepareConfigs.out.summary.view()
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
        --transrate-weight ${params.transrate_weight}
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
