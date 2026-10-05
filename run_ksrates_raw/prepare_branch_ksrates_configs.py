#!/usr/bin/env python3
"""Prepare per-dataset ksrates configs and ortholog-pair workloads.

Each dataset is an independent ksrates analysis - its own tree and its own focal-species
selection - whether that's one of many systematically-generated subtrees of one big backbone
tree, or a handful of entirely unrelated, independently-curated datasets (e.g. "grasses",
"asterids") that happen to share some species and want to use the same centralized database.

The script expects two small TSV files:

1. A species metadata table with at least:
      species    fasta_filename
   Optional columns:
      latin_name    gff_filename

2. A dataset table with at least:
      dataset_id    newick_tree
   Optional columns:
      focal_species    target_species

By default, every leaf in a dataset's tree becomes a focal species. If the
``focal_species`` column is present, only those comma/semicolon/space-separated
species are used as focal species for that dataset.

The normal workflow is:

    python prepare_branch_ksrates_configs.py \
        --datasets datasets.tsv \
        --species-metadata species_metadata.tsv \
        --out-dir ksrates_setup \
        --max-outgroups 4 \
        --run-init

The wrapper can either trust ksrates' complete ``ortholog_pairs_<focal>`` output
(``--pair-scope init-pairs``), or derive a smaller pair set from selected
``ortholog_trios_<focal>`` rows (``--pair-scope target-trios``).

This script is the thin orchestrator: the actual work is split across
ksrates_batch_common.py (shared types/helpers), dataset_inputs.py (reading the two input
TSVs), ksrates_config_writer.py (writing configs), ksrates_init_runner.py (running
'ksrates init'), ortholog_pair_collection.py (ortholog pair dedup), and
paralog_species_dedup.py (paralog species dedup).
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Tuple

from dataset_inputs import ensure_species_known, read_datasets, read_species_metadata
from ksrates_batch_common import ConfigRecord, PairOccurrence, canonical_pair, fail, safe_name, write_command_file, write_tsv
from ksrates_config_writer import write_ksrates_config
from ksrates_init_runner import build_init_command, run_init
from ortholog_pair_collection import (
    build_ortholog_analysis_command,
    build_ortholog_command,
    collect_pairs_for_config,
    dataset_pair_rows,
    global_pair_rows,
)
from paralog_species_dedup import build_paralog_command, species_still_needs_paralog_ks


YES_NO = ("yes", "no")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate one ksrates config per selected focal species in each "
            "independent dataset, run/plan ksrates init, and collect the "
            "ortholog Ks species pairs reported by ksrates."
        )
    )
    parser.add_argument(
        "--datasets",
        required=True,
        type=Path,
        help="TSV with dataset_id and newick_tree columns.",
    )
    parser.add_argument(
        "--species-metadata",
        required=True,
        type=Path,
        help="TSV with species and fasta_filename columns.",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        type=Path,
        help="Directory for generated configs, init outputs, and pair manifests.",
    )
    parser.add_argument(
        "--expert-config",
        type=Path,
        default=None,
        help="Optional ksrates expert config passed to init/orthologs-ks/paralogs-ks.",
    )
    parser.add_argument(
        "--paralog-database",
        type=Path,
        default=None,
        help=(
            "Path to the paralog Ks sqld server's address file. When given, every generated "
            "config gets a ks_list_paralog_database_path line pointed at it, and a "
            "paralogs_ks_commands.sh is written with one deduplicated 'paralogs-ks' command per "
            "unique species across all datasets (mirrors the existing ortholog-pair dedup). The "
            "paired expert config (--expert-config) must set use_paralog_ks_database = yes."
        ),
    )
    parser.add_argument(
        "--paralog-threads",
        type=int,
        default=4,
        help="Thread count to place in generated paralogs-ks command file.",
    )
    parser.add_argument(
        "--paralog-check-workers",
        type=int,
        default=16,
        help=(
            "Number of concurrent 'ksrates check-paralog-db' calls to run while filtering "
            "already-populated species out of paralogs_ks_commands.sh. Each call is a small "
            "network round trip to the sqld server, not CPU-bound, so this can comfortably "
            "exceed the machine's core count."
        ),
    )
    parser.add_argument(
        "--ksrates-command",
        default="ksrates",
        help=(
            "Command used to launch ksrates. May contain spaces, e.g. "
            "'apptainer exec image.sif ksrates'."
        ),
    )
    parser.add_argument(
        "--run-init",
        action="store_true",
        help="Actually run 'ksrates init' for every generated config.",
    )
    parser.add_argument(
        "--pair-scope",
        choices=("init-pairs", "target-trios"),
        default="init-pairs",
        help=(
            "init-pairs trusts ortholog_pairs_<focal> as written by ksrates. "
            "target-trios derives only focal-target-outgroup pairs from "
            "ortholog_trios_<focal> rows whose sister/comparison species is one of the "
            "species this dataset's target_species column says to actually compute "
            "ortholog pairs against."
        ),
    )
    parser.add_argument(
        "--max-outgroups",
        type=int,
        default=4,
        help="Value for ksrates max_number_outgroups.",
    )
    parser.add_argument(
        "--consensus-mode",
        default="mean among outgroups",
        help="Value for ksrates consensus_mode_for_multiple_outgroups.",
    )
    parser.add_argument(
        "--paranome",
        choices=YES_NO,
        default="no",
        help="Value for ksrates paranome setting in generated configs.",
    )
    parser.add_argument(
        "--collinearity",
        choices=YES_NO,
        default="no",
        help="Value for ksrates collinearity setting in generated configs.",
    )
    parser.add_argument(
        "--reciprocal-retention",
        choices=YES_NO,
        default="no",
        help="Value for ksrates reciprocal_retention setting in generated configs.",
    )
    parser.add_argument(
        "--gff-feature",
        default="mrna",
        help="Value for ksrates gff_feature.",
    )
    parser.add_argument(
        "--gff-attribute",
        default="id",
        help="Value for ksrates gff_attribute.",
    )
    parser.add_argument(
        "--num-bootstrap-iterations",
        type=int,
        default=200,
        help="Value for ksrates num_bootstrap_iterations.",
    )
    parser.add_argument(
        "--max-ks-paralogs",
        default="5",
        help="Value for ksrates max_ks_paralogs.",
    )
    parser.add_argument(
        "--max-ks-orthologs",
        default="10",
        help="Value for ksrates max_ks_orthologs.",
    )
    parser.add_argument(
        "--ortholog-threads",
        type=int,
        default=4,
        help="Thread count to place in generated orthologs-ks command file.",
    )
    parser.add_argument(
        "--keep-relative-paths",
        action="store_true",
        help=(
            "Keep FASTA/GFF paths exactly as provided. By default, relative "
            "metadata paths are resolved relative to the metadata file."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.datasets = args.datasets.resolve()
    args.species_metadata = args.species_metadata.resolve()
    args.out_dir = args.out_dir.resolve()
    if args.expert_config is not None:
        args.expert_config = args.expert_config.resolve()
    if args.paralog_database is not None:
        args.paralog_database = args.paralog_database.resolve()

    metadata = read_species_metadata(args.species_metadata, args.keep_relative_paths)
    datasets = read_datasets(args.datasets)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    # Each dataset's configs/, per-dataset ortholog_pairs.tsv, and per-focal init/pipeline
    # working directory all live together under datasets/<dataset_id>/ - one self-contained
    # home per dataset, instead of parallel top-level dirs sharded by dataset.
    datasets_root = args.out_dir / "datasets"
    database_root = args.out_dir / "databases"
    ortholog_work_dir = args.out_dir / "ortholog_runs"
    paralog_work_dir = args.out_dir / "paralog_runs"
    database_root.mkdir(parents=True, exist_ok=True)
    ortholog_work_dir.mkdir(parents=True, exist_ok=True)
    if args.paralog_database is not None:
        paralog_work_dir.mkdir(parents=True, exist_ok=True)

    peak_database_path = database_root / "ortholog_peak_db.tsv"
    ks_list_database_path = database_root / "ortholog_ks_list_db.tsv"

    configs: List[ConfigRecord] = []
    for dataset in datasets:
        ensure_species_known(dataset.leaves, metadata, f"dataset {dataset.dataset_id!r} Newick")
        ensure_species_known(dataset.focal_species, metadata, f"dataset {dataset.dataset_id!r} focal list")
        unknown_focals = sorted(set(dataset.focal_species) - set(dataset.leaves))
        if unknown_focals:
            fail(
                f"dataset {dataset.dataset_id!r} focal species are not in the dataset: "
                f"{', '.join(unknown_focals)}"
            )
        if dataset.target_species:
            ensure_species_known(
                dataset.target_species,
                metadata,
                f"dataset {dataset.dataset_id!r} target list",
            )
            unknown_targets = sorted(set(dataset.target_species) - set(dataset.leaves))
            if unknown_targets:
                fail(
                    f"dataset {dataset.dataset_id!r} target species are not in the dataset: "
                    f"{', '.join(unknown_targets)}"
                )

        dataset_dir_name = safe_name(dataset.dataset_id)
        dataset_dir = datasets_root / dataset_dir_name
        for focal_species in dataset.focal_species:
            config_path = (
                dataset_dir
                / "configs"
                / f"config_{dataset_dir_name}_{safe_name(focal_species)}.txt"
            )
            # Also the real Stage-4 main.nf pipeline run's own working directory for this
            # focal species, not just a throwaway init-probe location - ksrates' own hardcoded
            # rate_adjustment/<focal>/ subdirectory lands inside here either way.
            init_run_dir = dataset_dir / safe_name(focal_species)
            write_ksrates_config(
                config_path,
                dataset,
                focal_species,
                metadata,
                args,
                peak_database_path,
                ks_list_database_path,
                paralog_database_path=args.paralog_database,
            )
            configs.append(
                ConfigRecord(
                    dataset_id=dataset.dataset_id,
                    focal_species=focal_species,
                    config_path=config_path,
                    init_run_dir=init_run_dir,
                    target_species=dataset.target_species,
                )
            )

    if args.run_init:
        for config in configs:
            run_init(args, config)

    init_commands = [(config.init_run_dir, [build_init_command(args, config.config_path)]) for config in configs]
    write_command_file(args.out_dir / "init_commands.sh", init_commands)

    all_occurrences: List[PairOccurrence] = []
    missing_pair_outputs: List[ConfigRecord] = []
    for config in configs:
        occurrences = collect_pairs_for_config(args, config)
        if occurrences:
            all_occurrences.extend(occurrences)
        else:
            missing_pair_outputs.append(config)

    write_tsv(
        args.out_dir / "config_manifest.tsv",
        (
            "dataset_id",
            "focal_species",
            "config_path",
            "init_run_dir",
            "target_species",
        ),
        (
            (
                config.dataset_id,
                config.focal_species,
                config.config_path,
                config.init_run_dir,
                ",".join(config.target_species),
            )
            for config in configs
        ),
    )

    pair_to_occurrences: Dict[Tuple[str, str], List[PairOccurrence]] = defaultdict(list)
    dataset_to_occurrences: Dict[str, List[PairOccurrence]] = defaultdict(list)
    for occurrence in all_occurrences:
        pair = canonical_pair(occurrence.species_1, occurrence.species_2)
        pair_to_occurrences[pair].append(occurrence)
        dataset_to_occurrences[occurrence.dataset_id].append(occurrence)

    # Global cross-dataset manifests live under ortholog_runs/, next to the pair working files
    # and wgd output they inform - not scattered at the out-dir root separate from the work
    # they drive.
    write_tsv(
        ortholog_work_dir / "ortholog_pairs_by_source.tsv",
        (
            "species_1",
            "species_2",
            "dataset_id",
            "focal_species",
            "config_path",
            "init_run_dir",
            "source_kind",
            "source_file",
        ),
        (
            (
                occurrence.species_1,
                occurrence.species_2,
                occurrence.dataset_id,
                occurrence.focal_species,
                occurrence.config_path,
                occurrence.init_run_dir,
                occurrence.source_kind,
                occurrence.source_file,
            )
            for occurrence in sorted(
                all_occurrences,
                key=lambda item: (
                    item.dataset_id,
                    item.focal_species.lower(),
                    item.species_1.lower(),
                    item.species_2.lower(),
                ),
            )
        ),
    )

    write_tsv(
        ortholog_work_dir / "ortholog_pairs.tsv",
        (
            "species_1",
            "species_2",
            "datasets",
            "focal_species",
            "chosen_config",
            "chosen_init_run_dir",
        ),
        global_pair_rows(pair_to_occurrences),
    )

    for dataset_id, occurrences in dataset_to_occurrences.items():
        dataset_pair_map: Dict[Tuple[str, str], List[PairOccurrence]] = defaultdict(list)
        for occurrence in occurrences:
            dataset_pair_map[canonical_pair(occurrence.species_1, occurrence.species_2)].append(occurrence)
        write_tsv(
            datasets_root / safe_name(dataset_id) / "ortholog_pairs.tsv",
            (
                "species_1",
                "species_2",
                "dataset_id",
                "focal_species",
                "chosen_config",
            ),
            dataset_pair_rows(dataset_pair_map),
        )

    # Each entry chains orthologs-ks (raw Ks computation, local files only) with
    # orthologs-analysis (bootstrap peak estimation, writes to the shared peak/Ks-list TSVs) via
    # "&&", so the database actually gets populated - orthologs-ks alone never touches it. The
    # analysis step is scoped to a single-pair TSV (not the full per-focal ortholog_pairs file)
    # so it only ever processes the one pair just computed, not every pair that focal needs.
    ortholog_pairs_dir = ortholog_work_dir / "pairs"
    ortholog_commands = []
    for (species_1, species_2), occurrences in sorted(
        pair_to_occurrences.items(),
        key=lambda item: (item[0][0].lower(), item[0][1].lower()),
    ):
        chosen = sorted(
            occurrences,
            key=lambda item: (item.dataset_id, item.focal_species.lower(), str(item.config_path)),
        )[0]
        pair_tsv_path = ortholog_pairs_dir / f"{species_1}_{species_2}.tsv"
        write_tsv(pair_tsv_path, ("Species1", "Species2"), [(species_1, species_2)])
        ortholog_commands.append(
            (
                ortholog_work_dir,
                [
                    build_ortholog_command(args, chosen.config_path, species_1, species_2),
                    build_ortholog_analysis_command(args, chosen.config_path, pair_tsv_path),
                ],
            )
        )
    write_command_file(args.out_dir / "orthologs_ks_commands.sh", ortholog_commands)

    # Species-level dedup for paralog population, mirroring the pair-level dedup above: the same
    # real species can appear as focal_species in many different datasets (since species is
    # drawn from one global species_metadata.tsv, the string itself is already a stable
    # cross-dataset identity - no latin-name translation needed, unlike per-dataset informal
    # names elsewhere). Without this, feeding every (dataset, focal) config to paralogs-ks-multi
    # would launch one paralogs-ks call per dataset a species appears in, not once per species -
    # exactly the duplicate-computation race this preparation step exists to avoid.
    if args.paralog_database is not None:
        species_to_occurrences: Dict[str, List[ConfigRecord]] = defaultdict(list)
        for config in configs:
            species_to_occurrences[config.focal_species].append(config)

        chosen_config_by_species: Dict[str, Path] = {
            species: sorted(
                occurrences,
                key=lambda item: (item.dataset_id, str(item.config_path)),
            )[0].config_path
            for species, occurrences in species_to_occurrences.items()
        }

        # Each check is a small network round trip to the sqld server, not CPU-bound, so these
        # run concurrently - run serially this would dominate wall-clock time at ~2000-species
        # scale even though every individual check is cheap.
        ordered_species = sorted(species_to_occurrences, key=str.lower)
        with ThreadPoolExecutor(max_workers=args.paralog_check_workers) as executor:
            check_results = list(
                executor.map(
                    lambda species: species_still_needs_paralog_ks(args, chosen_config_by_species[species]),
                    ordered_species,
                )
            )
        queued_species: Dict[str, bool] = dict(zip(ordered_species, check_results))

        paralog_commands = []
        for species in ordered_species:
            if queued_species[species]:
                paralog_commands.append(
                    (paralog_work_dir, [build_paralog_command(args, chosen_config_by_species[species])])
                )
        write_command_file(args.out_dir / "paralogs_ks_commands.sh", paralog_commands)

        write_tsv(
            paralog_work_dir / "paralog_species.tsv",
            ("species", "datasets", "chosen_config", "queued"),
            (
                (
                    species,
                    ",".join(sorted({occ.dataset_id for occ in occurrences}, key=str.lower)),
                    sorted(occurrences, key=lambda item: (item.dataset_id, str(item.config_path)))[0].config_path,
                    "yes" if queued_species[species] else "no",
                )
                for species, occurrences in sorted(species_to_occurrences.items(), key=lambda item: item[0].lower())
            ),
        )

        already_populated = len(species_to_occurrences) - len(paralog_commands)
        print(f"Collected {len(species_to_occurrences)} unique paralog species "
              f"({len(paralog_commands)} queued, {already_populated} already fully populated)")
        print(f"Wrote paralog command file: {args.out_dir / 'paralogs_ks_commands.sh'}")

    print(f"Wrote {len(configs)} ksrates configs to {datasets_root}")
    print(f"Wrote init command file: {args.out_dir / 'init_commands.sh'}")
    if pair_to_occurrences:
        print(f"Collected {len(pair_to_occurrences)} unique ortholog pairs")
        print(f"Wrote ortholog command file: {args.out_dir / 'orthologs_ks_commands.sh'}")
    else:
        print(
            "No ortholog pairs were collected yet. Run with --run-init, or execute "
            f"{args.out_dir / 'init_commands.sh'} and rerun this script to collect outputs."
        )
    if missing_pair_outputs:
        print(
            f"WARNING: {len(missing_pair_outputs)} configs had no collectable ksrates "
            "pair/trio output.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
