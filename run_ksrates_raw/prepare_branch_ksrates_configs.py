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
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shlex
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


YES_NO = ("yes", "no")


@dataclass(frozen=True)
class SpeciesRecord:
    species: str
    latin_name: str
    fasta_filename: str
    gff_filename: str = ""


@dataclass(frozen=True)
class DatasetRecord:
    dataset_id: str
    newick_tree: str
    leaves: Tuple[str, ...]
    focal_species: Tuple[str, ...]
    target_species: Tuple[str, ...]


@dataclass(frozen=True)
class ConfigRecord:
    dataset_id: str
    focal_species: str
    config_path: Path
    init_run_dir: Path
    target_species: Tuple[str, ...]


@dataclass(frozen=True)
class PairOccurrence:
    species_1: str
    species_2: str
    dataset_id: str
    focal_species: str
    config_path: Path
    init_run_dir: Path
    source_file: Path
    source_kind: str


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


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def read_table(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open(newline="") as handle:
        first_line = handle.readline()
        if not first_line:
            fail(f"{path} is empty")
        delimiter = "\t" if "\t" in first_line else ","
        handle.seek(0)
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            fail(f"{path} has no header")
        fieldnames = [name.strip().lstrip("\ufeff") for name in reader.fieldnames]
        rows: List[Dict[str, str]] = []
        for raw_row in reader:
            row = {}
            for raw_name, value in raw_row.items():
                if raw_name is None:
                    continue
                name = raw_name.strip().lstrip("\ufeff")
                row[name] = (value or "").strip()
            rows.append(row)
    return fieldnames, rows


def column_name(fieldnames: Sequence[str], aliases: Sequence[str], required: bool) -> Optional[str]:
    by_normalized = {normalize_column_name(name): name for name in fieldnames}
    for alias in aliases:
        normalized = normalize_column_name(alias)
        if normalized in by_normalized:
            return by_normalized[normalized]
    if required:
        fail(f"Missing required column; expected one of: {', '.join(aliases)}")
    return None


def normalize_column_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", name.lower())


def normalize_data_path(value: str, base_dir: Path, keep_relative: bool) -> str:
    if not value:
        return ""
    path = Path(value)
    if keep_relative or path.is_absolute():
        return value
    return str((base_dir / path).resolve())


def read_species_metadata(path: Path, keep_relative_paths: bool) -> Dict[str, SpeciesRecord]:
    fieldnames, rows = read_table(path)
    species_col = column_name(fieldnames, ("species", "species_id", "species_name"), True)
    fasta_col = column_name(
        fieldnames,
        ("fasta_filename", "fasta", "cds_fasta", "transcript_fasta", "sequence_file"),
        True,
    )
    latin_col = column_name(fieldnames, ("latin_name", "latin_names", "scientific_name"), False)
    gff_col = column_name(fieldnames, ("gff_filename", "gff", "gff3", "annotation_gff"), False)

    base_dir = path.parent.resolve()
    metadata: Dict[str, SpeciesRecord] = {}
    assert species_col is not None
    assert fasta_col is not None
    for row_number, row in enumerate(rows, start=2):
        species = row.get(species_col, "").strip()
        if not species:
            fail(f"{path}:{row_number} has an empty species value")
        if species in metadata:
            fail(f"{path}:{row_number} repeats species {species!r}")
        fasta = normalize_data_path(row.get(fasta_col, ""), base_dir, keep_relative_paths)
        if not fasta:
            fail(f"{path}:{row_number} has no FASTA path for species {species!r}")
        latin = row.get(latin_col, "").strip() if latin_col else species
        gff = normalize_data_path(row.get(gff_col, ""), base_dir, keep_relative_paths) if gff_col else ""
        metadata[species] = SpeciesRecord(
            species=species,
            latin_name=latin or species,
            fasta_filename=fasta,
            gff_filename=gff,
        )
    return metadata


def read_datasets(path: Path) -> List[DatasetRecord]:
    fieldnames, rows = read_table(path)
    dataset_col = column_name(fieldnames, ("dataset_id", "branch_id", "branch", "node_id"), True)
    newick_col = column_name(fieldnames, ("newick_tree", "newick", "subtree_newick"), True)
    focal_col = column_name(fieldnames, ("focal_species", "focals", "focal_species_list"), False)
    target_col = column_name(fieldnames, ("target_species", "targets", "comparison_species"), False)

    datasets: List[DatasetRecord] = []
    assert dataset_col is not None
    assert newick_col is not None
    for row_number, row in enumerate(rows, start=2):
        dataset_id = row.get(dataset_col, "").strip()
        newick_tree = row.get(newick_col, "").strip()
        if not dataset_id:
            fail(f"{path}:{row_number} has an empty dataset_id")
        if not newick_tree:
            fail(f"{path}:{row_number} has an empty Newick tree")
        leaves = tuple(parse_newick_leaves(newick_tree))
        if not leaves:
            fail(f"{path}:{row_number} Newick tree has no detectable leaves")

        focal_species = tuple(split_species_list(row.get(focal_col, ""))) if focal_col else leaves
        if not focal_species:
            focal_species = leaves

        target_species = tuple(split_species_list(row.get(target_col, ""))) if target_col else tuple()
        datasets.append(
            DatasetRecord(
                dataset_id=dataset_id,
                newick_tree=newick_tree,
                leaves=leaves,
                focal_species=focal_species,
                target_species=target_species,
            )
        )
    return datasets


def split_species_list(value: str) -> List[str]:
    if not value:
        return []
    return [item for item in re.split(r"[,;\s]+", value.strip()) if item]


def parse_newick_leaves(newick: str) -> List[str]:
    leaves: List[str] = []
    seen: Set[str] = set()
    previous_token = ""
    i = 0

    while i < len(newick):
        char = newick[i]
        if char.isspace():
            i += 1
            continue
        if char == ":":
            i += 1
            while i < len(newick) and newick[i] not in ",);":
                i += 1
            continue
        if char in "(),;":
            previous_token = char
            i += 1
            continue
        if char == "'":
            i += 1
            label_chars: List[str] = []
            while i < len(newick):
                if newick[i] == "'":
                    i += 1
                    break
                label_chars.append(newick[i])
                i += 1
            label = "".join(label_chars).strip()
        else:
            start = i
            while i < len(newick) and newick[i] not in "(),:; \t\r\n":
                i += 1
            label = newick[start:i].strip()

        if label and previous_token in ("(", ",") and label not in seen:
            leaves.append(label)
            seen.add(label)
        previous_token = "label"

    return leaves


def ensure_species_known(species: Iterable[str], metadata: Dict[str, SpeciesRecord], context: str) -> None:
    missing = sorted(set(species) - set(metadata))
    if missing:
        fail(f"{context} references species missing from metadata: {', '.join(missing)}")


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9.+-]+", "_", value.strip())
    return cleaned.strip("_") or "unnamed"


def canonical_pair(species_1: str, species_2: str) -> Tuple[str, str]:
    if species_1.lower() <= species_2.lower():
        return species_1, species_2
    return species_2, species_1


def write_ksrates_config(
    config_path: Path,
    dataset: DatasetRecord,
    focal_species: str,
    metadata: Dict[str, SpeciesRecord],
    args: argparse.Namespace,
    peak_database_path: Path,
    ks_list_database_path: Path,
    paralog_database_path: Optional[Path] = None,
) -> None:
    records = [metadata[species] for species in dataset.leaves]
    focal_record = metadata[focal_species]

    lines = [
        "[SPECIES]",
        f"focal_species = {focal_species}",
        f"newick_tree = {dataset.newick_tree}",
        "",
        format_mapping("latin_names", [(record.species, record.latin_name) for record in records]),
        "",
        format_mapping("fasta_filenames", [(record.species, record.fasta_filename) for record in records]),
        "",
    ]
    if focal_record.gff_filename:
        lines.append(f"gff_filename = {focal_record.gff_filename}")
    lines.extend(
        [
            f"peak_database_path = {peak_database_path}",
            f"ks_list_database_path = {ks_list_database_path}",
        ]
    )
    if paralog_database_path is not None:
        lines.append(f"ks_list_paralog_database_path = {paralog_database_path}")
    lines.extend(
        [
            "",
            "",
            "[ANALYSIS SETTING]",
            f"paranome = {args.paranome}",
            f"collinearity = {args.collinearity}",
            f"reciprocal_retention = {args.reciprocal_retention}",
            "",
            f"gff_feature = {args.gff_feature}",
            f"gff_attribute = {args.gff_attribute}",
            "",
            f"max_number_outgroups = {args.max_outgroups}",
            f"consensus_mode_for_multiple_outgroups = {args.consensus_mode}",
            "",
            "",
            "[PARAMETERS]",
            "x_axis_max_limit_paralogs_plot = 5",
            "bin_width_paralogs = 0.1",
            "y_axis_max_limit_paralogs_plot = None",
            "",
            f"num_bootstrap_iterations = {args.num_bootstrap_iterations}",
            (
                "divergence_colors = Red, MediumBlue, Goldenrod, Crimson, "
                "ForestGreen, Gray, SaddleBrown, Black"
            ),
            "x_axis_max_limit_orthologs_plots = 5",
            "bin_width_orthologs = 0.1",
            "",
            f"max_ks_paralogs = {args.max_ks_paralogs}",
            f"max_ks_orthologs = {args.max_ks_orthologs}",
            "",
        ]
    )

    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text("\n".join(lines))


def format_mapping(name: str, values: Sequence[Tuple[str, str]]) -> str:
    if not values:
        return f"{name} ="
    key_width = max(len(key) for key, _ in values)
    prefix = f"{name} = "
    continuation = " " * len(prefix)
    lines = []
    for index, (key, value) in enumerate(values):
        trailer = "," if index + 1 < len(values) else ""
        left = prefix if index == 0 else continuation
        lines.append(f"{left}{key:<{key_width}} : {value}{trailer}")
    return "\n".join(lines)


def make_ksrates_command(args: argparse.Namespace, command_args: Sequence[str]) -> List[str]:
    return shlex.split(args.ksrates_command) + list(command_args)


def shell_join(parts: Sequence[str]) -> str:
    return " ".join(shlex.quote(part) for part in parts)


def build_init_command(args: argparse.Namespace, config_path: Path) -> List[str]:
    command = make_ksrates_command(args, ("init", str(config_path)))
    if args.expert_config is not None:
        command.extend(("--expert", str(args.expert_config)))
    return command


def build_ortholog_command(
    args: argparse.Namespace,
    config_path: Path,
    species_1: str,
    species_2: str,
) -> List[str]:
    command = make_ksrates_command(args, ("orthologs-ks", str(config_path)))
    if args.expert_config is not None:
        command.extend(("--expert", str(args.expert_config)))
    command.extend((species_1, species_2, "--n-threads", str(args.ortholog_threads)))
    return command


def build_ortholog_analysis_command(
    args: argparse.Namespace,
    config_path: Path,
    pair_tsv_path: Path,
) -> List[str]:
    """
    orthologs-analysis is what actually estimates the Ks distribution peak and writes it (plus
    the raw Ks list) into the shared peak_database_path/ks_list_database_path TSVs -
    orthologs-ks on its own only produces raw local Ks values, never touching those databases.
    Pointed at a single-pair TSV (pair_tsv_path) rather than a full per-focal ortholog_pairs file,
    so this only ever processes the one pair orthologs-ks was just run for.
    """
    command = make_ksrates_command(args, ("orthologs-analysis", str(config_path)))
    if args.expert_config is not None:
        command.extend(("--expert", str(args.expert_config)))
    command.extend(("--ortholog-pairs", str(pair_tsv_path)))
    return command


def build_paralog_command(
    args: argparse.Namespace,
    config_path: Path,
) -> List[str]:
    command = make_ksrates_command(args, ("paralogs-ks", str(config_path)))
    if args.expert_config is not None:
        command.extend(("--expert", str(args.expert_config)))
    command.extend(("--n-threads", str(args.paralog_threads)))
    return command


def build_check_paralog_db_command(
    args: argparse.Namespace,
    config_path: Path,
    analysis_type: str,
) -> List[str]:
    command = make_ksrates_command(args, ("check-paralog-db", str(config_path)))
    if args.expert_config is not None:
        command.extend(("--expert", str(args.expert_config)))
    command.extend(("--type", analysis_type))
    return command


def species_still_needs_paralog_ks(args: argparse.Namespace, config_path: Path) -> bool:
    """
    Checks the paralog Ks database (via 'ksrates check-paralog-db') for each analysis type
    requested by --paranome/--collinearity/--reciprocal-retention. This is what makes rerunning
    this script after changing those flags (e.g. turning collinearity on for an already-populated
    tree) cheap: without it, paralogs_ks_commands.sh would always list every unique species again,
    submitting a job for species that already have everything they need.

    :return: True if at least one requested analysis type is still missing for this species
             (a paralogs-ks command should be emitted for it), False only if every requested type
             is already present (nothing to do for this species).
    """
    requested = []
    if args.paranome == "yes":
        requested.append("paranome")
    if args.collinearity == "yes":
        requested.append("anchors")
    if args.reciprocal_retention == "yes":
        requested.append("recret")
    if not requested:
        return False

    for analysis_type in requested:
        command = build_check_paralog_db_command(args, config_path, analysis_type)
        try:
            # Any failure to run the check (database unreachable, ksrates missing, ...) is
            # treated the same as "missing": safe to over-include a species here (redundant
            # work, caught by ks_paralogs()'s own local-disk resume logic) rather than silently
            # drop one that actually still needs something.
            result = subprocess.run(command, capture_output=True, text=True)
        except OSError:
            return True
        if result.returncode != 0:
            return True
    return False


def run_init(args: argparse.Namespace, config: ConfigRecord) -> None:
    config.init_run_dir.mkdir(parents=True, exist_ok=True)
    command = build_init_command(args, config.config_path)
    print(
        f"[init] dataset={config.dataset_id} focal={config.focal_species}: "
        f"{shell_join(command)}",
        flush=True,
    )
    subprocess.run(command, cwd=str(config.init_run_dir), check=True)


def find_ksrates_file(init_run_dir: Path, focal_species: str, stem: str) -> Optional[Path]:
    rate_dir = init_run_dir / "rate_adjustment" / focal_species
    candidates = [
        rate_dir / f"{stem}_{focal_species}.tsv",
        rate_dir / f"{stem}_{focal_species}.txt",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    matches = sorted(rate_dir.glob(f"{stem}_*"))
    return matches[0] if matches else None


def read_ortholog_pairs(path: Path) -> List[Tuple[str, str]]:
    pairs: List[Tuple[str, str]] = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = re.split(r"[\t, ]+", stripped)
            if len(fields) < 2:
                continue
            first, second = fields[0], fields[1]
            if normalize_column_name(first) in {"species1", "speciesa", "focalspecies"}:
                continue
            if normalize_column_name(second) in {"species2", "speciesb", "sisterspecies"}:
                continue
            if first == second:
                print(
                    f"WARNING: skipping self-pair in {path}:{line_number}: {stripped}",
                    file=sys.stderr,
                )
                continue
            pairs.append(canonical_pair(first, second))
    return pairs


def derive_pairs_from_target_trios(
    path: Path,
    focal_species: str,
    target_species: Sequence[str],
) -> List[Tuple[str, str]]:
    targets = set(target_species)
    pairs: Set[Tuple[str, str]] = set()

    with path.open(newline="") as handle:
        sample = handle.readline()
        if not sample:
            return []
        delimiter = "\t" if "\t" in sample else ","
        handle.seek(0)
        reader = csv.reader(handle, delimiter=delimiter)
        rows = list(reader)

    if not rows:
        return []

    header = [normalize_column_name(value) for value in rows[0]]
    has_header = any(name in header for name in ("focalspecies", "sisterspecies", "outspecies"))
    data_rows = rows[1:] if has_header else rows

    if has_header:
        focal_idx = header.index("focalspecies") if "focalspecies" in header else None
        sister_idx = header.index("sisterspecies") if "sisterspecies" in header else None
        out_idx = header.index("outspecies") if "outspecies" in header else None
        if focal_idx is None or sister_idx is None or out_idx is None:
            fail(f"Cannot identify focal/sister/outgroup columns in {path}")

        for row in data_rows:
            if max(focal_idx, sister_idx, out_idx) >= len(row):
                continue
            focal = row[focal_idx].strip()
            sister = row[sister_idx].strip()
            outgroup = row[out_idx].strip()
            add_target_trio_pairs(pairs, focal_species, targets, focal, sister, outgroup)
    else:
        for row in data_rows:
            fields = [value.strip() for value in row if value.strip()]
            if len(fields) >= 4 and fields[1] == focal_species:
                focal, sister, outgroup = fields[1], fields[2], fields[3]
            elif len(fields) >= 3 and fields[0] == focal_species:
                focal, sister, outgroup = fields[0], fields[1], fields[2]
            else:
                continue
            add_target_trio_pairs(pairs, focal_species, targets, focal, sister, outgroup)

    return sorted(pairs, key=lambda pair: (pair[0].lower(), pair[1].lower()))


def add_target_trio_pairs(
    pairs: Set[Tuple[str, str]],
    expected_focal: str,
    target_species: Set[str],
    focal: str,
    sister: str,
    outgroup: str,
) -> None:
    if focal != expected_focal:
        return
    if target_species and sister not in target_species:
        return
    if not sister or not outgroup:
        return
    pairs.add(canonical_pair(focal, sister))
    pairs.add(canonical_pair(focal, outgroup))
    pairs.add(canonical_pair(sister, outgroup))


def collect_pairs_for_config(args: argparse.Namespace, config: ConfigRecord) -> List[PairOccurrence]:
    if args.pair_scope == "init-pairs":
        source_file = find_ksrates_file(config.init_run_dir, config.focal_species, "ortholog_pairs")
        source_kind = "ortholog_pairs"
        if source_file is None:
            return []
        pairs = read_ortholog_pairs(source_file)
    else:
        source_file = find_ksrates_file(config.init_run_dir, config.focal_species, "ortholog_trios")
        source_kind = "target_ortholog_trios"
        if source_file is None:
            return []
        pairs = derive_pairs_from_target_trios(
            source_file,
            config.focal_species,
            config.target_species,
        )

    return [
        PairOccurrence(
            species_1=species_1,
            species_2=species_2,
            dataset_id=config.dataset_id,
            focal_species=config.focal_species,
            config_path=config.config_path,
            init_run_dir=config.init_run_dir,
            source_file=source_file,
            source_kind=source_kind,
        )
        for species_1, species_2 in pairs
    ]


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow([str(value) for value in row])


def write_command_file(path: Path, entries: Sequence[Tuple[Path, Sequence[Sequence[str]]]]) -> None:
    """
    Write a bash script with one line per entry: "( cd <cwd> && <cmd1> && <cmd2> && ... )".
    Each entry's command list lets several ksrates invocations that must run in sequence (e.g.
    orthologs-ks then orthologs-analysis for the same pair) live on one script line, sharing one
    cwd, with "&&" ensuring a later command only runs if the earlier one succeeded.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("#!/usr/bin/env bash\n")
        handle.write("set -euo pipefail\n\n")
        for cwd, commands in entries:
            joined = " && ".join(shell_join(command) for command in commands)
            handle.write(f"( cd {shlex.quote(str(cwd))} && {joined} )\n")
    path.chmod(0o755)


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


def global_pair_rows(
    pair_to_occurrences: Dict[Tuple[str, str], List[PairOccurrence]]
) -> Iterable[Tuple[object, ...]]:
    for (species_1, species_2), occurrences in sorted(
        pair_to_occurrences.items(),
        key=lambda item: (item[0][0].lower(), item[0][1].lower()),
    ):
        datasets = sorted({occurrence.dataset_id for occurrence in occurrences})
        focals = sorted({occurrence.focal_species for occurrence in occurrences}, key=str.lower)
        chosen = sorted(
            occurrences,
            key=lambda item: (item.dataset_id, item.focal_species.lower(), str(item.config_path)),
        )[0]
        yield (
            species_1,
            species_2,
            ",".join(datasets),
            ",".join(focals),
            chosen.config_path,
            chosen.init_run_dir,
        )


def dataset_pair_rows(
    dataset_pair_map: Dict[Tuple[str, str], List[PairOccurrence]]
) -> Iterable[Tuple[object, ...]]:
    for (species_1, species_2), occurrences in sorted(
        dataset_pair_map.items(),
        key=lambda item: (item[0][0].lower(), item[0][1].lower()),
    ):
        dataset_id = sorted({occurrence.dataset_id for occurrence in occurrences})[0]
        focals = sorted({occurrence.focal_species for occurrence in occurrences}, key=str.lower)
        chosen = sorted(
            occurrences,
            key=lambda item: (item.focal_species.lower(), str(item.config_path)),
        )[0]
        yield (
            species_1,
            species_2,
            dataset_id,
            ",".join(focals),
            chosen.config_path,
        )


if __name__ == "__main__":
    main()
