#!/usr/bin/env python3
"""Prepare branch-centered ksrates configs and ortholog-pair workloads.

The script expects two small TSV files:

1. A species metadata table with at least:
      species    fasta_filename
   Optional columns:
      latin_name    gff_filename

2. A branch/subtree table with at least:
      branch_id    newick_tree
   Optional columns:
      focal_species    target_species

By default, every leaf in a subtree becomes a focal species. If the
``focal_species`` column is present, only those comma/semicolon/space-separated
species are used as focal species for that subtree.

The normal workflow is:

    python prepare_branch_ksrates_configs.py \
        --subtrees branch_subtrees.tsv \
        --species-metadata species_metadata.tsv \
        --out-dir ksrates_branch_setup \
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
class SubtreeRecord:
    branch_id: str
    newick_tree: str
    leaves: Tuple[str, ...]
    focal_species: Tuple[str, ...]
    target_species: Tuple[str, ...]


@dataclass(frozen=True)
class ConfigRecord:
    branch_id: str
    focal_species: str
    config_path: Path
    init_run_dir: Path
    target_species: Tuple[str, ...]


@dataclass(frozen=True)
class PairOccurrence:
    species_1: str
    species_2: str
    branch_id: str
    focal_species: str
    config_path: Path
    init_run_dir: Path
    source_file: Path
    source_kind: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate one ksrates config per selected focal species in each "
            "branch-centered subtree, run/plan ksrates init, and collect the "
            "ortholog Ks species pairs reported by ksrates."
        )
    )
    parser.add_argument(
        "--subtrees",
        required=True,
        type=Path,
        help="TSV with branch_id and newick_tree columns.",
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
            "unique species across all branches (mirrors the existing ortholog-pair dedup). The "
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
            "ortholog_trios_<focal> rows whose sister/comparison species is in "
            "the subtree table's target_species column."
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


def read_subtrees(path: Path) -> List[SubtreeRecord]:
    fieldnames, rows = read_table(path)
    branch_col = column_name(fieldnames, ("branch_id", "branch", "node_id", "test_branch"), True)
    newick_col = column_name(fieldnames, ("newick_tree", "newick", "subtree_newick"), True)
    focal_col = column_name(fieldnames, ("focal_species", "focals", "focal_species_list"), False)
    target_col = column_name(fieldnames, ("target_species", "targets", "comparison_species"), False)

    subtrees: List[SubtreeRecord] = []
    assert branch_col is not None
    assert newick_col is not None
    for row_number, row in enumerate(rows, start=2):
        branch_id = row.get(branch_col, "").strip()
        newick_tree = row.get(newick_col, "").strip()
        if not branch_id:
            fail(f"{path}:{row_number} has an empty branch_id")
        if not newick_tree:
            fail(f"{path}:{row_number} has an empty Newick tree")
        leaves = tuple(parse_newick_leaves(newick_tree))
        if not leaves:
            fail(f"{path}:{row_number} Newick tree has no detectable leaves")

        focal_species = tuple(split_species_list(row.get(focal_col, ""))) if focal_col else leaves
        if not focal_species:
            focal_species = leaves

        target_species = tuple(split_species_list(row.get(target_col, ""))) if target_col else tuple()
        subtrees.append(
            SubtreeRecord(
                branch_id=branch_id,
                newick_tree=newick_tree,
                leaves=leaves,
                focal_species=focal_species,
                target_species=target_species,
            )
        )
    return subtrees


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
    subtree: SubtreeRecord,
    focal_species: str,
    metadata: Dict[str, SpeciesRecord],
    args: argparse.Namespace,
    peak_database_path: Path,
    ks_list_database_path: Path,
    paralog_database_path: Optional[Path] = None,
) -> None:
    records = [metadata[species] for species in subtree.leaves]
    focal_record = metadata[focal_species]

    lines = [
        "[SPECIES]",
        f"focal_species = {focal_species}",
        f"newick_tree = {subtree.newick_tree}",
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


def build_paralog_command(
    args: argparse.Namespace,
    config_path: Path,
) -> List[str]:
    command = make_ksrates_command(args, ("paralogs-ks", str(config_path)))
    if args.expert_config is not None:
        command.extend(("--expert", str(args.expert_config)))
    command.extend(("--n-threads", str(args.paralog_threads)))
    return command


def run_init(args: argparse.Namespace, config: ConfigRecord) -> None:
    config.init_run_dir.mkdir(parents=True, exist_ok=True)
    command = build_init_command(args, config.config_path)
    print(
        f"[init] branch={config.branch_id} focal={config.focal_species}: "
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
            branch_id=config.branch_id,
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


def write_command_file(path: Path, commands: Sequence[Tuple[Path, Sequence[str]]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("#!/usr/bin/env bash\n")
        handle.write("set -euo pipefail\n\n")
        for cwd, command in commands:
            handle.write(f"( cd {shlex.quote(str(cwd))} && {shell_join(command)} )\n")
    path.chmod(0o755)


def main() -> None:
    args = parse_args()
    args.subtrees = args.subtrees.resolve()
    args.species_metadata = args.species_metadata.resolve()
    args.out_dir = args.out_dir.resolve()
    if args.expert_config is not None:
        args.expert_config = args.expert_config.resolve()
    if args.paralog_database is not None:
        args.paralog_database = args.paralog_database.resolve()

    metadata = read_species_metadata(args.species_metadata, args.keep_relative_paths)
    subtrees = read_subtrees(args.subtrees)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    config_root = args.out_dir / "configs"
    init_root = args.out_dir / "init_runs"
    branch_root = args.out_dir / "branches"
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
    for subtree in subtrees:
        ensure_species_known(subtree.leaves, metadata, f"branch {subtree.branch_id!r} Newick")
        ensure_species_known(subtree.focal_species, metadata, f"branch {subtree.branch_id!r} focal list")
        unknown_focals = sorted(set(subtree.focal_species) - set(subtree.leaves))
        if unknown_focals:
            fail(
                f"branch {subtree.branch_id!r} focal species are not in the subtree: "
                f"{', '.join(unknown_focals)}"
            )
        if subtree.target_species:
            ensure_species_known(
                subtree.target_species,
                metadata,
                f"branch {subtree.branch_id!r} target list",
            )
            unknown_targets = sorted(set(subtree.target_species) - set(subtree.leaves))
            if unknown_targets:
                fail(
                    f"branch {subtree.branch_id!r} target species are not in the subtree: "
                    f"{', '.join(unknown_targets)}"
                )

        branch_dir_name = safe_name(subtree.branch_id)
        for focal_species in subtree.focal_species:
            config_path = (
                config_root
                / branch_dir_name
                / f"config_{branch_dir_name}_{safe_name(focal_species)}.txt"
            )
            init_run_dir = init_root / branch_dir_name / safe_name(focal_species)
            write_ksrates_config(
                config_path,
                subtree,
                focal_species,
                metadata,
                args,
                peak_database_path,
                ks_list_database_path,
                paralog_database_path=args.paralog_database,
            )
            configs.append(
                ConfigRecord(
                    branch_id=subtree.branch_id,
                    focal_species=focal_species,
                    config_path=config_path,
                    init_run_dir=init_run_dir,
                    target_species=subtree.target_species,
                )
            )

    if args.run_init:
        for config in configs:
            run_init(args, config)

    init_commands = [(config.init_run_dir, build_init_command(args, config.config_path)) for config in configs]
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
            "branch_id",
            "focal_species",
            "config_path",
            "init_run_dir",
            "target_species",
        ),
        (
            (
                config.branch_id,
                config.focal_species,
                config.config_path,
                config.init_run_dir,
                ",".join(config.target_species),
            )
            for config in configs
        ),
    )

    pair_to_occurrences: Dict[Tuple[str, str], List[PairOccurrence]] = defaultdict(list)
    branch_to_occurrences: Dict[str, List[PairOccurrence]] = defaultdict(list)
    for occurrence in all_occurrences:
        pair = canonical_pair(occurrence.species_1, occurrence.species_2)
        pair_to_occurrences[pair].append(occurrence)
        branch_to_occurrences[occurrence.branch_id].append(occurrence)

    write_tsv(
        args.out_dir / "ortholog_pairs_by_source.tsv",
        (
            "species_1",
            "species_2",
            "branch_id",
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
                occurrence.branch_id,
                occurrence.focal_species,
                occurrence.config_path,
                occurrence.init_run_dir,
                occurrence.source_kind,
                occurrence.source_file,
            )
            for occurrence in sorted(
                all_occurrences,
                key=lambda item: (
                    item.branch_id,
                    item.focal_species.lower(),
                    item.species_1.lower(),
                    item.species_2.lower(),
                ),
            )
        ),
    )

    write_tsv(
        args.out_dir / "ortholog_pairs.tsv",
        (
            "species_1",
            "species_2",
            "branches",
            "focal_species",
            "chosen_config",
            "chosen_init_run_dir",
        ),
        global_pair_rows(pair_to_occurrences),
    )

    for branch_id, occurrences in branch_to_occurrences.items():
        branch_pair_map: Dict[Tuple[str, str], List[PairOccurrence]] = defaultdict(list)
        for occurrence in occurrences:
            branch_pair_map[canonical_pair(occurrence.species_1, occurrence.species_2)].append(occurrence)
        write_tsv(
            branch_root / safe_name(branch_id) / "ortholog_pairs.tsv",
            (
                "species_1",
                "species_2",
                "branch_id",
                "focal_species",
                "chosen_config",
            ),
            branch_pair_rows(branch_pair_map),
        )

    ortholog_commands = []
    for (species_1, species_2), occurrences in sorted(
        pair_to_occurrences.items(),
        key=lambda item: (item[0][0].lower(), item[0][1].lower()),
    ):
        chosen = sorted(
            occurrences,
            key=lambda item: (item.branch_id, item.focal_species.lower(), str(item.config_path)),
        )[0]
        ortholog_commands.append(
            (
                ortholog_work_dir,
                build_ortholog_command(args, chosen.config_path, species_1, species_2),
            )
        )
    write_command_file(args.out_dir / "orthologs_ks_commands.sh", ortholog_commands)

    # Species-level dedup for paralog population, mirroring the pair-level dedup above: the same
    # real species can appear as focal_species in many different branches' subtrees (since
    # species is drawn from one global species_metadata.tsv, the string itself is already a
    # stable cross-branch identity - no latin-name translation needed, unlike per-branch informal
    # names elsewhere). Without this, feeding every (branch, focal) config to paralogs-ks-multi
    # would launch one paralogs-ks call per branch a species appears in, not once per species -
    # exactly the duplicate-computation race this preparation step exists to avoid.
    if args.paralog_database is not None:
        species_to_occurrences: Dict[str, List[ConfigRecord]] = defaultdict(list)
        for config in configs:
            species_to_occurrences[config.focal_species].append(config)

        paralog_commands = []
        for species in sorted(species_to_occurrences, key=str.lower):
            occurrences = species_to_occurrences[species]
            chosen = sorted(
                occurrences,
                key=lambda item: (item.branch_id, str(item.config_path)),
            )[0]
            paralog_commands.append((paralog_work_dir, build_paralog_command(args, chosen.config_path)))
        write_command_file(args.out_dir / "paralogs_ks_commands.sh", paralog_commands)

        write_tsv(
            args.out_dir / "paralog_species.tsv",
            ("species", "branches", "chosen_config"),
            (
                (
                    species,
                    ",".join(sorted({occ.branch_id for occ in occurrences}, key=str.lower)),
                    sorted(occurrences, key=lambda item: (item.branch_id, str(item.config_path)))[0].config_path,
                )
                for species, occurrences in sorted(species_to_occurrences.items(), key=lambda item: item[0].lower())
            ),
        )

        print(f"Collected {len(species_to_occurrences)} unique paralog species")
        print(f"Wrote paralog command file: {args.out_dir / 'paralogs_ks_commands.sh'}")

    print(f"Wrote {len(configs)} ksrates configs to {config_root}")
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
        branches = sorted({occurrence.branch_id for occurrence in occurrences})
        focals = sorted({occurrence.focal_species for occurrence in occurrences}, key=str.lower)
        chosen = sorted(
            occurrences,
            key=lambda item: (item.branch_id, item.focal_species.lower(), str(item.config_path)),
        )[0]
        yield (
            species_1,
            species_2,
            ",".join(branches),
            ",".join(focals),
            chosen.config_path,
            chosen.init_run_dir,
        )


def branch_pair_rows(
    branch_pair_map: Dict[Tuple[str, str], List[PairOccurrence]]
) -> Iterable[Tuple[object, ...]]:
    for (species_1, species_2), occurrences in sorted(
        branch_pair_map.items(),
        key=lambda item: (item[0][0].lower(), item[0][1].lower()),
    ):
        branch_id = sorted({occurrence.branch_id for occurrence in occurrences})[0]
        focals = sorted({occurrence.focal_species for occurrence in occurrences}, key=str.lower)
        chosen = sorted(
            occurrences,
            key=lambda item: (item.focal_species.lower(), str(item.config_path)),
        )[0]
        yield (
            species_1,
            species_2,
            branch_id,
            ",".join(focals),
            chosen.config_path,
        )


if __name__ == "__main__":
    main()
