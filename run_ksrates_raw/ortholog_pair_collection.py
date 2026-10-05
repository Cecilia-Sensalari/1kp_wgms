"""Collecting, deduplicating, and building commands for ortholog species-pairs."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from ksrates_batch_common import (
    ConfigRecord,
    PairOccurrence,
    canonical_pair,
    fail,
    make_ksrates_command,
    normalize_column_name,
)


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
