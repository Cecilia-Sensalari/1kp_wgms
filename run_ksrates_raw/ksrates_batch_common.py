"""Shared dataclasses and generic helpers used across the batch_ksrates preparation scripts.

Nothing in this module is specific to datasets, ortholog pairs, or paralog species - just the
record types and generic TSV/path/command-string utilities the other modules build on.
"""

from __future__ import annotations

import argparse
import csv
import re
import shlex
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


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


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9.+-]+", "_", value.strip())
    return cleaned.strip("_") or "unnamed"


def canonical_pair(species_1: str, species_2: str) -> Tuple[str, str]:
    if species_1.lower() <= species_2.lower():
        return species_1, species_2
    return species_2, species_1


def make_ksrates_command(args: argparse.Namespace, command_args: Sequence[str]) -> List[str]:
    return shlex.split(args.ksrates_command) + list(command_args)


def shell_join(parts: Sequence[str]) -> str:
    return " ".join(shlex.quote(part) for part in parts)


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
