"""Shared dataclasses and generic helpers used across the batch_ksrates preparation scripts.

Nothing in this module is specific to datasets, ortholog pairs, or paralog species - just the
record types and generic TSV/path/command-string utilities the other modules build on.
"""

from __future__ import annotations

import argparse
import csv
import re
import shlex
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple, TypeVar


@dataclass(frozen=True)
class SpeciesRecord:
    species: str
    latin_name: str
    fasta_filename: str
    gff_filename: str = ""
    gff_feature: str = ""
    gff_attribute: str = ""


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


def ksrates_species_name(species_id: str) -> str:
    """Informal species name for ksrates, which allows no underscore or space in a name.

    A genome species id such as Phavu_v1.0 becomes Phavu (the version stays in its latin name); the 4-letter
    transcriptome ids are unchanged. Used for every species name read from the metadata and dataset tables, so
    that configs, pair files and commands all use the same names.
    """
    return species_id.split("_")[0]


def ksrates_newick(newick: str) -> str:
    """The Newick tree with every species name converted by ksrates_species_name."""
    return re.sub(r"[^(),;:]+", lambda match: ksrates_species_name(match.group(0)), newick)


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9.+-]+", "_", value.strip())
    return cleaned.strip("_") or "unnamed"


def canonical_pair(species_1: str, species_2: str) -> Tuple[str, str]:
    if species_1.lower() <= species_2.lower():
        return species_1, species_2
    return species_2, species_1


def make_ksrates_command(args: argparse.Namespace, command_args: Sequence[str]) -> List[str]:
    return shlex.split(args.ksrates_command) + list(command_args)


def ks_database_has_data(command: Sequence[str]) -> bool:
    """
    Runs one 'ksrates check-paralog-db' / 'check-ortholog-db' command and interprets its exit
    status: 0 = the data is in the Ks database (True), 1 = it is not, or the database is disabled
    in the expert config (False). Anything else stops the whole run, in particular 2 = the Ks
    database is requested but unreachable: treating "cannot tell" as "missing" would queue every
    species/pair for recomputation while the database server is down.
    """
    try:
        result = subprocess.run(list(command), capture_output=True, text=True)
    except OSError as error:
        fail(f"could not run [{' '.join(command)}]: {error}")
    if result.returncode == 0:
        return True
    if result.returncode == 1:
        return False
    fail(
        f"[{' '.join(command)}] exited with status {result.returncode} (2 = the Ks database is "
        f"unreachable): {result.stderr.strip() or result.stdout.strip()}. Is the Ks database "
        f"server running, and does the address file in the config point to it?"
    )
    return False  # unreachable: fail() always raises


T = TypeVar("T")
R = TypeVar("R")


def map_checks_concurrently(function: Callable[[T], R], items: Sequence[T], workers: int) -> List[R]:
    """
    Runs function on every item in parallel threads, keeping the input order. The checks are small
    network round trips to the sqld server, not CPU-bound, so they run concurrently - run serially
    they would dominate wall-clock time when upscaling to thousands of species. If one raises (e.g. database
    unreachable) the checks still queued are cancelled, since each would fail the same way.
    """
    with ThreadPoolExecutor(max_workers=workers) as executor:
        try:
            return list(executor.map(function, items))
        except BaseException:
            executor.shutdown(wait=False, cancel_futures=True)
            raise


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
