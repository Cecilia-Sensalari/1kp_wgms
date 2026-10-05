"""Reading species_metadata.tsv and the dataset-list TSV into typed records."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

from ksrates_batch_common import (
    DatasetRecord,
    SpeciesRecord,
    column_name,
    fail,
    normalize_data_path,
    read_table,
)


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
