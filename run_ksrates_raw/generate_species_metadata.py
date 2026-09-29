#!/usr/bin/env python3
"""Generate species_metadata.tsv for prepare_branch_ksrates_configs.py.

Joins the 1KP species table (sample ID -> latin name) against the CNGB transcriptome
FASTA files actually present on disk (sample IDs no longer come from CyVerse - see
run_paralog_batches.py, whose CyVerse-specific download/path logic is now unused for
this project since transcriptomes were instead already downloaded from CNGB).

Only species with BOTH a species-table entry AND an on-disk FASTA file are included -
species missing either are reported, not silently dropped.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


SPECIES_TABLE = Path(
    "/group/esb/cesen/1kp/source_data/1.species_dataset/1kp_paper_2019_suptab1_species.tsv"
)
TRANSCRIPTOMES_DIR = Path(
    "/group/esb/cesen/1kp/source_data/2.transcriptomes/1kp/cngb_transcriptomes"
)
FASTA_SUFFIX = "-translated-nucleotides.fa"
OUTPUT_PATH = Path("/group/esb/cesen/1kp/ks_analysis/1kp/ksrates_raw/species_metadata.tsv")


def read_species_table(path: Path) -> dict[str, str]:
    """Map 1KP Index ID -> latin species name."""
    mapping: dict[str, str] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample_id = row["1KP Index ID"].strip()
            species = row["Species"].strip()
            if sample_id and species:
                mapping[sample_id] = species
    return mapping


def find_transcriptome_fastas(directory: Path) -> dict[str, Path]:
    """Map sample ID -> FASTA path, for every '<ID>-translated-nucleotides.fa' file present."""
    fastas: dict[str, Path] = {}
    for fasta_path in directory.glob(f"*{FASTA_SUFFIX}"):
        sample_id = fasta_path.name[: -len(FASTA_SUFFIX)]
        fastas[sample_id] = fasta_path
    return fastas


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species-table", type=Path, default=SPECIES_TABLE)
    parser.add_argument("--transcriptomes-dir", type=Path, default=TRANSCRIPTOMES_DIR)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    species_by_id = read_species_table(args.species_table)
    fastas_by_id = find_transcriptome_fastas(args.transcriptomes_dir)

    known_ids = set(species_by_id)
    available_ids = set(fastas_by_id)
    usable_ids = sorted(known_ids & available_ids)
    missing_fasta = sorted(known_ids - available_ids)
    missing_table_entry = sorted(available_ids - known_ids)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["species", "fasta_filename", "latin_name"])
        for sample_id in usable_ids:
            # Append the 4-letter sample ID to the latin name (e.g. "Genus species AAXJ") so it
            # stays unique and traceable back to its source sample even when two different
            # samples share the same latin species name.
            latin_name = f"{species_by_id[sample_id]} {sample_id}"
            writer.writerow([sample_id, str(fastas_by_id[sample_id]), latin_name])

    print(f"Wrote {len(usable_ids)} species to {args.output}")
    if missing_fasta:
        print(f"WARNING: {len(missing_fasta)} species table entries have no CNGB FASTA file:")
        print(f"  {', '.join(missing_fasta)}")
    if missing_table_entry:
        print(f"WARNING: {len(missing_table_entry)} CNGB FASTA files have no species table entry:")
        print(f"  {', '.join(missing_table_entry)}")


if __name__ == "__main__":
    main()
