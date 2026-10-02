#!/usr/bin/env python3
"""Generate species_metadata.tsv for prepare_branch_ksrates_configs.py.

Joins the 1KP species table (sample ID -> latin name) against the CNGB transcriptome
FASTA files actually present on disk (sample IDs no longer come from CyVerse - see
run_paralog_batches.py, whose CyVerse-specific download/path logic is now unused for
this project since transcriptomes were instead already downloaded from CNGB).

Species with no transcriptome are filled in from the genome manifest
(genome_species_manifest.tsv, under source_data/2.genomes/1kp/) where available - a
species with both a transcriptome and a genome manifest entry always prefers the
transcriptome. Species with neither are reported, not silently dropped.
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
GENOME_MANIFEST = Path(
    "/group/esb/cesen/1kp/source_data/2.genomes/1kp/genome_species_manifest.tsv"
)
FASTA_SUFFIX = "-translated-nucleotides.fa"
OUTPUT_PATH = Path("/group/esb/cesen/1kp/ks_analysis/1kp/species_metadata.tsv")


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


def read_genome_manifest(path: Path) -> dict[str, tuple[str, str, str]]:
    """Map species_id -> (latin_name, cds_filepath, gff_filepath) from genome_species_manifest.tsv.

    Entries are included as-is even if cds_filepath doesn't exist on disk yet (e.g. a pending
    redownload) - the manifest's path is the one the real file is expected to land at, so once
    it does, the metadata is already correct with no regeneration needed.
    """
    manifest: dict[str, tuple[str, str, str]] = {}
    if not path.exists():
        return manifest
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            species_id = (row.get("species_id") or "").strip()
            if not species_id:
                continue
            latin_name = (row.get("latin_name") or "").strip()
            cds_filepath = (row.get("cds_filepath") or "").strip()
            gff_filepath = (row.get("gff_filepath") or "").strip()
            manifest[species_id] = (latin_name, cds_filepath, gff_filepath)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species-table", type=Path, default=SPECIES_TABLE)
    parser.add_argument("--transcriptomes-dir", type=Path, default=TRANSCRIPTOMES_DIR)
    parser.add_argument("--genome-manifest", type=Path, default=GENOME_MANIFEST)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    species_by_id = read_species_table(args.species_table)
    fastas_by_id = find_transcriptome_fastas(args.transcriptomes_dir)
    genome_manifest = read_genome_manifest(args.genome_manifest)

    known_ids = set(species_by_id)
    transcriptome_ids = set(fastas_by_id)
    transcriptome_usable_ids = sorted(known_ids & transcriptome_ids)
    missing_fasta = known_ids - transcriptome_ids
    missing_table_entry = sorted(transcriptome_ids - known_ids)

    # Genome-sourced species fill in exactly the gap left by missing transcriptomes - a species
    # with both a transcriptome and a genome manifest entry always prefers the transcriptome.
    genome_usable_ids = sorted(missing_fasta & set(genome_manifest))
    still_missing = sorted(missing_fasta - set(genome_manifest))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["species", "fasta_filename", "latin_name", "gff_filename"])
        for sample_id in transcriptome_usable_ids:
            # Append the 4-letter sample ID to the latin name (e.g. "Genus species AAXJ") so it
            # stays unique and traceable back to its source sample even when two different
            # samples share the same latin species name.
            latin_name = f"{species_by_id[sample_id]} {sample_id}"
            writer.writerow([sample_id, str(fastas_by_id[sample_id]), latin_name, ""])
        for sample_id in genome_usable_ids:
            latin_name, cds_filepath, gff_filepath = genome_manifest[sample_id]
            writer.writerow([sample_id, cds_filepath, latin_name, gff_filepath])

    total = len(transcriptome_usable_ids) + len(genome_usable_ids)
    print(
        f"Wrote {total} species to {args.output} "
        f"({len(transcriptome_usable_ids)} from transcriptomes, {len(genome_usable_ids)} from genomes)"
    )
    if still_missing:
        print(
            f"WARNING: {len(still_missing)} species table entries have no CNGB FASTA file "
            "or genome manifest entry:"
        )
        print(f"  {', '.join(still_missing)}")
    if missing_table_entry:
        print(f"WARNING: {len(missing_table_entry)} CNGB FASTA files have no species table entry:")
        print(f"  {', '.join(missing_table_entry)}")


if __name__ == "__main__":
    main()
