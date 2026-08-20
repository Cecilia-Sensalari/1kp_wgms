#!/usr/bin/env python3
"""Prepare and optionally run 1KP paralog Ks batch inputs for ksrates.

The script splits the 1KP species tree into adjacent leaf batches, writes one
CyVerse species-name list per batch, optionally calls the existing transcriptome
download/processing shell script, and writes one ksrates config per batch.

The ksrates configs are intended for later use with:

    ksrates paralogs-ks-multi configs/paralog_batches/config_paralogs_batch_001.txt

By default, transcriptomes are downloaded into this raw ksrates analysis area:

    /group/esb/cesen/1kp/source_data/2.transcriptomes/1kp/4.processed_transcriptomes/unfiltered


Regenerate configs/lists only
    python3 /group/esb/cesen/1kp/code/1kp_wgms/run_ksrates_raw/run_paralog_batches.py --skip-download

Download/process only batch 1
    python3 /group/esb/cesen/1kp/code/1kp_wgms/run_ksrates_raw/run_paralog_batches.py --start-batch 1 --end-batch 1

Later, from ksrates_raw, run one batch config
    cd /group/esb/cesen/1kp/ks_analysis/1kp/ksrates_raw
    ksrates paralogs-ks-multi configs/paralog_batches/config_paralogs_batch_001.txt
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path


PROJECT_ROOT = Path("/group/esb/cesen/1kp")

# Location 1KP data files
SOURCE_DATA_DIR = PROJECT_ROOT / "source_data"
# - Input tree, species table
TREE_FILE = SOURCE_DATA_DIR / (
    "3.phylogenetic_tree/1kp_trees/"
    "astral_trees_33_percent-FAA_estimated_species_tree.rooted.wgm_suptab3_mrca.nhx.manual_fixes.tree"
)
SPECIES_TABLE = SOURCE_DATA_DIR / "1.species_dataset/1kp_paper_2019_suptab1_species.tsv"
# - Output directories for transcriptomes
TRANSCRIPTOME_OUTPUT_DIR = SOURCE_DATA_DIR / "2.transcriptomes"

# Location paralog Ks analysis
KSRATES_RAW_DIR = PROJECT_ROOT / "ks_analysis/1kp/ksrates_raw"
CONFIG_DIR = KSRATES_RAW_DIR / "configs/paralog_batches"

# Location script to download transcriptomes from CyVerse
DOWNLOAD_SCRIPT = PROJECT_ROOT / "code/1kp_wgms/run_ksrates_raw/download_cyverse_transcriptomes.sh"


@dataclass
class Node:
    name: str = ""
    dist: str = ""
    children: list["Node"] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        return not self.children


def strip_nhx_comments(newick: str) -> str:
    """Remove NHX/Newick comments before parsing."""
    return re.sub(r"\[.*?\]", "", newick)


def parse_label_and_length(text: str, pos: int) -> tuple[str, str, int]:
    start = pos
    while pos < len(text) and text[pos] not in ":,();":
        pos += 1
    label = text[start:pos].strip()

    length = ""
    if pos < len(text) and text[pos] == ":":
        pos += 1
        start = pos
        while pos < len(text) and text[pos] not in ",();":
            pos += 1
        length = text[start:pos].strip()

    return label, length, pos


def parse_newick(newick: str) -> Node:
    """Parse enough Newick for this tree: labels, branch lengths, and topology."""
    text = strip_nhx_comments(newick).strip()
    if text.endswith(";"):
        text = text[:-1]

    def parse_subtree(pos: int) -> tuple[Node, int]:
        if text[pos] == "(":
            pos += 1
            children: list[Node] = []
            while True:
                child, pos = parse_subtree(pos)
                children.append(child)
                if pos >= len(text):
                    raise ValueError("Unexpected end of Newick while parsing children")
                if text[pos] == ",":
                    pos += 1
                    continue
                if text[pos] == ")":
                    pos += 1
                    break
                raise ValueError(f"Unexpected Newick character at position {pos}: {text[pos:pos + 40]!r}")

            label, length, pos = parse_label_and_length(text, pos)
            return Node(label, length, children), pos

        label, length, pos = parse_label_and_length(text, pos)
        if not label:
            raise ValueError(f"Missing leaf label at position {pos}")
        return Node(label, length), pos

    root, pos = parse_subtree(0)
    if pos != len(text):
        raise ValueError(f"Unexpected trailing Newick text at position {pos}: {text[pos:pos + 80]!r}")
    return root


def leaf_order(node: Node) -> list[str]:
    if node.is_leaf:
        return [node.name]

    leaves: list[str] = []
    for child in node.children:
        leaves.extend(leaf_order(child))
    return leaves


def prune_to_leaf_set(node: Node, keep: set[str]) -> Node | None:
    """Prune tree to selected leaves and suppress single-child internal nodes."""
    if node.is_leaf:
        if node.name in keep:
            return Node(node.name, node.dist)
        return None

    children: list[Node] = []
    for child in node.children:
        pruned = prune_to_leaf_set(child, keep)
        if pruned is not None:
            children.append(pruned)

    if not children:
        return None
    if len(children) == 1:
        return children[0]
    return Node(node.name, node.dist, children)


def node_to_newick(node: Node, keep_branch_lengths: bool) -> str:
    if node.is_leaf:
        text = node.name
    else:
        text = "(" + ",".join(node_to_newick(child, keep_branch_lengths) for child in node.children) + ")"
        if node.name:
            text += node.name

    if keep_branch_lengths and node.dist:
        text += f":{node.dist}"
    return text


def read_species_table(path: Path) -> dict[str, str]:
    """Return mapping from 1KP sample ID to Latin/species name."""
    mapping: dict[str, str] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample_id = row["1KP Index ID"].strip()
            species = row["Species"].strip()
            if sample_id and species:
                mapping[sample_id] = species
    return mapping


def species_name_to_cyverse_suffix(species_name: str) -> str:
    """Convert the species-table name to the expected CyVerse path suffix."""
    suffix = species_name.strip()
    suffix = re.sub(r"\s+", "_", suffix)
    return suffix


def cyverse_name(sample_id: str, species_by_id: dict[str, str]) -> str:
    try:
        species = species_by_id[sample_id]
    except KeyError as exc:
        raise KeyError(f"Species ID {sample_id!r} not found in {SPECIES_TABLE}") from exc
    return f"{sample_id}-{species_name_to_cyverse_suffix(species)}"


def processed_fasta_path(cyverse_dir_name: str, transcriptome_output_dir: Path, selection: str) -> Path:
    if selection == "filtered":
        return transcriptome_output_dir / "4.processed_transcriptomes/filtered" / f"{cyverse_dir_name}_filtered.FNA"
    if selection == "unfiltered":
        return transcriptome_output_dir / "4.processed_transcriptomes/unfiltered" / f"{cyverse_dir_name}_unfiltered.FNA"
    raise ValueError("Ksrates config can point to only one FASTA per species; use --config-fasta-selection filtered or unfiltered")


def format_mapping(entries: list[tuple[str, str]], initial_padding: int = 20) -> str:
    width = max(len(key) for key, _ in entries)
    lines: list[str] = []
    for i, (key, value) in enumerate(entries):
        comma = "," if i < len(entries) - 1 else ""
        prefix = "" if i == 0 else " " * initial_padding
        lines.append(f"{prefix}{key:<{width}} : {value}{comma}")
    return "\n".join(lines)


def write_ksrates_config(
    path: Path,
    batch_id: int,
    batch_species: list[str],
    subtree_newick: str,
    species_by_id: dict[str, str],
    transcriptome_output_dir: Path,
    config_fasta_selection: str,
) -> None:
    latin_entries = [(sample_id, species_by_id[sample_id]) for sample_id in batch_species]
    fasta_entries = []
    for sample_id in batch_species:
        name = cyverse_name(sample_id, species_by_id)
        fasta_entries.append((sample_id, str(processed_fasta_path(name, transcriptome_output_dir, config_fasta_selection))))

    focal_species = batch_species[0]
    text = f"""[SPECIES]
focal_species = {focal_species}
# first species in paralog batch {batch_id:03d}; required by ksrates config syntax

newick_tree = {subtree_newick};
# adjacent leaves from the 1KP species tree

latin_names =       {format_mapping(latin_entries)}
# informal names associated to their scientific names through a colon and separated by comma

fasta_filenames =   {format_mapping(fasta_entries)}
gff_filename =
# no GFF is needed because these batch configs use whole-paranome analysis only

peak_database_path = ortholog_peak_db.tsv
ks_list_database_path = ortholog_ks_list_db.tsv

[ANALYSIS SETTING]
paranome = yes
collinearity = no
reciprocal_retention = no

gff_feature =
gff_attribute =

max_number_outgroups = 4
consensus_mode_for_multiple_outgroups = mean among outgroups

[PARAMETERS]
x_axis_max_limit_paralogs_plot = 5
bin_width_paralogs = 0.1
y_axis_max_limit_paralogs_plot = None
num_bootstrap_iterations = 200
divergence_colors = Red, MediumBlue, DarkGoldenrod, ForestGreen, HotPink, DarkCyan, SaddleBrown, Black
x_axis_max_limit_orthologs_plots = 5
bin_width_orthologs = 0.1
max_ks_paralogs = 5
max_ks_orthologs = 10
"""
    path.write_text(text, encoding="utf-8")


def write_species_list(path: Path, batch_species: list[str], species_by_id: dict[str, str]) -> list[str]:
    names = [cyverse_name(sample_id, species_by_id) for sample_id in batch_species]
    path.write_text("\n".join(names) + "\n", encoding="utf-8")
    return names


def run_downloader(download_script: Path, species_list: Path, selection: str, output_dir: Path) -> None:
    command = [
        "bash",
        str(download_script),
        str(species_list),
        selection,
        str(output_dir),
    ]
    subprocess.run(command, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", type=Path, default=TREE_FILE)
    parser.add_argument("--species-table", type=Path, default=SPECIES_TABLE)
    parser.add_argument("--download-script", type=Path, default=DOWNLOAD_SCRIPT)
    parser.add_argument("--config-dir", type=Path, default=CONFIG_DIR)
    parser.add_argument("--transcriptome-output-dir", type=Path, default=TRANSCRIPTOME_OUTPUT_DIR)
    parser.add_argument("--batch-size", type=int, default=20)
    parser.add_argument("--download-selection", choices=["unfiltered", "filtered", "both"], default="unfiltered")
    parser.add_argument("--config-fasta-selection", choices=["unfiltered", "filtered"], default="unfiltered")
    parser.add_argument("--start-batch", type=int, default=1, help="First 1-based batch to process")
    parser.add_argument("--end-batch", type=int, default=None, help="Last 1-based batch to process")
    parser.add_argument("--skip-download", action="store_true", help="Only write species lists and ksrates configs")
    parser.add_argument("--drop-branch-lengths", action="store_true", help="Write topology-only Newick trees in configs")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.batch_size < 2:
        raise ValueError("--batch-size must be at least 2")

    root = parse_newick(args.tree.read_text(encoding="utf-8").strip())
    leaves = leaf_order(root)
    species_by_id = read_species_table(args.species_table)

    missing_from_table = [sample_id for sample_id in leaves if sample_id not in species_by_id]
    if missing_from_table:
        preview = ", ".join(missing_from_table[:10])
        raise RuntimeError(
            f"{len(missing_from_table)} tree leaves are missing from the species table. "
            f"First missing IDs: {preview}"
        )

    args.config_dir.mkdir(parents=True, exist_ok=True)

    num_batches = math.ceil(len(leaves) / args.batch_size)
    end_batch = args.end_batch or num_batches
    if args.start_batch < 1 or end_batch > num_batches or args.start_batch > end_batch:
        raise ValueError(f"Invalid batch range {args.start_batch}..{end_batch}; valid range is 1..{num_batches}")

    manifest_rows: list[dict[str, str]] = []

    for batch_id in range(1, num_batches + 1):
        start = (batch_id - 1) * args.batch_size
        batch_species = leaves[start : start + args.batch_size]
        selected = set(batch_species)
        subtree = prune_to_leaf_set(root, selected)
        if subtree is None:
            raise RuntimeError(f"Could not prune tree for batch {batch_id}")
        subtree_newick = node_to_newick(subtree, keep_branch_lengths=not args.drop_branch_lengths)

        species_list_path = args.config_dir / f"paralog_batch_{batch_id:03d}_species.txt"
        config_path = args.config_dir / f"config_paralogs_batch_{batch_id:03d}.txt"

        expected_names = write_species_list(species_list_path, batch_species, species_by_id)
        write_ksrates_config(
            config_path,
            batch_id,
            batch_species,
            subtree_newick,
            species_by_id,
            args.transcriptome_output_dir,
            args.config_fasta_selection,
        )

        if args.start_batch <= batch_id <= end_batch and not args.skip_download:
            run_downloader(args.download_script, species_list_path, args.download_selection, args.transcriptome_output_dir)

        fasta_paths = [
            str(processed_fasta_path(name, args.transcriptome_output_dir, args.config_fasta_selection))
            for name in expected_names
        ]
        manifest_rows.append(
            {
                "batch_id": f"{batch_id:03d}",
                "config_file": str(config_path),
                "species_list_file": str(species_list_path),
                "num_species": str(len(batch_species)),
                "first_species": batch_species[0],
                "last_species": batch_species[-1],
                "species_ids": ",".join(batch_species),
                "cyverse_names": ",".join(expected_names),
                "config_fasta_files": ",".join(fasta_paths),
            }
        )

    manifest_path = args.config_dir / "paralog_batch_manifest.tsv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]))
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"Tree leaves: {len(leaves)}")
    print(f"Batch size: {args.batch_size}")
    print(f"Config/species-list batches written: {num_batches}")
    print(f"Config directory: {args.config_dir}")
    print(f"Manifest: {manifest_path}")
    if args.skip_download:
        print("Download step: skipped")
    else:
        print(f"Download step: ran for batches {args.start_batch}..{end_batch}")


if __name__ == "__main__":
    main()
