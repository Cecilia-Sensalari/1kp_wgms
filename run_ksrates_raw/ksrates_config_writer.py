"""Writing per-(dataset, focal species) ksrates config files."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

from ksrates_batch_common import DatasetRecord, SpeciesRecord


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
