#!/usr/bin/env python3
"""Draw, for the branches served by the ksrates runs of one dataset, the dataset tree with the tested branch highlighted.

A ksrates run (a focal species of a dataset, with one list of divergence colours) can serve one or two tested branches,
and the focal branch of the tree is a property of the branch, so ksrates' own tree PDF of the run cannot show it. For
every branch listed in the column "branches" of config_manifest.tsv this script draws the tree of the branch in the
same style as ksrates' tree_<focal>.pdf, with
  - the tested branch (the branch above clade A: the branch of the species itself if A has one species) in the colour of
    the divergences inside A, thick;
  - the species names in the colour of the node that separates them from clade A: A in the colour of "inside A", B in
    the colour of "A|B", C in the colour of "(A+B)|C", D (outgroup) in black.
The file is written next to the plots of the run:
  <init_run_dir>/rate_adjustment/<focal>/tree_<focal>.b<branch>.pdf
It needs ete3 with its Qt rendering, i.e. the python3 of the ksrates container, and the folder of the run must already
exist (it is created by "ksrates init"); a run without it is reported and skipped. Existing files are overwritten.

Usage:  draw_branch_trees.py --dataset ID --branches focal_branch_to_dataset.tsv --manifest config_manifest.tsv
                             --species-metadata species_metadata.tsv --palette COLORS
(--palette is the divergence_palette of the wrapper: the colours of inside A, inside B, A|B, inside C, (A+B)|C.)
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path
from typing import Dict, List

from ksrates_batch_common import ksrates_species_name

# No display on the cluster nodes: Qt (used by ete3 to draw) has to render offscreen, as ksrates itself does
# (fc_manipulate_trees.py). Set before ete3 is imported.
os.environ["QT_QPA_PLATFORM"] = "offscreen"

PALETTE_SIZE = 5
OUTGROUP_COLOR = "black"
BRANCH_WIDTH = 2
FOCAL_BRANCH_WIDTH = 8


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dataset", required=True, help="dataset_id to draw")
    parser.add_argument("--branches", required=True, type=Path, help="focal_branch_to_dataset.tsv")
    parser.add_argument("--manifest", required=True, type=Path, help="config_manifest.tsv")
    parser.add_argument("--species-metadata", required=True, type=Path, help="species_metadata.tsv (latin names)")
    parser.add_argument("--palette", required=True, help=f"{PALETTE_SIZE} comma-separated colours")
    return parser.parse_args()


def split_list(value: str) -> List[str]:
    return [item for item in value.split(",") if item]


def read_latin_names(path: Path) -> Dict[str, str]:
    with path.open(newline="") as handle:
        return {row["species"]: row["latin_name"] for row in csv.DictReader(handle, delimiter="\t")}


def read_branches(path: Path, dataset: str) -> Dict[str, dict]:
    """The rows of the dataset in the branch table, by traversal_order."""
    with path.open(newline="") as handle:
        return {
            row["traversal_order"]: row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["dataset_id"] == dataset
        }


def draw_branch_tree(row: dict, latin_names: Dict[str, str], palette: List[str], out_path: Path) -> None:
    from ete3 import NodeStyle, TextFace, Tree, TreeStyle

    group_a, group_b, group_c = (split_list(row[key]) for key in ("A_species", "B_species", "C_species"))
    name_colors = {}
    name_colors.update({species: palette[0] for species in group_a})
    name_colors.update({species: palette[2] for species in group_b})
    name_colors.update({species: palette[4] for species in group_c})

    tree = Tree(row["newick_tree"], format=9)
    missing = [leaf.name for leaf in tree.get_leaves() if leaf.name not in latin_names]
    if missing:
        raise ValueError(f"no latin name in the species metadata for: {', '.join(missing)}")

    focal_node = tree & group_a[0] if len(group_a) == 1 else tree.get_common_ancestor(group_a)
    if set(focal_node.get_leaf_names()) != set(group_a):
        raise ValueError(f"no branch of the tree leads exactly to the A species {', '.join(group_a)}")

    # Same look as the tree that ksrates draws (fc_manipulate_trees.plot_uncorrected_phylogeny): equal branch lengths
    for node in tree.traverse():
        style = NodeStyle()
        style["hz_line_width"] = BRANCH_WIDTH
        style["vt_line_width"] = BRANCH_WIDTH
        style["size"] = 0
        if node is focal_node:
            style["hz_line_color"] = palette[0]
            style["hz_line_width"] = FOCAL_BRANCH_WIDTH
        node.set_style(style)
        node.dist = 0.125 if node.is_root() else 0.25
    for leaf in tree.get_leaves():
        label = TextFace(f"{latin_names[leaf.name]} ", ftype="Arial", fsize=18, fstyle="italic",
                         fgcolor=name_colors.get(leaf.name, OUTGROUP_COLOR))
        leaf.add_face(label, column=0, position="branch-right")

    tree_style = TreeStyle()
    tree_style.orientation = 1
    tree_style.branch_vertical_margin = 14
    tree_style.show_leaf_name = False
    tree_style.show_branch_length = False
    tree_style.margin_left = tree_style.margin_right = tree_style.margin_top = tree_style.margin_bottom = 25
    tree_style.scale = 200
    tree_style.show_scale = False
    tree.render(str(out_path), w=4.5, units="in", tree_style=tree_style)


def main() -> None:
    args = parse_args()
    palette = [color.strip() for color in args.palette.split(",")]
    if len(palette) != PALETTE_SIZE or not all(palette):
        sys.exit(f"ERROR: --palette needs {PALETTE_SIZE} colours")
    latin_names = read_latin_names(args.species_metadata)
    branches = read_branches(args.branches, args.dataset)

    drawn = 0
    problems: List[str] = []
    with args.manifest.open(newline="") as handle:
        for run in csv.DictReader(handle, delimiter="\t"):
            if run["dataset_id"] != args.dataset:
                continue
            focal = run["focal_species"]
            run_folder = Path(run["init_run_dir"]) / "rate_adjustment" / focal
            run_branches = split_list(run.get("branches", ""))
            if run_branches and not run_folder.is_dir():
                problems.append(f"{run_folder} does not exist (ksrates init not done?), no tree for branches {run_branches}")
                continue
            for branch in run_branches:
                row = branches.get(branch)
                if row is None:
                    problems.append(f"branch {branch} of {focal} is not in {args.branches}")
                    continue
                tested = split_list(row["A_species"]) + split_list(row["B_species"]) + split_list(row["C_species"])
                if focal not in {ksrates_species_name(species) for species in tested}:
                    problems.append(f"{focal} is not an A, B or C species of branch {branch}")
                    continue
                try:
                    draw_branch_tree(row, latin_names, palette, run_folder / f"tree_{focal}.b{branch}.pdf")
                    drawn += 1
                except Exception as error:  # one failing tree must not stop the others
                    problems.append(f"branch {branch}, run {focal}: {error}")

    print(f"dataset {args.dataset}: {drawn} branch trees drawn")
    for problem in problems:
        print(f"PROBLEM: {problem}", file=sys.stderr)
    if problems:
        sys.exit(1)


if __name__ == "__main__":
    main()
