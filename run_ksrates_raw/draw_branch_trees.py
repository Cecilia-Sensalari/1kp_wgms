#!/usr/bin/env python3
"""Draw, for the branches served by the ksrates runs of one dataset, the tree with the tested branch highlighted.

A ksrates run (a focal species of a dataset, with one list of divergence colours) can serve one or two tested branches,
and the focal branch of the tree is a property of the branch, so ksrates' own tree PDFs of the run cannot show it. For
every branch listed in the column "branches" of config_manifest.tsv this script draws one of two trees, chosen with --tree:
  equal      the dataset tree with all branches of equal length, in the style of ksrates' tree_<focal>.pdf; drawn here.
             Written to <init_run_dir>/rate_adjustment/<focal>/tree_<focal>.b<branch>.pdf. Needs only the folder of the run
             (created by "ksrates init"); needs --species-metadata.
  distances  the tree with branch lengths in Ks units, as ksrates' tree_<focal>_distances.pdf. ksrates' own plot_tree_rates()
             draws it, with a hook of this script that highlights the branch. Written to
             <init_run_dir>/rate_adjustment/<focal>/tree_<focal>_distances.b<branch>.pdf. Needs the rate-adjustment table and the
             ortholog peaks, i.e. the finished ksrates pipeline of the run, and a ksrates version whose plot_tree_rates() has
             the arguments tree_hook and output_suffix; needs --expert-config.
In both trees
  - the tested branch (the branch above clade A: the branch of the species itself if A has one species) is thick and in the
    colour of the divergences inside A;
  - the species names are in the colour of the node that separates them from clade A: A in the colour of "inside A", B in
    the colour of "A|B", C in the colour of "(A+B)|C", D (outgroup) in black.
It runs with the python3 of the ksrates container (ete3 with its Qt rendering). A run whose requirements are missing is
reported and skipped. Existing files are overwritten.

Usage:  draw_branch_trees.py --tree equal|distances --dataset ID --branches focal_branch_to_dataset.tsv
                             --manifest config_manifest.tsv --palette COLORS
                             [--species-metadata species_metadata.tsv]   (equal)
                             [--expert-config config_expert.txt]         (distances)
(--palette is the divergence_palette of the wrapper: the colours of inside A, inside B, A|B, inside C, (A+B)|C.)
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path
from typing import Callable, Dict, List

from ksrates_batch_common import ksrates_species_name

# No display on the cluster nodes: Qt (used by ete3 to draw) has to render offscreen, as ksrates itself does
# (fc_manipulate_trees.py). Set before ete3 is imported.
os.environ["QT_QPA_PLATFORM"] = "offscreen"

PALETTE_SIZE = 5
OUTGROUP_COLOR = "black"
BRANCH_WIDTH = 2
FOCAL_BRANCH_WIDTH = 8
ADJUSTMENT_TABLE = "adjustment_table_{}.tsv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tree", required=True, choices=["equal", "distances"], help="which tree to draw")
    parser.add_argument("--dataset", required=True, help="dataset_id to draw")
    parser.add_argument("--branches", required=True, type=Path, help="focal_branch_to_dataset.tsv")
    parser.add_argument("--manifest", required=True, type=Path, help="config_manifest.tsv")
    parser.add_argument("--species-metadata", type=Path, help="species_metadata.tsv (latin names), for --tree equal")
    parser.add_argument("--palette", required=True, help=f"{PALETTE_SIZE} comma-separated colours")
    parser.add_argument("--expert-config", type=Path, help="ksrates expert configuration file, for --tree distances")
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


def branch_colors(row: dict, palette: List[str], species_name: Callable[[str], str]) -> Dict[str, str]:
    """Colour of the name of each species of the branch (A, B, C), by species_name(species)."""
    colors: Dict[str, str] = {}
    for key, color in (("A_species", palette[0]), ("B_species", palette[2]), ("C_species", palette[4])):
        colors.update({species_name(species): color for species in split_list(row[key])})
    return colors


def draw_branch_tree(row: dict, latin_names: Dict[str, str], palette: List[str], out_path: Path) -> None:
    from ete3 import NodeStyle, TextFace, Tree, TreeStyle

    group_a = split_list(row["A_species"])
    name_colors = branch_colors(row, palette, str)

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


def make_tree_hook(tested_species: List[str], name_colors: Dict[str, str], branch_color: str) -> Callable:
    """Hook for ksrates' plotting_tree: called with the tree (informal leaf names) after styles and lengths are set."""
    from ete3 import NodeStyle

    def hook(tree) -> None:
        node = (tree & tested_species[0]) if len(tested_species) == 1 else tree.get_common_ancestor(tested_species)
        if set(node.get_leaf_names()) != set(tested_species):
            raise ValueError(f"no branch of the tree leads exactly to the A species {', '.join(tested_species)}")
        style = NodeStyle()  # a copy, ete3 can share one style object between nodes
        style.update(node.img_style)
        style["hz_line_color"] = branch_color
        style["hz_line_width"] = FOCAL_BRANCH_WIDTH
        node.set_style(style)
        for leaf in tree.get_leaves():
            if leaf.name in name_colors:
                leaf.add_feature("name_color", name_colors[leaf.name])  # read by ksrates when it labels the leaves

    return hook


def main() -> None:
    args = parse_args()
    palette = [color.strip() for color in args.palette.split(",")]
    if len(palette) != PALETTE_SIZE or not all(palette):
        sys.exit(f"ERROR: --palette needs {PALETTE_SIZE} colours")
    if args.tree == "equal" and not args.species_metadata:
        sys.exit("ERROR: --tree equal needs --species-metadata")
    if args.tree == "distances" and not args.expert_config:
        sys.exit("ERROR: --tree distances needs --expert-config")
    latin_names = read_latin_names(args.species_metadata) if args.tree == "equal" else {}
    branches = read_branches(args.branches, args.dataset)
    start_dir = Path.cwd()
    if args.tree == "distances":
        from ksrates.plot_tree import plot_tree_rates
        expert_config = str(args.expert_config.resolve())

    drawn = 0
    problems: List[str] = []
    with args.manifest.open(newline="") as handle:
        for run in csv.DictReader(handle, delimiter="\t"):
            if run["dataset_id"] != args.dataset:
                continue
            focal = run["focal_species"]
            run_dir = Path(run["init_run_dir"])
            run_folder = run_dir / "rate_adjustment" / focal
            run_branches = split_list(run.get("branches", ""))
            if args.tree == "equal" and run_branches and not run_folder.is_dir():
                problems.append(f"{run_folder} does not exist (ksrates init not done?), no tree for branches {run_branches}")
                continue
            if args.tree == "distances" and run_branches and not (run_folder / ADJUSTMENT_TABLE.format(focal)).is_file():
                problems.append(f"{focal}: no rate-adjustment table in {run_folder}, ksrates pipeline not finished?")
                continue
            for branch in run_branches:
                row = branches.get(branch)
                if row is None:
                    problems.append(f"branch {branch} of {focal} is not in {args.branches}")
                    continue
                if focal not in branch_colors(row, palette, ksrates_species_name):
                    problems.append(f"{focal} is not an A, B or C species of branch {branch}")
                    continue
                try:
                    if args.tree == "equal":
                        draw_branch_tree(row, latin_names, palette, run_folder / f"tree_{focal}.b{branch}.pdf")
                    else:
                        tested = [ksrates_species_name(species) for species in split_list(row["A_species"])]
                        hook = make_tree_hook(tested, branch_colors(row, palette, ksrates_species_name), palette[0])
                        os.chdir(run_dir)  # ksrates reads and writes relative to the run folder
                        try:
                            plot_tree_rates(run["config_path"], expert_config, None, True, hook, f".b{branch}")
                        finally:
                            os.chdir(start_dir)
                    drawn += 1
                except SystemExit as error:  # ksrates exits (code 1) when the adjustment data are incomplete
                    problems.append(f"branch {branch}, run {focal}: ksrates stopped (exit {error.code}), see its log lines")
                except Exception as error:  # one failing tree must not stop the others
                    problems.append(f"branch {branch}, run {focal}: {error}")

    print(f"dataset {args.dataset}: {drawn} {args.tree} branch trees drawn")
    for problem in problems:
        print(f"PROBLEM: {problem}", file=sys.stderr)
    if problems:
        sys.exit(1)


if __name__ == "__main__":
    main()
