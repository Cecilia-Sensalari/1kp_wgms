#!/usr/bin/env python3
"""Colours of the divergence lines of the Ks plots, by the role of each divergence node in the tested branch.

ksrates draws, in the plot of a focal species, one vertical line per divergence node on the path from the focal
species up to the root (the root itself is the split of the outgroup and is not drawn), and gives all the lines of a
node one colour of the config field divergence_colors: the first colour to the node closest to the focal species, the
second to the next one, and so on. This script writes that list for every focal species of every dataset, so that
the same divergence has the same colour in every plot of a branch. The colour depends on the role of the node in the
local tree of the branch (clades A, B, C, D of focal_branch_to_dataset.tsv; D is the outgroup, i.e. the root):

    node                  colour set in the config (divergence_palette)
    ------------------------------------------------------------
    inside A              Red
    inside B              RoyalBlue
    A|B                   MediumBlue
    inside C              Goldenrod
    (A+B)|C               #A0720A (dark gold)

A dataset can serve several branches (branches that select the same species), and the same node can have another role
in each of them, e.g. the node joining two species is inside A for one branch and inside B for the other. The
colour list of a focal species is then made per branch, and a (dataset, focal species) that gets different lists
is split in one run per list. The output has one row per distinct (dataset, focal species, colour list):

    dataset_id, focal_species, divergence_colors (comma separated), branches (traversal_order of the branches that use
    this list, comma separated), run_tag (empty if the focal species has a single list in the dataset, else b<first
    branch>; it is added to the name of the run folder)

Usage:   make_focal_colorings.py --branches focal_branch_to_dataset.tsv --out focal_colorings.tsv --palette COLORS
Needs ete3 (the same Python environment as make_branch_subtrees.py); it reads only the small per-branch trees of the table,
so it takes seconds.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

ROLES = ("inside A", "inside B", "A|B", "inside C", "(A+B)|C")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--branches", required=True, type=Path, help="focal_branch_to_dataset.tsv")
    parser.add_argument("--out", required=True, type=Path, help="output table")
    parser.add_argument(
        "--palette",
        required=True,
        help=f"{len(ROLES)} comma-separated colours, for: {', '.join(ROLES)}",
    )
    return parser.parse_args()


def split_species(value: str) -> Set[str]:
    return {item for item in value.split(",") if item}


def node_role(leaves: Set[str], group_a: Set[str], group_b: Set[str], group_c: Set[str]) -> int:
    """Index in ROLES of an internal node of the local tree, from the species below it."""
    if leaves <= group_a:
        return 0
    if leaves <= group_b:
        return 1
    if leaves == group_a | group_b:
        return 2
    if leaves <= group_c:
        return 3
    if leaves == group_a | group_b | group_c:
        return 4
    raise ValueError(f"node with species {sorted(leaves)} does not fit the groups A, B, C")


def focal_colors(tree, focal: str, groups: Tuple[Set[str], Set[str], Set[str]], palette: List[str]) -> List[str]:
    """Colour of each divergence node on the path from the focal species up to the root (root excluded)."""
    matches = tree.search_nodes(name=focal)
    if len(matches) != 1:
        raise ValueError(f"focal species {focal} is not exactly once in the tree")
    nodes = matches[0].get_ancestors()[:-1]
    return [palette[node_role(set(node.get_leaf_names()), *groups)] for node in nodes]


def main() -> None:
    args = parse_args()
    palette = [color.strip() for color in args.palette.split(",")]
    if len(palette) != len(ROLES) or not all(palette):
        sys.exit(f"ERROR: --palette needs {len(ROLES)} colours, for: {', '.join(ROLES)}")
    try:
        from ete3 import Tree
    except ImportError:
        sys.exit("ERROR: this script needs ete3 in the active Python environment")

    # (dataset_id, focal species) -> {colour list: [traversal_order of the branches]}, in order of appearance
    lists: Dict[Tuple[str, str], Dict[Tuple[str, ...], List[str]]] = defaultdict(dict)
    with args.branches.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            branch, dataset = row["traversal_order"], row["dataset_id"]
            groups = (split_species(row["A_species"]), split_species(row["B_species"]), split_species(row["C_species"]))
            tree = Tree(row["newick_tree"], format=9)
            for focal in row["focal_species"].split(","):
                try:
                    colors = tuple(focal_colors(tree, focal, groups, palette))
                except ValueError as error:
                    sys.exit(f"ERROR: branch {branch} (dataset {dataset}): {error}")
                if not colors:
                    sys.exit(
                        f"ERROR: branch {branch} (dataset {dataset}): focal species {focal} has no divergence node "
                        "(no outgroup in the tree), ksrates could not analyse it"
                    )
                lists[(dataset, focal)].setdefault(colors, []).append(branch)

    rows_out = []
    n_split = 0
    for (dataset, focal), by_colors in lists.items():
        if len(by_colors) > 1:
            n_split += 1
        for colors, branches in by_colors.items():
            tag = f"b{branches[0]}" if len(by_colors) > 1 else ""
            rows_out.append([dataset, focal, ",".join(colors), ",".join(branches), tag])

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["dataset_id", "focal_species", "divergence_colors", "branches", "run_tag"])
        writer.writerows(rows_out)
    print(f"{len(lists)} (dataset, focal species) pairs, {n_split} of them with more than one colour list: "
          f"{len(rows_out)} runs written to {args.out}")


if __name__ == "__main__":
    main()
