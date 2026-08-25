#!/usr/bin/env python3
"""Build compact branch-centered subtrees for ksrates.

For each internal non-root branch, the child node is treated as clade A.
The output subtree keeps:
  - up to 4 balanced representatives from A;
  - up to 3 progressively older sibling clades, named B/C/D;
  - up to 2 balanced representatives from each older clade.

The output TSV can be passed to prepare_branch_ksrates_configs.py.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence


TREE_PATH = Path(
    "/group/esb/cesen/1kp/source_data/3.phylogenetic_tree/1kp_trees/"
    "astral_trees_33_percent-FAA_estimated_species_tree.rooted."
    "wgm_suptab3_mrca.nhx.manual_fixes.tree"
)
DEFAULT_OUT = Path(__file__).with_name("branch_subtrees.tsv")
SCORE_TABLE = Path(
    "/group/esb/cesen/1kp/source_data/1.species_dataset/"
    "1kp_paper_2019_suptab1_species.tsv"
)
BUSCO_COL = "% BUSCOs (complete+fragmented)"
TRANSRATE_COL = "TransRate Score"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", type=Path, default=TREE_PATH)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--max-a", type=int, default=4)
    parser.add_argument("--max-older", type=int, default=2)
    parser.add_argument("--older-clades", type=int, default=3)
    parser.add_argument("--min-species", type=int, default=3)
    parser.add_argument("--score-table", type=Path, default=SCORE_TABLE)
    parser.add_argument("--busco-weight", type=float, default=0.8)
    parser.add_argument("--transrate-weight", type=float, default=0.2)
    return parser.parse_args()


def load_tree(path: Path):
    try:
        from ete3 import Tree
    except ImportError:
        sys.exit("ERROR: this script needs ete3 in the active Python environment")

    try:
        return Tree(str(path), format=1)
    except Exception as error:
        sys.exit(f"ERROR: could not read tree {path}: {error}")


def leaf_order(tree) -> Dict[str, int]:
    return {name: index for index, name in enumerate(tree.get_leaf_names())}


def leaves(node, order: Dict[str, int]) -> List[str]:
    return sorted(node.get_leaf_names(), key=order.get)


def first_leaf_index(node, order: Dict[str, int]) -> int:
    return min(order[name] for name in node.get_leaf_names())


def quality_key(name: str, order: Dict[str, int], scores: Dict[str, float]) -> tuple:
    score = scores.get(name)
    if score is None:
        return (1, 0.0, order[name], name)
    return (0, -score, order[name], name)


def ranked_leaves(node, n: int, order: Dict[str, int], scores: Dict[str, float]) -> List[str]:
    chosen = sorted(node.get_leaf_names(), key=lambda name: quality_key(name, order, scores))[:n]
    return sorted(chosen, key=order.get)


def balanced_species(
    node,
    n: int,
    order: Dict[str, int],
    scores: Dict[str, float],
) -> List[str]:
    """Choose up to n leaves while spreading picks across child subclades."""
    if n <= 0:
        return []
    node_leaves = leaves(node, order)
    if len(node_leaves) <= n:
        return node_leaves
    if n <= 1 or node.is_leaf():
        return ranked_leaves(node, 1, order, scores)

    children = [child for child in node.children if child.get_leaf_names()]
    children.sort(key=lambda child: first_leaf_index(child, order))

    slots = [[child, 0] for child in children]
    for slot in slots:
        if sum(count for _, count in slots) < n:
            slot[1] = 1

    while sum(count for _, count in slots) < n:
        candidates = [
            slot for slot in slots if slot[1] < len(slot[0].get_leaf_names())
        ]
        if not candidates:
            break
        slot = min(
            candidates,
            key=lambda slot: (
                slot[1] / len(slot[0].get_leaf_names()),
                first_leaf_index(slot[0], order),
            ),
        )
        slot[1] += 1

    picked: List[str] = []
    for child, count in slots:
        picked.extend(balanced_species(child, count, order, scores))
    return sorted(picked[:n], key=order.get)


def a_species(node, max_a: int, order: Dict[str, int], scores: Dict[str, float]) -> List[str]:
    """Choose A as up to two representatives from each child subclade."""
    if node.is_leaf():
        return leaves(node, order)

    picked: List[str] = []
    children = [child for child in node.children if child.get_leaf_names()]
    children.sort(key=lambda child: first_leaf_index(child, order))
    for child in children:
        picked.extend(balanced_species(child, min(2, len(child.get_leaf_names())), order, scores))
    return sorted(picked[:max_a], key=order.get)


def older_sibling_clades(a_node, max_clades: int) -> List:
    """Return sibling clades encountered while walking from A toward the root."""
    clades = []
    child_on_path = a_node
    parent = a_node.up

    while parent is not None and len(clades) < max_clades:
        siblings = [child for child in parent.children if child is not child_on_path]
        siblings.sort(key=lambda child: len(child.get_leaf_names()), reverse=True)
        for sibling in siblings:
            clades.append(sibling)
            if len(clades) >= max_clades:
                break
        child_on_path = parent
        parent = parent.up

    return clades


def pruned_newick(tree, species: Sequence[str]) -> str:
    subtree = tree.copy(method="deepcopy")
    subtree.prune(list(species), preserve_branch_length=True)
    return subtree.write(format=9).strip()


def branch_id(node, index: int) -> str:
    if node.is_leaf():
        safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", node.name or "").strip("_")
        return f"{index:06d}_{safe_name}"
    return f"{index:06d}_internal"


def to_float(value: str) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def percentile_ranks(values: Dict[str, float]) -> Dict[str, float]:
    ranked = sorted((value, species) for species, value in values.items())
    n = len(ranked)
    percentiles: Dict[str, float] = {}
    i = 0
    while i < n:
        j = i
        while j < n and ranked[j][0] == ranked[i][0]:
            j += 1
        average_rank = (i + 1 + j) / 2
        for _, species in ranked[i:j]:
            percentiles[species] = average_rank / n
        i = j
    return percentiles


def load_quality_scores(args: argparse.Namespace) -> Dict[str, float]:
    busco: Dict[str, float] = {}
    transrate: Dict[str, float] = {}

    with args.score_table.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            species = (row.get("1KP Index ID") or "").strip()
            if not species:
                continue
            busco_value = to_float(row.get(BUSCO_COL, ""))
            transrate_value = to_float(row.get(TRANSRATE_COL, ""))
            if busco_value is not None:
                busco[species] = busco_value
            if transrate_value is not None:
                transrate[species] = transrate_value

    busco_pct = percentile_ranks(busco)
    transrate_pct = percentile_ranks(transrate)
    scores = {}
    for species in set(busco_pct) | set(transrate_pct):
        scores[species] = (
            args.busco_weight * busco_pct.get(species, 0.0)
            + args.transrate_weight * transrate_pct.get(species, 0.0)
        )
    return scores


def rows_for_tree(
    tree,
    args: argparse.Namespace,
    scores: Dict[str, float],
) -> Iterable[List[str]]:
    order = leaf_order(tree)
    older_labels = ["B", "C", "D"]
    branch_index = 0

    for node in tree.traverse("preorder"):
        if node.is_root():
            continue

        branch_index += 1
        groups = {"A": a_species(node, args.max_a, order, scores)}
        for label, clade in zip(
            older_labels,
            older_sibling_clades(node, args.older_clades),
        ):
            groups[label] = balanced_species(clade, args.max_older, order, scores)

        selected: List[str] = []
        for label in ["A", *older_labels]:
            selected.extend(groups.get(label, []))
        selected = sorted(dict.fromkeys(selected), key=order.get)

        if len(selected) < args.min_species:
            continue

        yield [
            branch_id(node, branch_index),
            pruned_newick(tree, selected),
            ",".join(selected),
            ",".join(selected),
            ",".join(groups["A"]),
            ",".join(groups.get("B", [])),
            ",".join(groups.get("C", [])),
            ",".join(groups.get("D", [])),
            str(len(selected)),
        ]


def main() -> None:
    args = parse_args()

    print("Loading 1KP tree...")
    tree = load_tree(args.tree)
    scores = load_quality_scores(args)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "branch_id",
                "newick_tree",
                "focal_species",
                "target_species",
                "A_species",
                "B_species",
                "C_species",
                "D_species",
                "n_selected_species",
            ]
        )
        count = 0
        for row in rows_for_tree(tree, args, scores):
            writer.writerow(row)
            count += 1

    print(f"Wrote {count} branch subtrees to {args.out}")


if __name__ == "__main__":
    main()
