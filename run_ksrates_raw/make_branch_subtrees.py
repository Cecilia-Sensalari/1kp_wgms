#!/usr/bin/env python3
"""Build compact per-dataset subtrees for ksrates.

For each internal non-root branch, the child node is treated as clade A.
The output subtree keeps:
  - up to 4 representatives from A;
  - up to 3 progressively older sibling clades, named B/C/D;
  - up to 2 representatives from each older clade.

TWO SELECTION STRATEGIES are implemented, chosen with --selection:

  paranome (default)
      The current strategy, documented in docs/ortholog_ks_path_distance_calibration.md.
      Species are ranked by MEDIAN BIN HEIGHT - the median of the 20 weighted-node bar
      heights of their own paranome Ks distribution over Ks 0.5-2.5 - which measures
      whether a species' Ks distribution can actually carry the analysis. For the older
      clades B/C/D a candidate must also fall inside a Ks BUDGET, computed from the
      concatenation-tree path distance to the A species through the calibrated power law.
      No candidate is excluded for being too deep and no branch is skipped: instead each
      selection carries a RELIABILITY level, because the deepest branches are exactly the
      ones the 1KP paper made WGM calls on and a method that declines to analyse them
      cannot check them.

  quality
      The earlier strategy: a weighted percentile of BUSCO completeness and TransRate
      score, which is a property of the ASSEMBLY rather than of the Ks distribution.
      Retained because the decision to drop it has not been made. It needs --score-table,
      --genome-busco-table, --busco-weight and --transrate-weight.

Two TSVs are written:
  - focal_branch_to_dataset.tsv (--out): one row per branch, with its A/B/C/D species. This is the link between
    a branch (traversal_order) and its dataset_id, used to interpret the results of a branch.
  - subtree_datasets.tsv (--datasets-out): one row per distinct dataset_id, input of
    prepare_branch_ksrates_configs.py. Several branches can select the very same species set (for example
    the two sister branches below a node): they have the same tree and the same ksrates runs, so they
    are one dataset. Its focal species are the union of those of its branches and its 'branches' column
    lists their traversal_order values.
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence


BUSCO_COL = "% BUSCOs (complete+fragmented)"
TRANSRATE_COL = "TransRate Score"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True, type=Path, help="rooted species tree (Newick)")
    parser.add_argument("--out", required=True, type=Path, help="per-branch table (focal_branch_to_dataset.tsv)")
    parser.add_argument(
        "--datasets-out",
        type=Path,
        default=None,
        help="per-dataset table, one row per distinct dataset_id; default: subtree_datasets.tsv next to --out",
    )
    parser.add_argument("--max-a", required=True, type=int, help="maximum species taken from clade A")
    parser.add_argument("--max-older", required=True, type=int, help="maximum species taken from each older sibling clade")
    parser.add_argument("--older-clades", required=True, type=int, help="number of older sibling clades (B, C, D...)")
    parser.add_argument("--min-species", required=True, type=int, help="datasets with fewer species are skipped")
    parser.add_argument(
        "--max-datasets",
        type=int,
        default=None,
        help=(
            "Stop the tree traversal after writing this many dataset rows, instead of "
            "walking the full tree - for quick testing without paying for the full "
            "traversal/pruning cost. Unset: full run, all qualifying branches."
        ),
    )
    parser.add_argument(
        "--selection",
        choices=("paranome", "quality"),
        default="paranome",
        help=(
            "How to rank candidate species. 'paranome' (default): median bin height of the "
            "species' own paranome Ks distribution, plus the expected-Ks budget for B/C/D. "
            "'quality': the earlier BUSCO/TransRate assembly score."
        ),
    )
    parser.add_argument(
        "--calibration-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "others" / "calibration",
        help="directory holding calib_common.py (the calibration constants and path distances)",
    )
    parser.add_argument(
        "--bin-heights",
        type=Path,
        default=None,
        help="paralog_ks_bin_heights.tsv; default: the one in the calibration data directory",
    )
    parser.add_argument(
        "--budget-floor",
        type=float,
        default=2.0,
        help="lower bound of the expected-Ks budget for B/C/D (default 2.0)",
    )
    parser.add_argument(
        "--budget-delta",
        type=float,
        default=0.5,
        help=(
            "width of the expected-Ks budget above the closest candidate of the clade "
            "(default 0.5, about two SD of the calibration error)"
        ),
    )
    parser.add_argument("--score-table", type=Path, help="1KP species table with the BUSCO and TransRate columns (--selection quality)")
    parser.add_argument(
        "--genome-busco-table",
        type=Path,
        help=(
            "BUSCO (complete+fragmented, run with collect_busco_scores.py) of the genome-based species, "
            "which are not in the 1KP score table. Columns species_id, busco_complete_fragmented."
        ),
    )
    parser.add_argument("--busco-weight", type=float, help="weight of the BUSCO percentile in the quality score (--selection quality)")
    parser.add_argument("--transrate-weight", type=float, help="weight of the TransRate percentile in the quality score (--selection quality)")

    args = parser.parse_args()
    if args.selection == "quality":
        missing = [
            name
            for name, value in (
                ("--score-table", args.score_table),
                ("--genome-busco-table", args.genome_busco_table),
                ("--busco-weight", args.busco_weight),
                ("--transrate-weight", args.transrate_weight),
            )
            if value is None
        ]
        if missing:
            parser.error("--selection quality needs " + ", ".join(missing))
    return args


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


def dataset_id(selected: Sequence[str]) -> str:
    """
    Content-derived identity for a branch/dataset: a short hash of its own full selected
    species set, not the traversal order it happened to be visited in. This makes the id stable
    across reruns even if the tree or traversal changes elsewhere - a clade whose species set is
    unchanged always gets the same id, while a clade whose set actually changes gets a new one
    (rather than silently reusing/renaming an existing dataset directory with different content).
    """
    return hashlib.sha1("_".join(sorted(selected)).encode("utf-8")).hexdigest()[:8]


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

    # Genome species have no TransRate (it needs reads). They are ranked by BUSCO together with the 1KP
    # transcriptomes and scored by the BUSCO percentile alone (weight rescaled to 1), not penalised for the
    # missing TransRate.
    genome_busco: Dict[str, float] = {}
    with args.genome_busco_table.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            species = (row.get("species_id") or "").strip()
            value = to_float(row.get("busco_complete_fragmented", ""))
            if not species or value is None:
                continue
            if species in busco:
                sys.exit(f"ERROR: genome species {species} is also in the 1KP score table")
            genome_busco[species] = value
    busco.update(genome_busco)

    busco_pct = percentile_ranks(busco)
    transrate_pct = percentile_ranks(transrate)
    scores = {}
    for species in set(busco_pct) | set(transrate_pct):
        scores[species] = (
            args.busco_weight * busco_pct.get(species, 0.0)
            + args.transrate_weight * transrate_pct.get(species, 0.0)
        )
    for species in genome_busco:
        scores[species] = busco_pct[species]
    return scores


# --------------------------------------------------------------------- paranome selection
# The current strategy. calib_common.py holds the calibrated power law, the two reliability
# thresholds and the memoised concatenation-tree path distances; importing it keeps exactly one
# copy of those constants in the repository, instead of a second one drifting out of step here.


class ParanomeSelector:
    """Rank species by median bin height, and B/C/D candidates by an expected-Ks budget."""

    def __init__(self, args: argparse.Namespace):
        sys.path.insert(0, str(args.calibration_dir))
        try:
            import calib_common
        except ImportError as error:
            sys.exit(
                f"ERROR: --selection paranome needs calib_common.py in {args.calibration_dir} "
                f"({error}). Pass --calibration-dir, or use --selection quality."
            )
        self.C = calib_common
        bin_heights = args.bin_heights or Path(calib_common.BIN_HEIGHTS_TSV)
        if not Path(bin_heights).exists():
            sys.exit(
                f"ERROR: {bin_heights} not found - run others/calibration/paralog_ks_bin_heights.py "
                f"first, or pass --bin-heights"
            )
        self.scores = calib_common.median_bin_heights(str(bin_heights))
        _, self.dist = calib_common.load_concat_tree()
        self.floor, self.delta = args.budget_floor, args.budget_delta

    # -- A ---------------------------------------------------------------------------
    def a_species(self, node, max_a: int, order: Dict[str, int]) -> List[str]:
        """The best A_PER_CHILD species of each direct child subclade, capped at max_a."""
        children = [child for child in node.children if child.get_leaf_names()] or [node]
        picked: List[str] = []
        for child in sorted(children, key=lambda c: first_leaf_index(c, order)):
            ranked = sorted(child.get_leaf_names(), key=lambda s: (-self.scores.get(s, 0.0), order[s]))
            picked.extend(ranked[: self.C.A_PER_CHILD])
        return sorted(picked[:max_a], key=order.get)

    # -- B, C, D ---------------------------------------------------------------------
    def older_species(self, clade, a_species: Sequence[str], n: int, order: Dict[str, int]):
        """Up to n representatives of an older clade, plus the reliability of the best of them.

        Every candidate is scored, none is discarded for depth: the budget decides which are
        close enough to be worth ranking on bin height, and if nothing is inside it the budget
        has already opened around the closest candidate, so a representative is always returned.
        """
        paths = {}
        for name in clade.get_leaf_names():
            path = self.C.median_path(name, a_species, self.dist)
            if path is not None:
                paths[name] = path
        if not paths:
            return [], None, None

        expected = {name: self.C.path_to_ks(path) for name, path in paths.items()}
        budget = max(self.floor, min(expected.values()) + self.delta)
        inside = [name for name, ks in expected.items() if ks <= budget] or list(expected)
        chosen = sorted(inside, key=lambda s: (-self.scores.get(s, 0.0), order[s]))[:n]

        # the flag describes the closest chosen representative: that is the pair whose ortholog
        # Ks peak anchors the rate adjustment, and the one most likely to be measurable
        best_path = min(paths[name] for name in chosen)
        return sorted(chosen, key=order.get), best_path, self.C.path_reliability(best_path)


WORST_FIRST = {"unreliable": 0, "borderline": 1, "reliable": 2}


def overall_reliability(levels: Iterable[Optional[str]]) -> str:
    """The best level among B/C/D: one usable outgroup is enough to anchor the adjustment."""
    present = [level for level in levels if level]
    if not present:
        return "none"
    return max(present, key=lambda level: WORST_FIRST[level])


def rows_for_tree(
    tree,
    args: argparse.Namespace,
    scores: Dict[str, float],
    selector: Optional["ParanomeSelector"] = None,
) -> Iterable[List[str]]:
    order = leaf_order(tree)
    older_labels = ["B", "C", "D"]
    branch_index = 0

    for node in tree.traverse("preorder"):
        if node.is_root():
            continue

        branch_index += 1
        paths: Dict[str, Optional[float]] = {}
        levels: Dict[str, Optional[str]] = {}
        if selector is None:
            groups = {"A": a_species(node, args.max_a, order, scores)}
            for label, clade in zip(
                older_labels,
                older_sibling_clades(node, args.older_clades),
            ):
                groups[label] = balanced_species(clade, args.max_older, order, scores)
        else:
            groups = {"A": selector.a_species(node, args.max_a, order)}
            for label, clade in zip(
                older_labels,
                older_sibling_clades(node, args.older_clades),
            ):
                picked, path, level = selector.older_species(
                    clade, groups["A"], args.max_older, order
                )
                groups[label] = picked
                paths[label], levels[label] = path, level

        selected: List[str] = []
        for label in ["A", *older_labels]:
            selected.extend(groups.get(label, []))
        selected = sorted(dict.fromkeys(selected), key=order.get)
        focal_species = []
        for label in ["A", "B", "C"]:
            focal_species.extend(groups.get(label, []))
        focal_species = sorted(dict.fromkeys(focal_species), key=order.get)

        if len(selected) < args.min_species:
            continue

        row = [
            dataset_id(selected),
            branch_index,
            pruned_newick(tree, selected),
            ",".join(focal_species),
            ",".join(selected),
            ",".join(groups["A"]),
            ",".join(groups.get("B", [])),
            ",".join(groups.get("C", [])),
            ",".join(groups.get("D", [])),
            str(len(selected)),
        ]
        # the reliability columns are empty under --selection quality, which knows nothing
        # about path distance
        for label in older_labels:
            path = paths.get(label)
            row.append("" if path is None else f"{path:.4f}")
            row.append(levels.get(label) or "")
        row.append(overall_reliability(levels.values()) if selector is not None else "")
        yield row


BRANCH_COLUMNS = [
    "dataset_id",
    "traversal_order",
    "newick_tree",
    "focal_species",
    "target_species",
    "A_species",
    "B_species",
    "C_species",
    "D_species",
    "n_selected_species",
    "B_path",
    "B_reliability",
    "C_path",
    "C_reliability",
    "D_path",
    "D_reliability",
    "rate_adjustment_reliability",
]
DATASET_COLUMNS = [
    "dataset_id",
    "newick_tree",
    "focal_species",
    "target_species",
    "n_selected_species",
    "branches",
]


def merge_branches_into_datasets(
    branch_rows: Sequence[List[str]], order: Dict[str, int]
) -> List[List[str]]:
    """One row per dataset_id, in order of first appearance.

    Rows with the same dataset_id have the same selected species, so the same pruned tree and target
    list; this is checked. Their focal species can differ (e.g. when only the outgroup D differs), and
    the dataset gets the union of them.
    """
    merged: Dict[str, dict] = {}
    for row in branch_rows:
        dataset, newick, focal, target, n_selected = row[0], row[2], row[3], row[4], row[9]
        # columns past n_selected_species are per-branch flags and do not define a dataset
        traversal_order = str(row[1])
        entry = merged.get(dataset)
        if entry is None:
            merged[dataset] = {
                "newick": newick,
                "target": target,
                "n": n_selected,
                "focal": set(focal.split(",")),
                "branches": [traversal_order],
            }
            continue
        if entry["newick"] != newick or entry["target"] != target:
            sys.exit(
                f"ERROR: branches {entry['branches'][0]} and {traversal_order} have the same "
                f"dataset_id {dataset} but different trees or species"
            )
        entry["focal"].update(focal.split(","))
        entry["branches"].append(traversal_order)

    return [
        [
            dataset,
            entry["newick"],
            ",".join(sorted(entry["focal"], key=order.get)),
            entry["target"],
            entry["n"],
            ",".join(entry["branches"]),
        ]
        for dataset, entry in merged.items()
    ]


def main() -> None:
    args = parse_args()

    print("Loading 1KP tree...", flush=True)
    tree = load_tree(args.tree)
    print(f"Tree loaded: {len(tree.get_leaf_names())} species", flush=True)
    
    selector = None
    scores: Dict[str, float] = {}
    if args.selection == "quality":
        print("Loading quality scores (1KP table and genome BUSCO)...", flush=True)
        scores = load_quality_scores(args)
        print(f"Quality scores loaded for {len(scores)} species", flush=True)
    else:
        print("Loading paranome bin heights and concatenation-tree path distances...", flush=True)
        selector = ParanomeSelector(args)
        print(
            f"Median bin height for {len(selector.scores)} species; "
            f"expected Ks = (path/{selector.C.POWER_A})^{selector.C.POWER_EXPONENT}; "
            f"budget = max({args.budget_floor}, min + {args.budget_delta}); "
            f"reliable to path {selector.C.RELIABLE_PATH}, borderline to {selector.C.BORDERLINE_PATH}",
            flush=True,
        )

    datasets_out = args.datasets_out or args.out.with_name("subtree_datasets.tsv")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    datasets_out.parent.mkdir(parents=True, exist_ok=True)

    print("Selecting species and pruning a subtree for every branch...", flush=True)
    branch_rows: List[List[str]] = []
    for row in rows_for_tree(tree, args, scores, selector):
        branch_rows.append(row)
        if len(branch_rows) % 100 == 0:
            print(f"  {len(branch_rows)} branches done", flush=True)
        if args.max_datasets is not None and len(branch_rows) >= args.max_datasets:
            break

    print("Merging branches with the same species set into datasets...", flush=True)
    dataset_rows = merge_branches_into_datasets(branch_rows, leaf_order(tree))

    print("Writing the tables...", flush=True)
    for path, columns, rows in (
        (args.out, BRANCH_COLUMNS, branch_rows),
        (datasets_out, DATASET_COLUMNS, dataset_rows),
    ):
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(columns)
            writer.writerows(rows)

    print(f"Wrote {len(branch_rows)} branches to {args.out}")
    print(f"Wrote {len(dataset_rows)} distinct datasets to {datasets_out}")

    if selector is not None:
        counts = collections.Counter(row[-1] for row in branch_rows)
        total = sum(counts.values()) or 1
        print("Rate-adjustment reliability of the branches (best of B/C/D):")
        for level in ("reliable", "borderline", "unreliable", "none"):
            if counts[level]:
                print(f"  {level:<11}{counts[level]:>6}  ({100 * counts[level] / total:4.1f}%)")


if __name__ == "__main__":
    main()
