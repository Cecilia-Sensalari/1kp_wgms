#!/usr/bin/env python3
"""Build datasets that deliberately contain both measurable and saturated subtracted terms.

`outgroup_saturation_effect.py` asks whether a saturated Ks(sister, outgroup) damages the rate
adjustment, by comparing outgroups WITHIN one run: focal species, sister and all other input
fixed, only the subtracted outgroup varying. That comparison needs a group holding both kinds of
outgroup at once, and the real datasets almost never do - the selection rule exists to prevent
exactly that mixture. Over the 2354 branches of the 1KP tree only 22 hold both a comfortably
measurable term (< Ks 2.0) and a saturated one (> Ks 2.72), and none of them is in the
flowering plants. This script finds the tree positions where such a dataset CAN be built, and
writes the wrapper tables for them.

TWO CONSTRAINTS SHAPE THE SEARCH, and both were learned the hard way.

1. THE DATASET SHAPE IS FIXED. make_focal_colorings.py walks the ancestors of the focal species
   (`get_ancestors()[:-1]`) and requires every one to be exactly A, inside A, inside B, A|B,
   inside C, or A|B|C. So the tree must be (((A,B),C),D) and nothing else: four nested clades,
   with D the sibling of A|B|C. A ladder of five or six successive outgroups - the obvious
   design, and the first one tried - is rejected before any Ks is computed.

   The consequence is that only three outgroup levels are available (B, C, D), and the depth
   D can reach is bounded by needing three nested clades below it. Inside the angiosperms that
   ceiling is about Ks 1.5, which is why the test cannot be run there: 363 angiosperm positions
   have a valid nesting and not one offers a saturated D. The gymnosperms do, because their
   sibling clade is the flowering plants.

2. TOPOLOGY AND DISTANCE COME FROM DIFFERENT TREES. calib_common.load_analysis_tree() is the
   rooted ASTRAL species tree the datasets are built from; load_concat_tree() carries the branch
   lengths the path distances are read off. Deriving the induced subtree from the concatenation
   tree puts the gymnosperm inside the angiosperms - silently, because the newick is still
   well formed.

    python3 design_saturation_test.py --out-dir <dir> [--lineage gymnosperm] [--n 2]

Writes <out-dir>/subtree_datasets.tsv and <out-dir>/focal_branch_to_dataset.tsv, the pair the
wrapper needs beside each other. Point a run's `dataset_list` at the first, in its own base_dir
so the datasets directory of the real analysis stays clean, and leave draw_branch_trees off: the
branch table here is a stub, carrying only what the colourings step reads.

Runs nothing and needs no Ks data; it only walks the tree.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import calib_common as C

# A subtracted term this far inside the limit is a trustworthy anchor; one beyond it is the case
# under test. The gap between them is deliberate - the fern attempt failed because its "healthy"
# outgroups sat at Ks 2.1-2.6, barely inside, and agreed with each other only to ~3%.
HEALTHY = 2.0
SATURATED = 2.72


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", required=True, help="directory to write the two tables into")
    ap.add_argument("--lineage", default="gymnosperm",
                    help="lineage the focal clade must belong to (default: gymnosperm, the only "
                         "one that offers the contrast; 'any' to search the whole tree)")
    ap.add_argument("--n", type=int, default=2, help="how many datasets to write")
    ap.add_argument("--min-nodes", type=int, default=3000,
                    help="minimum paranome nodes for the two A representatives")
    ap.add_argument("--healthy", type=float, default=HEALTHY)
    ap.add_argument("--saturated", type=float, default=SATURATED)
    args = ap.parse_args()

    analysis = C.load_analysis_tree()          # topology: rooted species tree
    _, dist = C.load_concat_tree()             # distances: concatenation tree
    meta = C.species_metadata()
    groups = C.groups(meta)

    quality = {}
    quality_path = os.path.join(C.ROOT, "ks_analysis/1kp/paralog_ks_quality_nodes.tsv")
    with open(quality_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            quality[row["species"].split()[-1]] = int(row["nodes"])

    names = set(analysis.leaves())
    parent, nodes = {}, []

    def walk(node):
        nodes.append(node)
        for child in node.children:
            parent[id(child)] = node
            walk(child)

    walk(analysis)

    def code(leaf):
        return leaf.split("_")[-1] if "_" in leaf else leaf

    def ks(one, two):
        tip_a, tip_b = C.resolve_tip(code(one), meta), C.resolve_tip(code(two), meta)
        if not tip_a or not tip_b:
            return None
        path = dist(tip_a, tip_b)
        return None if path is None else C.path_to_ks(path)

    def sibling(node):
        up = parent.get(id(node))
        if not up:
            return None
        others = [c for c in up.children if c is not node]
        return others[0] if len(others) == 1 else None

    def best_two(node):
        """The two richest paranomes in a clade - they only have to represent it."""
        return sorted(node.leaves(), key=lambda leaf: -quality.get(code(leaf), 0))[:2]

    found = []
    for node in nodes:
        if len(node.leaves()) < 2:
            continue
        b_clade = sibling(node)
        if b_clade is None:
            continue
        up = parent[id(node)]
        c_clade = sibling(up)
        if c_clade is None:
            continue
        grand = parent.get(id(up))
        if grand is None:
            continue
        d_clade = sibling(grand)
        if d_clade is None:
            continue

        a_reps = best_two(node)
        if len(a_reps) < 2 or min(quality.get(code(x), 0) for x in a_reps) < args.min_nodes:
            continue
        if args.lineage != "any":
            if groups.get(C.resolve_tip(code(a_reps[0]), meta)) != args.lineage:
                continue

        b_reps, c_reps, d_reps = best_two(b_clade), best_two(c_clade), best_two(d_clade)
        if not (b_reps and c_reps and d_reps):
            continue
        # every term is read from the SISTER, because Ks(sister, outgroup) is what is subtracted
        sister = a_reps[1]
        terms = {role: [ks(sister, x) for x in reps]
                 for role, reps in (("B", b_reps), ("C", c_reps), ("D", d_reps))}
        if any(v is None for values in terms.values() for v in values):
            continue
        if not (max(terms["B"]) < args.healthy and max(terms["C"]) < args.healthy
                and min(terms["D"]) > args.saturated):
            continue
        found.append((min(terms["D"]) - max(terms["C"]), a_reps, b_reps, c_reps, d_reps, terms))

    if not found:
        sys.exit(f"no tree position in '{args.lineage}' offers both a measurable "
                 f"(< {args.healthy}) and a saturated (> {args.saturated}) subtracted term")

    found.sort(key=lambda item: -item[0])
    print(f"{len(found)} usable tree positions; writing the {min(args.n, len(found))} widest\n")

    os.makedirs(args.out_dir, exist_ok=True)
    dataset_rows, branch_rows = [], []
    # Neighbouring positions often differ only by swapping which clade is A and which is B, which
    # gives the same species set and so recomputes every ortholog pair in a second directory. One
    # dataset per species set is enough: its focal_species lists both A representatives, so the
    # _all table already reports a group for each of them.
    seen_sets, index = set(), 0
    for _, a_reps, b_reps, c_reps, d_reps, terms in found:
        key = frozenset(code(x) for x in a_reps + b_reps + c_reps + d_reps)
        if key in seen_sets:
            continue
        seen_sets.add(key)
        index += 1
        if index > args.n:
            break
        label = f"sat{index}"
        picks = {}
        for role, reps in (("A", a_reps), ("B", b_reps), ("C", c_reps), ("D", d_reps)):
            for leaf in reps:
                picks[leaf] = code(leaf)
        newick = induced(analysis, picks) + ";"
        groups_abc = (set(code(x) for x in a_reps), set(code(x) for x in b_reps),
                      set(code(x) for x in c_reps))
        check_shape(newick, groups_abc, [code(x) for x in a_reps])

        codes = [code(x) for x in a_reps + b_reps + c_reps + d_reps]
        focal = ",".join(code(x) for x in a_reps)
        print(f"{label}: {newick}")
        print(f"   focal {focal}   sister {code(a_reps[1])}")
        for role in ("B", "C", "D"):
            reps = {"B": b_reps, "C": c_reps, "D": d_reps}[role]
            for leaf, value in zip(reps, terms[role]):
                tag = "SATURATED" if value > args.saturated else "measurable"
                print(f"     {role} {code(leaf):10} Ks(sister,outgroup) = {value:5.2f}  {tag}")
        print()

        dataset_rows.append(dict(dataset_id=label, newick_tree=newick, focal_species=focal,
                                 target_species=",".join(codes), n_selected_species=len(codes),
                                 branches=1))
        branch_rows.append(dict(
            dataset_id=label, traversal_order=f"900{index}", newick_tree=newick,
            focal_species=focal, target_species=",".join(codes),
            A_species=",".join(code(x) for x in a_reps),
            B_species=",".join(code(x) for x in b_reps),
            C_species=",".join(code(x) for x in c_reps),
            D_species=",".join(code(x) for x in d_reps),
            n_selected_species=len(codes),
            B_path=f"{dist(C.resolve_tip(code(a_reps[0]), meta), C.resolve_tip(code(a_reps[1]), meta)):.4f}",
            B_reliability="reliable", B_outgroup_term="test-dataset",
            C_path="", C_reliability="reliable", C_outgroup_term="test-dataset",
            D_path="", D_reliability="reliable", D_outgroup_term="test-dataset",
            rate_adjustment_reliability="reliable"))

    write(os.path.join(args.out_dir, "subtree_datasets.tsv"), dataset_rows)
    write(os.path.join(args.out_dir, "focal_branch_to_dataset.tsv"), branch_rows)
    print(f"wrote both tables into {args.out_dir}")


def induced(node, picks):
    """Newick of the subtree induced on `picks`, collapsing nodes left with one child."""
    if not node.children:
        return picks.get(node.name)
    kept = [induced(child, picks) for child in node.children]
    kept = [item for item in kept if item]
    if not kept:
        return None
    return kept[0] if len(kept) == 1 else "(" + ",".join(kept) + ")"


def check_shape(newick, groups_abc, focal_codes):
    """Reject a tree make_focal_colorings.py would reject, before the wrapper ever sees it."""
    import re

    tokens = re.findall(r"[(),]|[^(),;]+", newick.rstrip(";"))
    position = [0]

    def parse():
        if tokens[position[0]] == "(":
            position[0] += 1
            children = []
            while True:
                children.append(parse())
                if tokens[position[0]] == ",":
                    position[0] += 1
                else:
                    break
            position[0] += 1
            return children
        name = tokens[position[0]]
        position[0] += 1
        return name

    tree = parse()
    group_a, group_b, group_c = groups_abc

    def leaves(item):
        return {item} if isinstance(item, str) else set().union(*[leaves(c) for c in item])

    def ancestors(item, target, path=()):
        if isinstance(item, str):
            return list(path) if item == target else None
        for child in item:
            found = ancestors(child, target, path + (item,))
            if found is not None:
                return found
        return None

    for focal in focal_codes:
        for node in ancestors(tree, focal)[1:]:     # the root is not checked, as in the wrapper
            below = leaves(node)
            if not (below == group_a or below < group_a or below <= group_b
                    or below == group_a | group_b or below <= group_c
                    or below == group_a | group_b | group_c):
                raise SystemExit(f"shape check failed for focal {focal}: node {sorted(below)} "
                                 f"fits none of the roles make_focal_colorings.py accepts")


def write(path, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
