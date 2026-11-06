#!/usr/bin/env python3
"""How often does the rate adjustment subtract a term it cannot measure?

The adjustment of a focal species A by a sister B with outgroup O is

    Ks_adjusted = ( Ks(A,B) + Ks(A,O) - Ks(B,O) ) / 2

so three ortholog peaks have to be measurable, not one. The selection rule constrains only
Ks(A,B): a candidate B is kept when its expected Ks to A is inside
`max(BUDGET_FLOOR, min + BUDGET_DELTA)`. Nothing in it looks at Ks(B,O), the term that is
SUBTRACTED, although that term moves when the candidate changes.

The sister-rate test (ks_analysis/1kp/sister_rate_test) showed why this matters: in branch 667
every candidate passed the budget while three of them had an outgroup term beyond 2.7, where an
ortholog distribution no longer has a peak (see ortholog_ks_path_distance_calibration.md, 3.5),
and that branch had by far the largest spread between the estimates its sisters produced.

The offset is predictable, in PATH units. Path distances add along the tree, so with S the node
where A and B split and a, b, c the paths from S to A, B and O,

    path(A,B) = a + b,   path(B,O) = b + c
    delta := path(B,O) - path(A,B) = c - a          <- b cancels EXACTLY

delta is therefore a property of the focal species and the outgroup alone, and a constraint on
the subtracted term is a shifted constraint on the one the rule already applies:

    path(B,O) <= RELIABLE_PATH   <=>   path(A,B) <= RELIABLE_PATH - delta

The same difference taken in Ks is NOT constant, because Ks = (path/A)^c is convex: the Ks gap
grows with b (0.96 at b=0.10, 1.10 at b=0.45, 1.25 at b=1.00 for one real branch). That is the
drift of up to 0.13 Ks seen across the candidates of a sister-rate branch. The arithmetic
belongs in path space; Ks values are reported for reading only.

This script measures, over all branches of the 1KP tree, how often delta is positive (the only
case in which the current rule is blind to an unmeasurable subtraction), and how many branches
actually have one.

    python3 outgroup_term_survey.py [--table ...] [--limit 2.7]

Reads the selection table written by `make_branch_subtrees.py --selection paranome` and the
concatenation tree; runs nothing and needs no Ks data.
"""

from __future__ import annotations

import argparse
import csv
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import calib_common as C

DEFAULT_TABLE = os.path.join(C.ROOT, "ks_analysis/1kp/focal_branch_to_dataset_paranome.tsv")


def split_species(field):
    return [s for s in (field or "").split(",") if s]


def expected_ks(dist, tip_a, tip_b):
    """Expected ortholog Ks between two tips, or None when the path is unavailable."""
    path = dist(tip_a, tip_b)
    return None if path is None else C.path_to_ks(path)


def worst_pair(dist, tips_1, tips_2):
    """Largest expected Ks over the cross pairs of two clades, and the pair that gives it."""
    best = None
    for t1 in tips_1:
        for t2 in tips_2:
            ks = expected_ks(dist, t1, t2)
            if ks is not None and (best is None or ks > best[0]):
                best = (ks, t1, t2)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", default=DEFAULT_TABLE)
    ap.add_argument("--limit", type=float, default=2.7,
                    help="expected Ks beyond which a peak is not measurable (default 2.7, the "
                         "Ks of RELIABLE_PATH)")
    ap.add_argument("--out", help="write the per-branch table here")
    args = ap.parse_args()

    _, dist = C.load_concat_tree()
    meta = C.species_metadata()
    resolved = {}

    def tips(field):
        out = []
        for name in split_species(field):
            if name not in resolved:
                resolved[name] = C.resolve_tip(name.split()[-1], meta)
            if resolved[name]:
                out.append(resolved[name])
        return out

    rows = []
    with open(args.table) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            a, b = tips(row["A_species"]), tips(row["B_species"])
            if not a or not b:
                continue
            ab = worst_pair(dist, a, b)
            if ab is None:
                continue
            for role in ("C", "D"):
                o = tips(row[f"{role}_species"])
                if not o:
                    continue
                bo = worst_pair(dist, b, o)
                ao = worst_pair(dist, a, o)
                if bo is None or ao is None:
                    continue
                # delta on the pair that actually binds: the worst B-O against the worst A-B.
                # Exact in path units; the Ks difference is kept only because it reads easily.
                path_ab = dist(ab[1], ab[2])
                path_bo = dist(bo[1], bo[2])
                rows.append(dict(dataset=row["dataset_id"], order=row["traversal_order"],
                                 role=role, ks_ab=ab[0], ks_ao=ao[0], ks_bo=bo[0],
                                 path_ab=path_ab, path_bo=path_bo,
                                 delta_path=path_bo - path_ab, delta_ks=bo[0] - ab[0],
                                 ceiling_path=C.RELIABLE_PATH - (path_bo - path_ab),
                                 sister=ab[2], outgroup=bo[2]))

    if not rows:
        print("no rows: is the selection table the paranome one?")
        return

    limit = args.limit
    print(f"{len(rows)} (branch, outgroup role) combinations from {args.table}")
    print(f"measurability limit: expected Ks {limit}\n")

    print(f"{'role':>5}{'n':>7}{'delta>0':>10}{'median delta':>14}"
          f"{'Ks(B,O) > limit':>18}{'blind':>8}")  # delta in path units
    for role in ("C", "D", None):
        group = [r for r in rows if role is None or r["role"] == role]
        if not group:
            continue
        positive = sum(1 for r in group if r["delta_path"] > 0)
        over = [r for r in group if r["ks_bo"] > limit]
        # blind: the rule sees nothing wrong (A-B measurable) but the subtracted term is not
        blind = [r for r in group if r["ks_bo"] > limit >= r["ks_ab"]]
        label = role or "both"
        print(f"{label:>5}{len(group):>7}{100 * positive / len(group):>9.1f}%"
              f"{statistics.median(r['delta_path'] for r in group):>+14.2f}"
              f"{100 * len(over) / len(group):>17.1f}%{100 * len(blind) / len(group):>7.1f}%")

    print("\ndelta (path units) distribution over all combinations:")
    deltas = sorted(r["delta_path"] for r in rows)
    for q in (5, 25, 50, 75, 95):
        print(f"  p{q:<3} {deltas[min(len(deltas) - 1, q * len(deltas) // 100)]:+.2f}")

    worst = sorted(rows, key=lambda r: -r["delta_path"])[:10]
    print("\nlargest positive offsets (the cases the current rule cannot see):")
    print(f"  {'dataset':<10}{'role':>5}{'Ks(A,B)':>9}{'Ks(B,O)':>9}{'delta path':>12}"
          f"  sister / outgroup")
    for r in worst:
        print(f"  {r['dataset']:<10}{r['role']:>5}{r['ks_ab']:>9.2f}{r['ks_bo']:>9.2f}"
              f"{r['delta_path']:>+12.3f}  {r['sister']} / {r['outgroup']}")

    if args.out:
        with open(args.out, "w", newline="") as fh:
            writer = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nper-branch table: {args.out}")


if __name__ == "__main__":
    main()
