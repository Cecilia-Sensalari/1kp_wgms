#!/usr/bin/env python3
"""Do the outgroup-term flags say something true about the data?

make_branch_subtrees.py decides, from the TREE ALONE, whether the subtracted term Ks(B,O) of the
rate adjustment

    Ks_adjusted = ( Ks(A,B) + Ks(A,O) - Ks(B,O) ) / 2

will be measurable, and records the verdict per branch in B_outgroup_term / C_outgroup_term /
D_outgroup_term as

    ok | saturated:<roles> | saturated:<roles>;partly:<roles>

where a `saturated:` role fails for EVERY representative picked from that clade and a `partly:`
role for only some of them. "Fails" means a concat-tree path distance above --outgroup-term-limit
(1.0, expected Ks 2.72): beyond that an ortholog distribution no longer has a usable peak
(ortholog_ks_path_distance_calibration.md, 3.5).

That is a PREDICTION. This script tests it against the peaks actually stored in the Ks database:
for every (role, outgroup) pair of every branch it compares what the flag claims with the measured
mode, and reports where the two disagree. Everything downstream - comparing the spread of the
adjusted estimates between flagged and unflagged branches - only means something if the flag is
describing the data it claims to describe.

    python3 check_outgroup_term_flags.py --table <focal_branch_to_dataset.tsv> \
        --peaks <...ortholog_peaks_*.tsv> [--branches 664-673] [--limit 2.72]

The peak export is the one `ksrates inspect-ks-db <address file> --type ortholog` writes. Its
species are Latin names ending in the 1KP code ("Danaea nodosa DFHO"), which is how they are
matched back to the codes in the branch table. Reads no Ks data of its own and runs nothing.
"""

from __future__ import annotations

import argparse
import csv
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import calib_common as C


def parse_flag(value):
    """'saturated:C,D;partly:B' -> ({'C','D'}, {'B'}). Empty or 'ok' -> (set(), set())."""
    saturated, partly = set(), set()
    for part in (value or "").split(";"):
        part = part.strip()
        if part.startswith("saturated:"):
            saturated = {r for r in part[len("saturated:"):].split(",") if r}
        elif part.startswith("partly:"):
            partly = {r for r in part[len("partly:"):].split(",") if r}
    return saturated, partly


def load_peaks(path):
    """{(code1, code2) sorted: mode}. The export has an unnamed index column first."""
    peaks = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                mode = float(row["Mode"])
            except (TypeError, ValueError):
                continue
            code1 = row["Species1"].split()[-1]
            code2 = row["Species2"].split()[-1]
            peaks[tuple(sorted((code1, code2)))] = mode
    return peaks


def codes(field):
    return [s.split()[-1] for s in (field or "").split(",") if s]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True, help="focal_branch_to_dataset.tsv with the flags")
    ap.add_argument("--peaks", required=True, help="ortholog peak export of the Ks database")
    ap.add_argument("--branches", help="traversal orders to check, e.g. 664-673 or 664,668 "
                                       "(default: every branch in the table)")
    ap.add_argument("--limit", type=float, default=None,
                    help="measured Ks above which a peak counts as saturated "
                         "(default: the Ks of RELIABLE_PATH)")
    args = ap.parse_args()

    limit = args.limit if args.limit is not None else C.path_to_ks(C.RELIABLE_PATH)
    peaks = load_peaks(args.peaks)

    wanted = None
    if args.branches:
        wanted = set()
        for piece in args.branches.split(","):
            if "-" in piece:
                lo, hi = piece.split("-")
                wanted.update(str(n) for n in range(int(lo), int(hi) + 1))
            else:
                wanted.add(piece.strip())

    print(f"{len(peaks)} ortholog peaks from {os.path.basename(args.peaks)}")
    print(f"a peak counts as saturated at measured Ks > {limit:.2f}\n")

    agree = disagree = missing = 0
    rows = []
    with open(args.table) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if wanted is not None and row["traversal_order"] not in wanted:
                continue
            picks = {r: codes(row[f"{r}_species"]) for r in ("B", "C", "D")}
            for role in ("B", "C"):
                saturated, partly = parse_flag(row[f"{role}_outgroup_term"])
                for outgroup in ("C", "D"):
                    if outgroup <= role or not picks[role] or not picks[outgroup]:
                        continue
                    measured = []
                    for a in picks[role]:
                        for b in picks[outgroup]:
                            mode = peaks.get(tuple(sorted((a, b))))
                            measured.append((a, b, mode))
                    present = [m for _, _, m in measured if m is not None]
                    if not present:
                        missing += 1
                        rows.append((row["traversal_order"], role, outgroup,
                                     row[f"{role}_outgroup_term"] or "ok", "-", "-", "NO DATA"))
                        continue
                    # the flag speaks about every pick of the clade, so compare the extremes
                    n_over = sum(1 for m in present if m > limit)
                    if n_over == len(present):
                        observed = "all over"
                    elif n_over:
                        observed = "some over"
                    else:
                        observed = "none over"
                    claimed = ("all over" if outgroup in saturated else
                               "some over" if outgroup in partly else "none over")
                    verdict = "agree" if claimed == observed else "DISAGREE"
                    if verdict == "agree":
                        agree += 1
                    else:
                        disagree += 1
                    rows.append((row["traversal_order"], role, outgroup,
                                 row[f"{role}_outgroup_term"] or "ok",
                                 f"{min(present):.2f}-{max(present):.2f}",
                                 f"{n_over}/{len(present)}", verdict))

    header = ("branch", "role", "outgrp", "flag", "measured Ks", "over", "verdict")
    print(f"{header[0]:>7} {header[1]:>5} {header[2]:>7}  {header[3]:<22}"
          f"{header[4]:>14}{header[5]:>7}  {header[6]}")
    for r in rows:
        print(f"{r[0]:>7} {r[1]:>5} {r[2]:>7}  {r[3]:<22}{r[4]:>14}{r[5]:>7}  {r[6]}")

    total = agree + disagree
    print(f"\n{total} (branch, role, outgroup) combinations with data: "
          f"{agree} agree, {disagree} disagree"
          + (f"; {missing} had no peak in the database" if missing else ""))

    over = [float(r[4].split("-")[1]) for r in rows if r[4] != "-" and "DISAGREE" not in r[6]]
    if over:
        print(f"median of the largest measured Ks per combination: {statistics.median(over):.2f}")


if __name__ == "__main__":
    main()
