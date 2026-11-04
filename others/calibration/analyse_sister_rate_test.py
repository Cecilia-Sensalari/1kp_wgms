#!/usr/bin/env python3
"""Read the sister-rate test back: does the rate adjustment return the same divergence?

Every sister of a given branch split from the focal species at the same node, so every row of
that branch's `adjustment_table_<focal>.tsv` is an estimate of ONE quantity. Any spread between
them is error, and the question is what the error tracks:

  * the sister's own expected Ks (its rate), which is what the budget rule reads;
  * the subtracted term Ks(sister, outgroup), which is what saturates first;
  * the sister's median bin height, which is what the rule falls back on inside the budget.

If the spread is flat in expected Ks, a further-but-richer representative costs nothing and the
budget floor can stay at 2.0. If it grows, the floor is buying bias, and the rung at which it
starts growing is where the floor belongs.

    python analyse_sister_rate_test.py --dir /group/esb/cesen/1kp/ks_analysis/1kp/sister_rate_test
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import calib_common as C

OUT_COLUMNS = ["branch_index", "focal_species", "species", "informal_name", "expected_ks",
               "expected_ks_to_outgroup", "bin_height", "lineage", "adjusted_mode",
               "adjusted_mode_sd", "ks_focal", "ks_sister", "original_mode",
               "deviation_from_branch_median", "relative_deviation"]


def read_design(directory):
    path = os.path.join(directory, "sister_rate_test_design.tsv")
    with open(path) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["role"] == "sister"]
    by_branch = {}
    for row in rows:
        by_branch.setdefault((row["branch_index"], row["focal_species"]), {})[row["informal_name"]] = row
    return by_branch


def read_adjustment(directory, branch_index, focal_informal):
    pattern = os.path.join(directory, "runs", f"b{branch_index}_{focal_informal}",
                           "rate_adjustment", focal_informal,
                           f"adjustment_table_{focal_informal}.tsv")
    found = glob.glob(pattern)
    if not found:
        return None
    with open(found[0]) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def spearman(xs, ys):
    """Rank correlation without scipy; ties get average ranks."""
    def ranks(values):
        order = sorted(range(len(values)), key=lambda i: values[i])
        out = [0.0] * len(values)
        i = 0
        while i < len(order):
            j = i
            while j < len(order) and values[order[j]] == values[order[i]]:
                j += 1
            average = (i + j + 1) / 2
            for k in order[i:j]:
                out[k] = average
            i = j
        return out

    if len(xs) < 3:
        return None
    rx, ry = ranks(xs), ranks(ys)
    mx, my = statistics.mean(rx), statistics.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return None if den == 0 else num / den


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=os.path.join(C.ROOT, "ks_analysis/1kp/sister_rate_test"))
    ap.add_argument("--out", help="write the per-sister table here (default: <dir>/sister_rate_test_results.tsv)")
    args = ap.parse_args()

    design = read_design(args.dir)
    out_rows, missing = [], []

    print(f"{'branch':>7}{'focal':>16}{'n':>4}{'median':>9}{'spread':>9}{'rel':>7}"
          f"{'reportedSD':>12}  per-sister adjusted mode")
    for (branch_index, focal_species), sisters in design.items():
        focal_informal = focal_species.split("_", 1)[0]
        table = read_adjustment(args.dir, branch_index, focal_informal)
        if table is None:
            missing.append(f"b{branch_index}_{focal_informal}")
            continue

        values, reported, per_sister = [], [], []
        for row in table:
            informal = row["Sister_Species"].split()[-1]
            if informal not in sisters:
                continue
            adjusted = float(row["Adjusted_Mode_Mean"])
            values.append(adjusted)
            reported.append(float(row["Adjusted_Mode_Mean_SD"]))
            per_sister.append((informal, adjusted, row))
        if not values:
            missing.append(f"b{branch_index}_{focal_informal} (no matching sisters)")
            continue

        median = statistics.median(values)
        spread = max(values) - min(values)
        shown = "  ".join(f"{name}:{value:.2f}" for name, value in
                          sorted(((n, v) for n, v, _ in per_sister), key=lambda t: t[1]))
        print(f"{branch_index:>7}{focal_species:>16}{len(values):>4}{median:>9.2f}{spread:>9.2f}"
              f"{100 * spread / median:>6.0f}%{statistics.mean(reported):>12.2f}  {shown}")

        for informal, adjusted, row in per_sister:
            d = sisters[informal]
            out_rows.append([branch_index, focal_species, d["species"], informal,
                             d["expected_ks"], d["expected_ks_to_outgroup"], d["bin_height"],
                             d["lineage"], f"{adjusted:.4f}", row["Adjusted_Mode_Mean_SD"],
                             row["Ks_Focal_Mean"], row["Ks_Sister_Mean"], row["Original_Mode"],
                             f"{adjusted - median:+.4f}",
                             f"{(adjusted - median) / median:+.4f}"])

    if missing:
        print("\nnot finished yet: " + ", ".join(missing))
    if not out_rows:
        return

    path = args.out or os.path.join(args.dir, "sister_rate_test_results.tsv")
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(OUT_COLUMNS)
        writer.writerows(out_rows)
    print(f"\nper-sister table: {path}")

    # What does the deviation track? Signed, because a one-sided bias is the thing that would
    # make a deeper representative dangerous; |deviation| alone would hide its direction.
    deviation = [float(r[13]) for r in out_rows]
    for label, column in (("sister expected Ks", 4), ("Ks(sister, outgroup)", 5),
                          ("sister bin height", 6)):
        values = [float(r[column]) for r in out_rows if r[column] != ""]
        if len(values) != len(deviation):
            continue
        rho = spearman(values, deviation)
        print(f"  Spearman(deviation, {label:<22}) = "
              + ("n/a" if rho is None else f"{rho:+.2f}"))
    print(f"  n = {len(out_rows)} sisters over {len(design)} branches")


if __name__ == "__main__":
    main()
