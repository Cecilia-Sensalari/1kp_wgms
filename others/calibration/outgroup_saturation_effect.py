#!/usr/bin/env python3
"""Does a saturated subtracted term actually damage the rate adjustment?

The selection rule assumes it does: a candidate is rejected, or flagged, when the ortholog Ks
peak between it and an outgroup would sit beyond the measurable range, because that peak enters

    Ks_adjusted = ( Ks(A,B) + Ks(A,O) - Ks(B,O) ) / 2

with a minus sign. `check_outgroup_term_flags.py` shows the rule predicts WHICH peaks will be
unmeasurable. This script asks the next question: when one of them is used anyway, is the answer
worse?

WHY THIS DESIGN. The obvious test - recompute a branch with the species the old rule chose and
compare - does not work, because the comparison then differs in the species as well as in the
saturation. Two representatives that happen to be near-identical lineages agree almost by
construction; two that are far apart disagree for reasons that have nothing to do with the
subtracted term. Both confounds were met in practice on branches 668 and 669.

This test avoids them by never changing a species. ksrates dates each divergence once per
available outgroup and writes every one of them to adjustment_table_<focal>_all.tsv. Within one
such group the focal species, the sister and all the input data are fixed; the only thing that
varies between rows is which outgroup was subtracted. So:

    group   = (run, node, focal species, sister)   - several outgroups, one divergence
    split   = outgroups whose Ks(sister, outgroup) is measurable, versus those where it is not
    measure = how far each outgroup's estimate falls from the consensus of the measurable ones

The measurable rows are scored LEAVE-ONE-OUT, against the median of the *other* measurable rows,
so that they are not being compared with a median they themselves define.

Rows inside a group are not independent, so the headline test is paired: one comparison per
group, Wilcoxon signed-rank over the per-group differences. The pooled row counts are reported
too, but they overstate the evidence.

    python3 outgroup_saturation_effect.py [--runs <datasets dir> ...] [--limit 2.72] [--examples 3]

Needs an ortholog peak export (`ksrates inspect-ks-db <address file> --type ortholog`) to decide
which sister-outgroup peaks are saturated, and the _all.tsv adjustment tables of finished runs.
Runs nothing itself.

ON THE FERN SHARD (10 branches, 2026-11) this gave 22 usable groups: median deviation 3.2% for
measurable outgroups against 10.5% for saturated ones, but only 14 of 22 groups in that direction
and a Wilcoxon p of 0.53 - the right sign, nowhere near conclusive. 166 fern groups have three or
more outgroups and only 22 of those are mixed, so the full 2354-branch run should yield several
hundred comparisons rather than twenty.
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


def load_peaks(path):
    """{(code1, code2) sorted: mode}. The export carries an unnamed index column first."""
    peaks = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                peaks[tuple(sorted((row["Species1"].split()[-1],
                                    row["Species2"].split()[-1])))] = float(row["Mode"])
            except (TypeError, ValueError, KeyError):
                continue
    return peaks


def collect(run_dirs):
    """{(table, node, focal, sister): [(outgroup, adjusted Ks), ...]} from the _all tables."""
    groups = {}
    for root in run_dirs:
        pattern = os.path.join(root, "*", "*", "rate_adjustment", "*",
                               "adjustment_table_*_all.tsv")
        for path in glob.glob(pattern):
            with open(path) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    try:
                        key = (path, row["Node"], row["Focal_Species"], row["Sister_Species"])
                        groups.setdefault(key, []).append(
                            (row["Out_Species"], float(row["Adjusted_Mode"])))
                    except (TypeError, ValueError, KeyError):
                        continue
    return groups


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", nargs="+",
                    default=[os.path.join(C.ROOT, "ks_analysis/1kp/datasets")],
                    help="dataset directories of finished runs")
    ap.add_argument("--peaks", help="ortholog peak export (default: the newest in ks_database/)")
    ap.add_argument("--limit", type=float, default=None,
                    help="measured Ks above which the subtracted term counts as saturated "
                         "(default: the Ks of RELIABLE_PATH)")
    ap.add_argument("--examples", type=int, default=3,
                    help="how many of the clearest groups to print in full")
    args = ap.parse_args()

    peaks_path = args.peaks
    if not peaks_path:
        found = sorted(glob.glob(os.path.join(C.ROOT, "ks_analysis/1kp/ks_database",
                                              "*ortholog_peaks_*.tsv")))
        if not found:
            sys.exit("no ortholog peak export found; pass --peaks")
        peaks_path = found[-1]
    limit = args.limit if args.limit is not None else C.path_to_ks(C.RELIABLE_PATH)

    peaks = load_peaks(peaks_path)
    groups = collect(args.runs)
    print(f"{len(peaks)} ortholog peaks from {os.path.basename(peaks_path)}")
    print(f"{len(groups)} (focal, sister) groups; a subtracted term is saturated above "
          f"Ks {limit:.2f}\n")

    pooled_ok, pooled_bad, pairs, shown = [], [], [], []
    missing_peaks, skipped_for_missing = set(), 0
    for (path, node, focal, sister), rows in groups.items():
        labelled = []
        for outgroup, value in rows:
            key = tuple(sorted((sister.split()[-1], outgroup.split()[-1])))
            peak = peaks.get(key)
            if peak is None:
                # Almost always a stale export rather than a genuinely absent pair: the peak is
                # written to the Ks database server, and this file is only a snapshot of it.
                missing_peaks.add(key)
                labelled = None
                skipped_for_missing += 1
                break
            labelled.append((outgroup, value, peak))
        if not labelled or len(labelled) < 3:
            continue
        good = [r for r in labelled if r[2] <= limit]
        bad = [r for r in labelled if r[2] > limit]
        if len(good) < 2 or not bad:
            continue                      # need a consensus to judge against, and something to judge

        values = [r[1] for r in good]
        reference = statistics.median(values)
        if reference <= 0:
            continue
        loo = []                          # leave-one-out, so the median is not self-defined
        for i, value in enumerate(values):
            other = statistics.median(values[:i] + values[i + 1:])
            if other > 0:
                loo.append(abs(value - other) / other * 100)
        dev_bad = [abs(r[1] - reference) / reference * 100 for r in bad]
        if not loo or not dev_bad:
            continue

        pooled_ok += loo
        pooled_bad += dev_bad
        a, b = statistics.median(loo), statistics.median(dev_bad)
        pairs.append((a, b))
        shown.append((b - a, focal, sister, labelled, reference))

    if missing_peaks:
        print(f"WARNING: {skipped_for_missing} groups were skipped because {len(missing_peaks)} "
              f"sister-outgroup pairs are not in this peak export, e.g. "
              f"{', '.join('-'.join(p) for p in sorted(missing_peaks)[:3])}.")
        print("         The export is a snapshot; refresh it with")
        print("           ksrates inspect-ks-db <address file> --type ortholog")
        print("         and run again, or pass a newer file with --peaks.\n")

    if not pairs:
        print("no group had both kinds of outgroup; nothing to compare"
              + (" - but see the warning above, which is the more likely explanation"
                 if missing_peaks else ""))
        return

    print(f"{len(pairs)} groups have both kinds of outgroup and at least two measurable ones\n")
    print(f"  {'subtracted term':26}{'n':>5}{'median':>9}{'mean':>8}{'p90':>8}")
    for label, data in (("measurable", pooled_ok), ("saturated", pooled_bad)):
        data = sorted(data)
        print(f"  {label:26}{len(data):>5}{statistics.median(data):>8.1f}%"
              f"{statistics.mean(data):>7.1f}%{data[int(0.9 * len(data))]:>7.1f}%")
    print("  (pooled rows are not independent - the paired test below is the real one)")

    worse = sum(1 for a, b in pairs if b > a)
    print(f"\npaired, one comparison per group:")
    print(f"  saturated outgroup further from the consensus: {worse} of {len(pairs)} groups")
    print(f"  median per group: measurable {statistics.median(a for a, b in pairs):.1f}%, "
          f"saturated {statistics.median(b for a, b in pairs):.1f}%")
    try:
        from scipy.stats import wilcoxon
        print(f"  Wilcoxon signed-rank: p = {wilcoxon([b - a for a, b in pairs]).pvalue:.4f}")
    except Exception as err:
        print(f"  (scipy unavailable, no signed-rank test: {err})")

    for _, focal, sister, labelled, reference in sorted(shown, reverse=True)[:args.examples]:
        print(f"\nclearest group - focal {focal}, sister {sister}")
        print(f"  {'outgroup':32}{'Ks(sister,outgroup)':>20}{'':2}{'dated at':>9}{'off by':>9}")
        for outgroup, value, peak in sorted(labelled, key=lambda r: r[2]):
            tag = "saturated " if peak > limit else "measurable"
            print(f"  {outgroup:32}{peak:>15.2f} {tag}{value:>9.3f}"
                  f"{abs(value - reference) / reference * 100:>8.1f}%")
        print(f"  consensus of the measurable outgroups: {reference:.3f}")


if __name__ == "__main__":
    main()
