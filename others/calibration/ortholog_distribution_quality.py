#!/usr/bin/env python3
"""Score how measurable the peak of an ortholog Ks distribution actually is.

The expected-Ks calibration predicts, from the concatenation-tree path distance, how deep an
A-B split is. This script answers the complementary question AFTER a run: given the ortholog
Ks distribution that came out, is its peak a measurement or an artefact of the
max_ks_orthologs cutoff?

The statistic is the CONCENTRATION: the fraction of retained Ks values lying within
+-HALF_WIDTH of the distribution mode. It is read against the CHANCE LEVEL, the value a flat
distribution of the same width would give: the width is the 5-95 percentile span of the
retained values, and a 2*HALF_WIDTH-wide band placed anywhere on a flat distribution of that
span holds 2*HALF_WIDTH/span of it. The chance level is not a constant - a shallow pair spans
~3 Ks and a deep one ~7 - so it is computed per distribution and reported beside the
concentration. Over the 86 distributions measured when this was written:

    path    0-0.6  0.6-0.8  0.8-1.0  1.0-1.2  1.2-1.5  1.5-1.8  1.8-2.1  2.1-2.4  2.4-3.6
    conc      68%      55%      41%      21%      20%      17%      16%      16%      16%
    chance    32%      22%      17%      15%      15%      15%      15%      15%      15%
    ratio    2.1x     2.5x     2.4x    1.37x    1.33x    1.18x    1.09x    1.06x    1.08x

The mode holds 2.1-2.5x chance up to path 1.0 and 1.37x straight after it, so the thresholds
below split "the peak is a measurement" from "the peak is where the KDE happened to land on a
flat smear". An earlier version of this script compared every concentration against a single
2*HALF_WIDTH/max_ks = 10% floor, the value for a distribution filling the whole 0-10 window.
No real distribution is that wide, so that floor was too low and flattered the deep pairs.

Unlike the path-based flag this is measured on the run's own output, so it is the one to quote
in the results.

    python3 ortholog_distribution_quality.py <run directory> [...]

Each run directory is one holding an `ortholog_distributions/wgd_<a>_<b>/` tree.
"""
import os, sys, glob, argparse
from statistics import median
import numpy as np
from scipy.stats import gaussian_kde

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import calib_common as C

HALF_WIDTH = 0.5
GOOD, DEGRADED = 0.40, 0.20          # concentration thresholds; see the table above
GRID_POINTS = 2000


def concentration(ks_values, max_ks):
    """Mode, concentration, chance level for this distribution's own width, and n."""
    kept = np.asarray([k for k in ks_values if 0 < k <= max_ks])
    if len(kept) < 100:
        return None
    grid = np.linspace(0.01, max_ks, GRID_POINTS)
    mode = grid[np.argmax(gaussian_kde(kept)(grid))]
    conc = float(((kept > mode - HALF_WIDTH) & (kept < mode + HALF_WIDTH)).mean())
    lo, hi = np.percentile(kept, [5, 95])
    chance = float(2 * HALF_WIDTH / (hi - lo))
    return mode, conc, chance, len(kept)


def verdict(conc):
    return "measured" if conc >= GOOD else "degraded" if conc >= DEGRADED else "no peak"


def read_ks(path):
    out = []
    with open(path) as fh:
        column = fh.readline().rstrip("\n").split("\t").index("Ks")
        for line in fh:
            try:
                out.append(float(line.rstrip("\n").split("\t")[column]))
            except (ValueError, IndexError):
                pass
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dirs", nargs="+")
    ap.add_argument("--max-ks", type=float, default=10.0,
                    help="max_ks_orthologs of the run that produced these (default 10)")
    args = ap.parse_args()

    _, dist = C.load_concat_tree()
    meta = C.species_metadata()

    print(f"{'pair':<18}{'path':>6}{'flag':>12}{'n':>7}{'mode':>8}{'conc%':>8}"
          f"{'chance%':>9}{'ratio':>7}{'verdict':>10}")
    rows = []
    for run in args.run_dirs:
        base = os.path.join(run, "ortholog_distributions")
        base = base if os.path.isdir(base) else run
        for folder in sorted(glob.glob(os.path.join(base, "wgd_*"))):
            pair = os.path.basename(folder)[4:]
            ks_file = os.path.join(folder, pair + ".ks.tsv")
            if not os.path.exists(ks_file):
                continue
            result = concentration(read_ks(ks_file), args.max_ks)
            if result is None:
                continue
            mode, conc, chance, n = result
            a, b = pair.split("_", 1)
            tip_a, tip_b = C.resolve_tip(a, meta), C.resolve_tip(b, meta)
            path = dist(tip_a, tip_b) if tip_a and tip_b else None
            rows.append((pair, path, mode, conc, chance, n))

    for pair, path, mode, conc, chance, n in sorted(rows, key=lambda r: (r[1] is None, r[1])):
        flag = C.path_reliability(path) or "no path"
        shown_path = f"{path:>6.2f}" if path is not None else f"{'NA':>6}"
        print(f"{pair:<18}{shown_path}{flag:>12}{n:>7}{mode:>8.2f}"
              f"{100 * conc:>8.1f}{100 * chance:>9.1f}{conc / chance:>7.2f}"
              f"{verdict(conc):>10}")

    agree = sum(1 for _, p, _, c, _, _ in rows
                if p is not None and (C.path_reliability(p) == "reliable") == (c >= GOOD))
    with_path = sum(1 for _, p, _, _, _, _ in rows if p is not None)
    print(f"\n{len(rows)} distributions; median chance level = "
          f"{100 * median(r[4] for r in rows):.0f}%")
    if with_path:
        print(f"path flag and measured verdict agree on {agree}/{with_path} "
              f"({100 * agree / with_path:.0f}%)")


if __name__ == "__main__":
    main()
