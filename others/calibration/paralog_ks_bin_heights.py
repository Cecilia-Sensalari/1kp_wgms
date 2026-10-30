#!/usr/bin/env python3
"""Weighted-node bar heights of every species' paranome Ks distribution.

For each species, the paralog pairs of its paranome run are binned in intervals of 0.1 Ks
over the window Ks 0.5-2.5 and summed using the ksrates node weights, so that each
duplication node contributes a total weight of 1 regardless of how many pairs it produced.
The 20 resulting bar heights are what the median-bin-height indicator is computed from.

Writes ks_analysis/1kp/calibration/paralog_ks_bin_heights.tsv:
    species_tip, bin_0.5_0.6 ... bin_2.4_2.5, median, total, n_pairs_in_window

    python3 paralog_ks_bin_heights.py [--paralog-dir DIR] [--out FILE]
"""
import os, sys, csv, glob, argparse, statistics

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import calib_common as C

# ksrates stores one weight per paralog pair; the "outliers included" column is the one the
# plotted distribution uses, so it is the one the indicator has to be built from.
WEIGHT_COLUMN = "WeightOutliersIncluded"


def bars(path):
    out = [0.0] * C.N_BINS
    n = 0
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                ks = float(row["Ks"])
                weight = float(row[WEIGHT_COLUMN])
            except (TypeError, ValueError, KeyError):
                continue
            if C.BIN_LO <= ks < C.BIN_HI:
                out[int((ks - C.BIN_LO) / C.BIN_WIDTH)] += weight
                n += 1
    return out, n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paralog-dir", default=C.PARALOG_DIR)
    ap.add_argument("--out", default=C.BIN_HEIGHTS_TSV)
    args = ap.parse_args()

    meta = C.species_metadata()
    code_to_tip = {info["code"]: tip for tip, info in meta.items()}

    rows, missing = [], []
    for folder in sorted(glob.glob(os.path.join(args.paralog_dir, "wgd_*"))):
        code = os.path.basename(folder)[4:]
        path = os.path.join(folder, code + ".ks.tsv")
        if not os.path.exists(path) or os.path.getsize(path) < 1000:
            missing.append(code)
            continue
        values, n = bars(path)
        rows.append((code_to_tip.get(code, code), values, n))

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    edges = [C.BIN_LO + i * C.BIN_WIDTH for i in range(C.N_BINS)]
    with open(args.out, "w") as fh:
        fh.write("species_tip\t" + "\t".join(f"bin_{e:.1f}_{e + C.BIN_WIDTH:.1f}" for e in edges)
                 + "\tmedian\ttotal\tn_pairs_in_window\n")
        for tip, values, n in sorted(rows):
            fh.write(tip + "\t" + "\t".join(f"{v:.1f}" for v in values)
                     + f"\t{statistics.median(values):.1f}\t{sum(values):.1f}\t{n}\n")
    print(f"{len(rows)} species -> {args.out}")
    if missing:
        print(f"no usable paranome Ks file for {len(missing)} species: {missing[:10]}")


if __name__ == "__main__":
    main()
