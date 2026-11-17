#!/usr/bin/env python3
"""Weighted-node bar heights of every species' paranome Ks distribution.

For each species, the paralog pairs of its paranome run are binned in intervals of 0.1 Ks
over the window Ks 0.5-2.5 and summed using node weights, so that each duplication node
contributes a total weight of 1 regardless of how many pairs it produced. The weights are
recomputed (calib_common.weighted_paralog_ks) over the Ks range of the plots, [0.005, --max-ks],
so the bars are the ones ksrates draws; the weight columns stored in the .ks.tsv are not used.
The 20 resulting bar heights are what the median-bin-height indicator is computed from.

Writes ks_analysis/1kp/calibration/paralog_ks_bin_heights.tsv:
    species_tip, bin_0.5_0.6 ... bin_2.4_2.5, median, total, n_pairs_in_window

    python3 paralog_ks_bin_heights.py [--paralog-dir DIR] [--out FILE] [--max-ks 5]
"""
import os, sys, glob, argparse, statistics

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import calib_common as C


def bars(path, max_ks):
    out = [0.0] * C.N_BINS
    n = 0
    for _, _, ks, weight in C.weighted_paralog_ks(path, max_ks=max_ks):
        if C.BIN_LO <= ks < C.BIN_HI:
            out[int((ks - C.BIN_LO) / C.BIN_WIDTH)] += weight
            n += 1
    return out, n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paralog-dir", default=C.PARALOG_DIR)
    ap.add_argument("--out", default=C.BIN_HEIGHTS_TSV)
    ap.add_argument("--max-ks", type=float, default=C.PARALOG_MAX_KS,
                    help="max Ks of the plots; the node weights sum to 1 within [0.005, max Ks]")
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
        values, n = bars(path, args.max_ks)
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
