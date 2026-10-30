#!/usr/bin/env python3
"""How well does each summary of the paranome bin heights predict a usable distribution?

Scores every candidate statistic of the twenty Ks 0.5-2.5 bars against the blind visual
labels in the filled review sheets (ks_analysis/1kp/ks_quality_review_sheet_*.tsv), using the
Mann-Whitney AUC with `usable = yes` as the positive class and everything else (no, doubtful,
borderline) as the negative one.

This is the measurement behind the choice of the median bar height as the A-species
indicator, and behind the answer to "why the median and not the total?".

    python3 indicator_auc.py
"""
import os, sys, csv, glob, argparse, statistics

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import calib_common as C

SHEETS = {"angiosperms": "ks_quality_review_sheet_angiosperms.tsv",
          "bryophytes": "ks_quality_review_sheet_bryophytes.tsv",
          "lycophytes": "ks_quality_review_sheet_lycophytes.tsv"}
SHEET_DIR = os.path.join(C.ROOT, "ks_analysis/1kp")

STATISTICS = {"median (the indicator)": statistics.median,
              "total (= 20 x mean)": sum,
              "mean": lambda v: sum(v) / len(v),
              "tallest bar": max,
              "lowest bar": min}


def auc(positive, negative):
    """P(a random positive scores above a random negative); ties count a half."""
    if not positive or not negative:
        return float("nan")
    hits = sum(1.0 if p > q else 0.5 if p == q else 0.0 for p in positive for q in negative)
    return hits / (len(positive) * len(negative))


def labels(path):
    """code -> True if judged usable. Tolerates the typos in the filled sheets."""
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            verdict = (row.get("usable") or "").strip().lower()
            if verdict:
                out[row["code"].strip()] = verdict.startswith("y")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin-heights", default=C.BIN_HEIGHTS_TSV)
    args = ap.parse_args()

    bars = C.bin_heights(args.bin_heights)
    meta = C.species_metadata()
    code_to_tip = {info["code"]: tip for tip, info in meta.items()}

    # the review sheets were drawn up per lineage; bryophytes and lycophytes were one round
    sets = {"angiosperms": labels(os.path.join(SHEET_DIR, SHEETS["angiosperms"]))}
    combined = labels(os.path.join(SHEET_DIR, SHEETS["bryophytes"]))
    combined.update(labels(os.path.join(SHEET_DIR, SHEETS["lycophytes"])))
    sets["bryo+lyco"] = combined

    for name, labelled in sets.items():
        usable = sum(1 for v in labelled.values() if v)
        print(f"{name}: {len(labelled)} labelled species ({usable} usable, "
              f"{len(labelled) - usable} not)")
    print(f"total labelled: {sum(len(v) for v in sets.values())}\n")

    print(f"{'statistic':<24}" + "".join(f"{k:>14}" for k in sets) + f"{'pooled':>10}")
    for name, function in STATISTICS.items():
        pooled_pos, pooled_neg = [], []
        row = f"{name:<24}"
        for labelled in sets.values():
            positive, negative = [], []
            for code, ok in labelled.items():
                tip = code_to_tip.get(code, code)
                if tip not in bars:
                    continue
                (positive if ok else negative).append(function(bars[tip]))
            pooled_pos += positive
            pooled_neg += negative
            row += f"{auc(positive, negative):>14.3f}"
        print(row + f"{auc(pooled_pos, pooled_neg):>10.3f}")


if __name__ == "__main__":
    main()
