#!/usr/bin/env python3
"""Is there a better predictor of ortholog Ks than the concatenation path distance?

Uses the raw codeml output still on disk for the calibration pairs - every *.ks.tsv carries
Ka, Ks and Omega per gene pair - to ask three things:

  1. does a Ka ruler built from the SAME codeml run beat the concat path?
  2. does a per-species omega taken from the paranome runs transfer to the pair's ortholog
     omega, and does correcting by it tighten the calibration?
  3. is omega an independent lineage property, or just a restatement of depth?

    python3 alternative_predictors.py
"""
import os, sys, csv, glob, argparse, statistics, collections

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import calib_common as C

ORTHOLOG_GLOB = os.path.join(C.ROOT, "ks_analysis/**/ortholog_distributions/wgd_*_*")
MIN_PAIRS = 50
BANDS = [(0, 1), (1, 2), (2, 3), (3, 99)]


def band_label(lo, hi):
    return f"{lo}-{hi}" if hi < 99 else f"{lo}+"


def cv(values):
    mean = sum(values) / len(values)
    return 100 * (sum((x - mean) ** 2 for x in values) / len(values)) ** .5 / mean


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sx = sum((a - mx) ** 2 for a in x) ** .5
    sy = sum((b - my) ** 2 for b in y) ** .5
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / (sx * sy)


def medians(path):
    """Median Ka, Ks and Omega over the gene pairs of one wgd run."""
    ka, ks, omega = [], [], []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                a, s, o = float(row["Ka"]), float(row["Ks"]), float(row["Omega"])
            except (TypeError, ValueError, KeyError):
                continue
            if 0 < s < 10 and 0 < a < 10:
                ka.append(a); ks.append(s); omega.append(o)
    if len(ks) < MIN_PAIRS:
        return None
    return dict(n=len(ks), ka=statistics.median(ka), ks=statistics.median(ks),
                om=statistics.median(omega))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", default=C.CALIBRATION_PAIRS)
    args = ap.parse_args()

    meta = C.species_metadata()
    index = {}
    for folder in glob.glob(ORTHOLOG_GLOB, recursive=True):
        if "calibration_run" in folder and not os.path.exists(
                os.path.join(folder, os.path.basename(folder)[4:] + ".ks.tsv")):
            continue
        name = os.path.basename(folder)[4:]
        path = os.path.join(folder, name + ".ks.tsv")
        if os.path.exists(path) and os.path.getsize(path) > 1000:
            index.setdefault(tuple(sorted(name.split("_"), key=str.casefold)), path)
    print(f"on-disk ortholog ks.tsv files: {len(index)}")

    rows = []
    with open(args.pairs) as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue
            f = line.rstrip("\n").split("\t")
            t1 = C.resolve_tip(f[0].split()[-1], meta)
            t2 = C.resolve_tip(f[1].split()[-1], meta)
            if not t1 or not t2:
                continue
            key = tuple(sorted([meta[t1]["code"], meta[t2]["code"]], key=str.casefold))
            path = index.get(key)
            if not path:
                continue
            stats = medians(path)
            if stats:
                rows.append(dict(t1=t1, t2=t2, mode=float(f[2]), path=float(f[4]), **stats))
    print(f"calibration pairs with raw codeml output on disk: {len(rows)}\n")

    print("=== (1) how well does each ruler track the stored ortholog Ks mode? ===")
    print(f"{'band':<8}{'n':>4}{'CV of Ks/path':>18}{'CV of Ks/Ka(same run)':>26}")
    for lo, hi in BANDS:
        g = [r for r in rows if lo <= r["mode"] < hi]
        if len(g) < 3:
            continue
        print(f"{band_label(lo, hi):<8}{len(g):>4}"
              f"{cv([r['mode'] / r['path'] for r in g]):>17.1f}%"
              f"{cv([r['mode'] / r['ka'] for r in g]):>25.1f}%")

    print("\n=== (2) does a paranome omega stand in for the pair's ortholog omega? ===")
    needed = {t for r in rows for t in (r["t1"], r["t2"])}
    paranome = {}
    for tip in needed:
        path = os.path.join(C.PARALOG_DIR, f"wgd_{meta[tip]['code']}",
                            f"{meta[tip]['code']}.ks.tsv")
        if not os.path.exists(path):
            continue
        values = []
        with open(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                try:
                    s, o = float(row["Ks"]), float(row["Omega"])
                except (TypeError, ValueError, KeyError):
                    continue
                if 0.1 < s < 2 and 0 < o < 2:
                    values.append(o)
        if len(values) >= 100:
            paranome[tip] = statistics.median(values)
    usable = [r for r in rows if r["t1"] in paranome and r["t2"] in paranome]
    print(f"species with a usable paranome omega: {len(paranome)} of {len(needed)}; "
          f"pairs where both have one: {len(usable)}")
    if len(usable) >= 10:
        x = [(paranome[r["t1"]] + paranome[r["t2"]]) / 2 for r in usable]
        y = [r["om"] for r in usable]
        print(f" r(mean paranome omega, pair ortholog omega) = {pearson(x, y):+.2f}; "
              f"paranome omega is {statistics.median([a / b for a, b in zip(x, y)]):.2f}x "
              f"the ortholog omega")
        for lo, hi in BANDS:
            g = [r for r in usable if lo <= r["mode"] < hi]
            if len(g) < 5:
                continue
            raw = [r["mode"] / r["path"] for r in g]
            corrected = [r["mode"] / (r["path"] / ((paranome[r["t1"]] + paranome[r["t2"]]) / 2))
                         for r in g]
            print(f"  Ks {band_label(lo, hi):<4} (n={len(g):>3}): CV of Ks/path {cv(raw):>5.1f}% "
                  f"-> with paranome-omega correction {cv(corrected):>5.1f}%")

    print("\n=== (3) is omega independent of depth? ===")
    print(f"r(omega, path distance) over all {len(rows)} pairs = "
          f"{pearson([r['path'] for r in rows], [r['om'] for r in rows]):+.2f}")
    group = C.groups(meta)
    by = collections.defaultdict(list)
    for r in rows:
        by[tuple(sorted((group.get(r["t1"]), group.get(r["t2"]))))].append(r["om"])
    allom = [r["om"] for r in rows]

    def var(v):
        m = sum(v) / len(v)
        return sum((x - m) ** 2 for x in v) / len(v)

    sized = [v for v in by.values() if len(v) > 1]
    within = sum(var(v) * len(v) for v in sized) / sum(len(v) for v in sized)
    print(f"lineage pair explains {100 * (1 - within / var(allom)):.0f}% of the variance in omega "
          f"- but see the correlation above")
    print(f"\n{'lineage pair':<28}{'n':>4}{'median omega':>14}")
    for key, values in sorted(by.items(), key=lambda kv: -statistics.median(kv[1])):
        if len(values) > 1:
            print(f"{key[0] + ' / ' + key[1]:<28}{len(values):>4}{statistics.median(values):>14.3f}")


if __name__ == "__main__":
    main()
