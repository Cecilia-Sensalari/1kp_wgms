#!/usr/bin/env python3
"""Build the ortholog-Ks / path-distance calibration table.

Joins the ortholog peak exports of the project's two Ks databases (written by
`ksrates inspect-ks-db <address file> --type ortholog`) to the tip-to-tip path distance on
the concatenation ML tree, and writes one row per unique species pair.

Writes ks_analysis/1kp/calibration/calibration_pairs.tsv:
    sp1, sp2, ortholog_ks_mode, mode_sd, path_distance, source

    python3 build_calibration_pairs.py
"""
import os, sys, glob, argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import calib_common as C

DATABASES = [("1kp", os.path.join(C.ROOT, "ks_analysis/1kp/ks_database")),
             ("master1", os.path.join(C.ROOT, "ks_analysis/master_1_2026/ks_analysis/ks_database"))]


def newest_export(folder):
    """Most recent ortholog peak export; the filenames carry a sortable YYYYMMDD_HHMMSS stamp."""
    found = sorted(glob.glob(os.path.join(folder, "ks_database_server_address_ortholog_peaks_*.tsv")))
    if not found:
        raise SystemExit(f"no ortholog peaks export in {folder} - run ksrates inspect-ks-db first")
    return found[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=C.CALIBRATION_PAIRS)
    args = ap.parse_args()

    _, dist = C.load_concat_tree()
    meta = C.species_metadata()
    print(f"concat tree: {len(dist.index)} tips", file=sys.stderr)

    pairs, unmapped = {}, set()
    for source, folder in DATABASES:
        path = newest_export(folder)
        print(f"using {source}: {os.path.basename(path)}", file=sys.stderr)
        with open(path) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 5:
                    continue
                name1, name2, mode, sd = f[1], f[2], f[3], f[4]
                tip1 = C.resolve_tip(name1.split()[-1], meta)
                tip2 = C.resolve_tip(name2.split()[-1], meta)
                if tip1 is None or tip2 is None:
                    unmapped.add(name1 if tip1 is None else name2)
                    continue
                d = dist(tip1, tip2)
                if d is None:
                    continue
                pairs.setdefault(tuple(sorted([tip1, tip2])),
                                 (float(mode), float(sd), d, source, name1, name2))

    if unmapped:
        print(f"unmapped species: {sorted(unmapped)}", file=sys.stderr)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("sp1\tsp2\tortholog_ks_mode\tmode_sd\tpath_distance\tsource\n")
        for _, (mode, sd, d, source, name1, name2) in sorted(pairs.items(), key=lambda kv: kv[1][0]):
            fh.write(f"{name1}\t{name2}\t{mode:.4f}\t{sd:.4f}\t{d:.6f}\t{source}\n")
    print(f"{len(pairs)} unique pairs -> {args.out}")


if __name__ == "__main__":
    main()
