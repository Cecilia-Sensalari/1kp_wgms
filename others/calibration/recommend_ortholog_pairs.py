#!/usr/bin/env python3
"""Choose the species pairs worth running to extend the calibration.

Partitions the tree into lineage buckets, takes the best few species of each ranked by
published genome first and then BUSCO completeness (an ortholog Ks estimate needs a complete
gene set, not many paralogs), and selects pairs evenly in PATH DISTANCE - the predictor -
with every depth band covered by several independent lineage combinations.

Writes ks_analysis/1kp/calibration/recommended_ortholog_pairs.tsv.

    python3 recommend_ortholog_pairs.py [--per-band 4] [--min-path 0.45] [--max-path 2.35]
"""
import os, sys, argparse, collections

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import calib_common as C

PER_BUCKET = 4
MAX_PER_SPECIES = 7
MAX_PER_BUCKET_PAIR = 2
HIGH_PRIORITY_PATH = 0.95      # below this the calibration is already well populated

GENOME_BUCKET = {
    "Arath_TAIR10": "eudicot-rosid", "Orysa_v7.0": "monocot", "Vitvi_Genoscope.12X": "eudicot-rosid",
    "Solly_iTAGv2.3": "eudicot-asterid", "Ambtr_v1.0.27": "angiosperm-basal",
    "Nelnu_v1.0": "eudicot-basal", "Aquco_v1.1": "eudicot-basal", "Musac_v1.0": "monocot",
    "Sorbi_v2.1": "monocot", "Poptr_v3.0": "eudicot-rosid", "Pinta_v2.0": "gymnosperm-conifer",
    "Selmo_v1.0": "lycophyte", "Phypa_v3.0": "moss", "Klefl_v1.0": "alga-strep",
    "Chlre_v5.5": "alga-chloro", "Mimgu_v2.0": "eudicot-asterid", "Betvu_v1.1": "eudicot-caryo",
    "Theca_v1.1": "eudicot-rosid", "Spipo_v2": "monocot", "Elagu_v2.0": "monocot",
}
_KEYS = (("Lycopodiopsida", "lycophyte"), ("Polypodiopsida", "fern"),
         ("Pinopsida", "gymnosperm-conifer"), ("Cupressaceae", "gymnosperm-conifer"),
         ("Pinaceae", "gymnosperm-conifer"), ("Araucaria", "gymnosperm-conifer"),
         ("Podocarp", "gymnosperm-conifer"), ("Acrogymnospermae", "gymnosperm-other"),
         ("Bryophyta", "moss"), ("Marchantiophyta", "liverwort"),
         ("Anthocerotophyta", "hornwort"), ("Liliopsida", "monocot"),
         ("Caryophyllales", "eudicot-caryo"), ("asterids", "eudicot-asterid"),
         ("rosids", "eudicot-rosid"), ("Ranunculales", "eudicot-basal"),
         ("Proteales", "eudicot-basal"), ("eudicotyledons", "eudicot-rosid"),
         ("Magnoliopsida", "angiosperm-basal"),
         ("Charophyceae", "alga-strep"), ("Zygnemato", "alga-strep"),
         ("Klebsormidio", "alga-strep"), ("Coleochaet", "alga-strep"),
         ("Mesostigma", "alga-strep"))


def bucket(tip, meta, lin):
    if tip in GENOME_BUCKET:
        return GENOME_BUCKET[tip]
    lineage = lin.get(meta[tip]["code"], "")
    return next((b for key, b in _KEYS if key in lineage), "alga-chloro")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-band", type=int, default=4)
    ap.add_argument("--min-path", type=float, default=0.45)
    ap.add_argument("--max-path", type=float, default=2.35)
    ap.add_argument("--exclude", default="", help="comma-separated tree tips to keep out of "
                    "the pool (e.g. species whose all-vs-all blastp is pathologically slow)")
    ap.add_argument("--out", default=os.path.join(C.DATA_DIR, "recommended_ortholog_pairs.tsv"))
    args = ap.parse_args()

    _, dist = C.load_concat_tree()
    meta = C.species_metadata()
    lin = C.lineages()
    busco = C.busco()
    median_bin = C.median_bin_heights()

    excluded = {x.strip() for x in args.exclude.split(",") if x.strip()}
    if excluded:
        print(f"excluded from the pool: {sorted(excluded)}")
    pool = collections.defaultdict(list)
    for tip in dist.index:
        if tip not in meta or tip in excluded:
            continue
        is_genome = tip in GENOME_BUCKET
        pool[bucket(tip, meta, lin)].append(
            (tip, median_bin.get(tip, 0), busco.get(meta[tip]["code"], 100.0 if is_genome else 0.0),
             is_genome))
    chosen = {}
    for name, entries in pool.items():
        entries.sort(key=lambda x: (-x[3], -x[2], -x[1]))
        chosen[name] = entries[:PER_BUCKET]
    print("pool (<=%d best per lineage, ranked by genome > BUSCO > paranome depth):" % PER_BUCKET)
    for name in sorted(chosen):
        print(f"  {name:<20}" + ", ".join(
            f"{t}{'(G)' if g else ''} busco {b:.0f} bin {q:.0f}" for t, q, b, g in chosen[name]))

    already = set()
    if os.path.exists(C.CALIBRATION_PAIRS):
        with open(C.CALIBRATION_PAIRS) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    continue
                f = line.split("\t")
                a = C.resolve_tip(f[0].split()[-1], meta)
                b = C.resolve_tip(f[1].split()[-1], meta)
                if a and b:
                    already.add(tuple(sorted([a, b])))
    print(f"\npairs already in the calibration: {len(already)}")

    candidates = []
    names = sorted(chosen)
    for i, first in enumerate(names):
        for second in names[i:]:
            for t1, q1, _, g1 in chosen[first]:
                for t2, q2, _, g2 in chosen[second]:
                    if t1 >= t2 or tuple(sorted([t1, t2])) in already:
                        continue
                    d = dist(t1, t2)
                    if d is None or not (args.min_path <= d <= args.max_path):
                        continue
                    candidates.append((d, t1, t2, first, second, min(q1, q2), g1 + g2))
    print(f"candidate pairs in path {args.min_path}-{args.max_path}: {len(candidates)}")

    used, bucket_pairs, selected = collections.Counter(), collections.Counter(), []
    step = 0.1
    bands = [(round(args.min_path + i * step, 2), round(args.min_path + (i + 1) * step, 2))
             for i in range(int(round((args.max_path - args.min_path) / step)))]
    for low, high in bands:
        in_band = sorted([c for c in candidates if low <= c[0] < high],
                         key=lambda c: (bucket_pairs[tuple(sorted((c[3], c[4])))],
                                        used[c[1]] + used[c[2]], -c[6], -c[5]))
        taken = 0
        for c in in_band:
            if taken >= args.per_band:
                break
            if used[c[1]] >= MAX_PER_SPECIES or used[c[2]] >= MAX_PER_SPECIES:
                continue
            key = tuple(sorted((c[3], c[4])))
            if bucket_pairs[key] >= MAX_PER_BUCKET_PAIR:
                continue
            selected.append(c)
            used[c[1]] += 1
            used[c[2]] += 1
            bucket_pairs[key] += 1
            taken += 1
    selected.sort()

    with open(args.out, "w") as fh:
        fh.write("species_1\tspecies_2\tlineage_1\tlineage_2\tpath_distance\t"
                 "predicted_ks_linear185\tmin_median_bin_height\tn_genomes_in_pair\tpriority\n")
        for d, t1, t2, b1, b2, q, g in selected:
            fh.write(f"{t1}\t{t2}\t{b1}\t{b2}\t{d:.4f}\t{1.85 * d:.2f}\t{q:.0f}\t{g}\t"
                     f"{'high' if d >= HIGH_PRIORITY_PATH else 'low'}\n")
    print(f"\nselected {len(selected)} pairs -> {args.out}")
    print(f"{'path band':<12}{'Ks (1.85)':>13}{'new':>5}   lineage pairs added")
    for low, high in bands:
        group = [c for c in selected if low <= c[0] < high]
        if group:
            print(f"{f'{low:.2f}-{high:.2f}':<12}{f'{1.85 * low:.2f}-{1.85 * high:.2f}':>13}"
                  f"{len(group):>5}   " + ", ".join(sorted({f'{c[3]}/{c[4]}' for c in group})))
    print(f"\nspecies used: {len(used)}; busiest: {used.most_common(5)}")


if __name__ == "__main__":
    main()
