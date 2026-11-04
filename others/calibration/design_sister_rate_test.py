#!/usr/bin/env python3
"""Design the sister-rate test: does a further-but-richer B beat a nearer-but-poorer one?

Every species in a branch's sister clade B split from A at the SAME node, so they all estimate
the same divergence. What varies between them is not age but substitution rate - which is
exactly what the expected-Ks budget reads when it ranks candidates within a clade. So a ladder
of candidates from one sister clade, spanning expected Ks from ~1.0 to ~2.7, is a controlled
rate series with the divergence held fixed, and the rate adjustment is supposed to return the
same `Adjusted_Mode_Mean` for every rung.

The spread across the ladder is the error the budget rule is trading against distribution
quality. If it is flat, a further-but-richer representative is free and the budget floor can
stay at 2.0 (or rise); if it grows with expected Ks, the floor is buying bias and should come
down.

This generalises ks_analysis/1kp/rate_test (one branch, six seed-plant sisters, PKOX focal),
whose sisters sat at Ks ~3.6-4.6 - already saturated - and whose fast/slow contrast was
confounded with conifer-vs-angiosperm. Here the ladders sit inside the calibrated range and
each ladder is a single sister clade.

Writes, per selected branch, a ksrates config carrying the topology induced from the analysis
tree. It has to be strictly bifurcating - ksrates rejects a polytomy outright - so a flat sister
clade is not an option, and the nesting inside the clade costs a few ortholog pairs beyond the
trios (ksrates computes them for its distance tree: rate_test used 18 for 13 trio pairs). The
nesting does not affect the comparison: the divergence node is the focal species' sister clade
as a whole, so every sister still sits on the same node.

    python design_sister_rate_test.py --out-dir /group/esb/cesen/1kp/ks_analysis/1kp/sister_rate_test
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import calib_common as C

# One candidate per rung, so the ladder spans the whole range the budget floor can reach:
# 2.0 is the floor itself, 2.7 the edge of the calibrated range (path 1.0).
RUNGS = (1.0, 1.35, 1.7, 2.0, 2.35, 2.7)
RUNG_TOLERANCE = 0.18
MIN_RUNGS = 5
MIN_FOCAL_BIN_HEIGHT = 60
MIN_OUTGROUP_BIN_HEIGHT = 60

# Chosen with --survey from the branches that admit a ladder, and ranked by how saturated the
# SUBTRACTED term Ks(sister, outgroup) is, because that term, not the sister's own depth, is
# what the adjustment can least afford to lose. It is not simply the deepest term: the outgroup
# is older in time but can be much NEARER in substitutions than the sisters are to the focal,
# which is what happens whenever the focal itself is on a fast branch. Five of these six come
# out comfortably measurable on that term (0.65-2.38); the fern branch 667 is the one where it
# saturates (1.68-3.29), and it is kept deliberately as the contrasting arm.
DEFAULT_BRANCHES = (1323, 2221, 1096, 1322, 1883, 667)

TEMPLATE = """[SPECIES]
focal_species = {focal}
newick_tree = {newick}

latin_names = {latin}

fasta_filenames = {fasta}

ks_database_server_address_path = {db_address}


[ANALYSIS SETTING]
paranome = yes
collinearity = no
reciprocal_retention = no

gff_feature = mrna
gff_attribute = id

max_number_outgroups = 4
consensus_mode_for_multiple_outgroups = mean among outgroups


[PARAMETERS]
x_axis_max_limit_paralogs_plot = 5
bin_width_paralogs = 0.1
y_axis_max_limit_paralogs_plot = None

num_bootstrap_iterations = 200
divergence_colors = LightCoral, Red, MediumBlue, #A0720A
x_axis_max_limit_orthologs_plots = 5
bin_width_orthologs = 0.1

max_ks_paralogs = 5
max_ks_orthologs = 10
"""


def short(name: str) -> str:
    """ksrates informal names take no underscore: Nelnu_v1.0 -> Nelnu (as the batch wrapper does)."""
    return name.split("_", 1)[0]


def read_metadata(path):
    with open(path) as fh:
        return {row["species"]: row for row in csv.DictReader(fh, delimiter="\t")}


def prune(node, keep):
    """Induced topology of `keep` as nested tuples, collapsing nodes left with one child."""
    if not node.children:
        return node.name if node.name in keep else None
    kids = [k for k in (prune(c, keep) for c in node.children) if k is not None]
    if not kids:
        return None
    return kids[0] if len(kids) == 1 else tuple(kids)


def subtree_leaves(tree):
    if isinstance(tree, str):
        return [tree]
    return [leaf for child in tree for leaf in subtree_leaves(child)]


def resolve_polytomies(tree, dist):
    """Make every node bifurcating, joining the closest two children first.

    ksrates refuses a tree with unresolved relationships, and the 1KP analysis tree does
    contain polytomies. Where one survives into a subtree, the children are joined in order of
    their shortest path distance, so the invented structure at least follows substitution
    distance rather than the order the leaves happened to be written in. It changes nothing in
    the comparison: the divergence node is the focal species' sister clade as a whole, so every
    sister sits on the same node whatever the shape inside the clade.
    """
    if isinstance(tree, str):
        return tree
    kids = [resolve_polytomies(child, dist) for child in tree]
    while len(kids) > 2:
        best = None
        for i in range(len(kids)):
            for j in range(i + 1, len(kids)):
                gap = min((dist(a, b) for a in subtree_leaves(kids[i])
                           for b in subtree_leaves(kids[j])
                           if dist(a, b) is not None), default=None)
                if gap is not None and (best is None or gap < best[0]):
                    best = (gap, i, j)
        if best is None:
            best = (0.0, 0, 1)
        _, i, j = best
        joined = (kids[i], kids[j])
        kids = [k for n, k in enumerate(kids) if n not in (i, j)] + [joined]
    return tuple(kids)


def render(tree, names):
    if isinstance(tree, str):
        return names[tree]
    return "(" + ",".join(render(child, names) for child in tree) + ")"


def ladder(candidates):
    """The richest candidate within RUNG_TOLERANCE of each rung, each species used once."""
    picked, used = [], set()
    for rung in RUNGS:
        near = [c for c in candidates if abs(c["expected_ks"] - rung) <= RUNG_TOLERANCE
                and c["species"] not in used]
        if near:
            best = max(near, key=lambda c: c["bin_height"])
            used.add(best["species"])
            picked.append((rung, best))
    return picked


def branch_design(index, node, dist, bins, lineages):
    a_species = C.pick_A(node, bins)
    if not a_species:
        return None
    focal = max(a_species, key=lambda s: bins.get(s, 0))
    if bins.get(focal, 0) < MIN_FOCAL_BIN_HEIGHT:
        return None
    sisters = C.sister_clades(node)
    if len(sisters) < 2:
        return None

    candidates = []
    for species in sisters[0].leaves():
        path = C.median_path(species, a_species, dist)
        if path is None or species not in bins:
            continue
        candidates.append({"species": species, "path": path,
                           "expected_ks": C.path_to_ks(path),
                           "bin_height": bins[species],
                           "lineage": lineages.get(species, "")})
    rungs = ladder(candidates)
    if len(rungs) < MIN_RUNGS:
        return None

    # The rate adjustment subtracts Ks(sister, outgroup), so it is the SISTER-outgroup
    # distances - not the focal-outgroup one - that decide how saturated the subtracted term
    # is. Take the outgroup that keeps the worst of them shortest, among those whose own
    # paranome is good enough for the species to be trusted at all.
    chosen, rung_species = None, [c["species"] for _, c in rungs]
    for species in sisters[1].leaves():
        if bins.get(species, 0) < MIN_OUTGROUP_BIN_HEIGHT:
            continue
        to_sisters = [dist(s, species) for s in rung_species]
        if any(d is None for d in to_sisters):
            continue
        worst = C.path_to_ks(max(to_sisters))
        if chosen is None or worst < chosen["worst_sister_ks"]:
            out_path = C.median_path(species, a_species, dist)
            chosen = {"species": species, "path": out_path,
                      "expected_ks": None if out_path is None else C.path_to_ks(out_path),
                      "bin_height": bins[species], "lineage": lineages.get(species, ""),
                      "worst_sister_ks": worst}
    if chosen is None:
        return None

    return {"index": index, "focal": focal, "focal_bin_height": bins[focal],
            "a_species": a_species, "rungs": rungs, "outgroup": chosen}


def write_config(design, meta, out_dir, db_address, analysis_tree, dist):
    focal, outgroup = design["focal"], design["outgroup"]["species"]
    sisters = [c["species"] for _, c in design["rungs"]]
    members = [outgroup, focal, *sisters]

    names = {s: short(s) for s in members}
    if len(set(names.values())) != len(names):
        raise SystemExit(f"branch {design['index']}: informal names collide: {sorted(names.values())}")

    # The real induced topology, not a flattened sister clade: ksrates requires every internal
    # node to have exactly two children and rejects a polytomy outright. The nesting inside the
    # sister clade costs a few extra ortholog pairs (ksrates computes them for its distance
    # tree) but does not affect the comparison, since all the sisters still share one
    # divergence node with the focal species.
    topology = resolve_polytomies(prune(analysis_tree, set(members)), dist)
    newick = render(topology, names) + ";"
    if set(subtree_leaves(topology)) != set(members):
        raise SystemExit(f"branch {design['index']}: pruned tree lost species")
    focal_side = next((child for child in topology if outgroup not in subtree_leaves(child)), None)
    if focal_side is None or set(subtree_leaves(focal_side)) != {focal, *sisters}:
        raise SystemExit(f"branch {design['index']}: the outgroup is not outside the focal+sister clade")
    pad = max(len(n) for n in names.values())
    latin = ",\n              ".join(
        f"{names[s]:<{pad}} : {meta[s]['latin_name']}" for s in members)
    fasta = ",\n                  ".join(
        f"{names[s]:<{pad}} : {meta[s]['fasta_filename']}" for s in members)

    path = os.path.join(out_dir, "configs", f"config_sisterrate_b{design['index']}_{names[focal]}.txt")
    with open(path, "w") as fh:
        fh.write(TEMPLATE.format(focal=names[focal], newick=newick, latin=latin,
                                 fasta=fasta, db_address=db_address))
    return path, names


ROW_COLUMNS = ["branch_index", "focal_species", "focal_bin_height", "role", "rung",
               "species", "informal_name", "path", "expected_ks", "reliability",
               "expected_ks_to_outgroup", "bin_height", "lineage"]


def design_rows(design, names, dist):
    focal, fbh = design["focal"], design["focal_bin_height"]
    out = design["outgroup"]
    for rung, cand in design["rungs"]:
        to_out = dist(cand["species"], out["species"])
        yield [design["index"], focal, f"{fbh:.1f}", "sister", f"{rung:.2f}",
               cand["species"], names[cand["species"]], f"{cand['path']:.4f}",
               f"{cand['expected_ks']:.3f}", C.path_reliability(cand["path"]),
               "" if to_out is None else f"{C.path_to_ks(to_out):.3f}",
               f"{cand['bin_height']:.1f}", cand["lineage"]]
    yield [design["index"], focal, f"{fbh:.1f}", "outgroup", "",
           out["species"], names[out["species"]],
           "" if out["path"] is None else f"{out['path']:.4f}",
           "" if out["expected_ks"] is None else f"{out['expected_ks']:.3f}",
           C.path_reliability(out["path"]) or "", "", f"{out['bin_height']:.1f}", out["lineage"]]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default=os.path.join(C.ROOT, "ks_analysis/1kp/sister_rate_test"))
    ap.add_argument("--branches", help="comma-separated traversal indices; default: the chosen five")
    ap.add_argument("--min-outgroup-bin-height", type=float, default=MIN_OUTGROUP_BIN_HEIGHT,
                    help="an outgroup with a poorer paranome than this is not trusted as a species")
    ap.add_argument("--survey", action="store_true",
                    help="list every branch that admits a ladder and exit, without writing anything")
    ap.add_argument("--species-metadata",
                    default=os.path.join(C.ROOT, "ks_analysis/1kp/species_metadata.tsv"))
    ap.add_argument("--db-address",
                    default=os.path.join(C.ROOT, "ks_analysis/1kp/ks_database/ks_database_server_address.txt"))
    args = ap.parse_args()

    _, dist = C.load_concat_tree()
    bins = C.median_bin_heights()
    analysis_tree = C.load_analysis_tree()
    nodes = C.branch_nodes(analysis_tree)
    lineages = C.groups()
    globals()["MIN_OUTGROUP_BIN_HEIGHT"] = args.min_outgroup_bin_height

    if args.survey:
        print(f"{'branch':>7}{'focal':>16}{'binH':>7}{'rungs':>7}{'outgroup':>10}{'outBinH':>9}"
              f"{'worstSisOutKs':>15}  ladder (species expKs/binHeight)")
        for index, node in enumerate(nodes):
            design = branch_design(index, node, dist, bins, lineages)
            if design is None:
                continue
            rung_text = "  ".join(f"{c['species']}:{c['expected_ks']:.2f}/{c['bin_height']:.0f}"
                                  for _, c in design["rungs"])
            out = design["outgroup"]
            print(f"{index:>7}{design['focal']:>16}{design['focal_bin_height']:>7.0f}"
                  f"{len(design['rungs']):>7}{out['species']:>10}{out['bin_height']:>9.0f}"
                  f"{out['worst_sister_ks']:>15.2f}  {rung_text}")
        return

    wanted = (tuple(int(x) for x in args.branches.split(",")) if args.branches
              else DEFAULT_BRANCHES)
    meta = read_metadata(args.species_metadata)
    os.makedirs(os.path.join(args.out_dir, "configs"), exist_ok=True)

    rows, configs = [], []
    for index in wanted:
        design = branch_design(index, nodes[index], dist, bins, lineages)
        if design is None:
            raise SystemExit(f"branch {index} does not admit a ladder")
        path, names = write_config(design, meta, args.out_dir, args.db_address,
                                   analysis_tree, dist)
        configs.append(path)
        rows.extend(design_rows(design, names, dist))

    table = os.path.join(args.out_dir, "sister_rate_test_design.tsv")
    with open(table, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(ROW_COLUMNS)
        writer.writerows(rows)

    # An explicit manifest, not a glob: configs from an earlier choice of branches stay on disk
    # and would otherwise shift the SLURM array indices onto the wrong experiment.
    manifest = os.path.join(args.out_dir, "configs.txt")
    with open(manifest, "w") as fh:
        fh.write("\n".join(configs) + "\n")

    # Per branch ksrates needs, for k sisters and one outgroup: k focal-sister pairs,
    # k sister-outgroup pairs and 1 focal-outgroup pair - plus roughly one more per internal
    # node of the sister clade, which it computes for its distance tree.
    sisters = sum(1 for r in rows if r[3] == "sister")
    trio_pairs = 2 * sisters + len(configs)
    print(f"{len(configs)} configs, {sisters} sister rungs, "
          f"{trio_pairs} trio ortholog pairs (expect roughly "
          f"{trio_pairs + sisters - 2 * len(configs)} in total)")
    print(f"design table: {table}")
    print(f"manifest:     {manifest}   (launch with --array=0-{len(configs) - 1})")
    for path in configs:
        print(f"  {path}")


if __name__ == "__main__":
    main()
