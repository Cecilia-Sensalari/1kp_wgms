#!/usr/bin/env python3
"""Split the branch and dataset tables into one self-contained pair per lineage.

The 1KP analysis is large enough to be worth running on several clusters at once, one lineage
each. This writes, for every lineage, a directory holding the two tables the wrapper needs:

    <out-dir>/<lineage>/focal_branch_to_dataset.tsv
    <out-dir>/<lineage>/subtree_datasets.tsv

and a run of that shard is then just

    dataset_list = "<out-dir>/<lineage>/subtree_datasets.tsv"

in the analysis nextflow.config. The second file has to sit beside the first under exactly that
name, because main.nf derives the branch table as `<parent of dataset_list>/focal_branch_to_
dataset.tsv` rather than taking it as a parameter.

A branch belongs to the lineage of its first A species. A dataset goes to every lineage that has
a branch using it; measured on the 2026-11-06 tables that duplicated one dataset in 1755, because
a dataset is built around its own clade and rarely serves two lineages.

SHARDS DO NOT NEED A SHARED Ks DATABASE. Each lineage can run against its own server and the
results be merged afterwards: over the six lineages the shards ask for 9291 ortholog pairs where
8 distinct ones would do (9101), and 1228 paralog species where 1178 would do - 2% and 4% of
duplicated work, because the species a dataset reaches for are nearly always inside its own
lineage. Sharing one database across clusters removes even that, if the clusters can reach it.

    python3 split_tables_by_lineage.py --branches <focal_branch_to_dataset.tsv> --out-dir <dir>

The dataset table is read from beside the branch table unless --datasets is given.
"""

from __future__ import annotations

import argparse
import collections
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "calibration"))
import calib_common as C


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--branches", required=True,
                    help="focal_branch_to_dataset.tsv of the run to split")
    ap.add_argument("--datasets", default=None,
                    help="subtree_datasets.tsv (default: beside --branches)")
    ap.add_argument("--out-dir", required=True, help="directory to write the shards into")
    ap.add_argument("--lineages", nargs="*", default=None,
                    help="only these lineages (default: all that occur)")
    args = ap.parse_args()

    datasets_path = args.datasets or os.path.join(os.path.dirname(os.path.abspath(args.branches)),
                                                  "subtree_datasets.tsv")
    meta = C.species_metadata()
    groups = C.groups(meta)
    resolved = {}

    def lineage_of(name):
        if name not in resolved:
            tip = C.resolve_tip(name.split()[-1], meta)
            resolved[name] = groups.get(tip) if tip else None
        return resolved[name]

    with open(args.branches) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        branch_columns, branch_rows = reader.fieldnames, list(reader)
    with open(datasets_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        dataset_columns, dataset_rows = reader.fieldnames, list(reader)
    datasets_by_id = {row["dataset_id"]: row for row in dataset_rows}

    by_lineage = collections.defaultdict(list)
    unassigned = 0
    for row in branch_rows:
        first_a = next((s for s in row["A_species"].split(",") if s), None)
        lineage = lineage_of(first_a) if first_a else None
        if lineage is None:
            unassigned += 1
            continue
        by_lineage[lineage].append(row)

    wanted = args.lineages or sorted(by_lineage)
    print(f"{'lineage':<14}{'branches':>10}{'datasets':>10}  directory")
    seen_datasets = collections.Counter()
    for lineage in wanted:
        rows = by_lineage.get(lineage)
        if not rows:
            print(f"{lineage:<14}{'-':>10}{'-':>10}  (no branch)")
            continue
        ids = []
        for row in rows:
            if row["dataset_id"] not in ids:
                ids.append(row["dataset_id"])
        seen_datasets.update(ids)

        folder = os.path.join(args.out_dir, lineage)
        os.makedirs(folder, exist_ok=True)
        with open(os.path.join(folder, "focal_branch_to_dataset.tsv"), "w", newline="") as fh:
            writer = csv.DictWriter(fh, delimiter="\t", fieldnames=branch_columns)
            writer.writeheader()
            writer.writerows(rows)
        with open(os.path.join(folder, "subtree_datasets.tsv"), "w", newline="") as fh:
            writer = csv.DictWriter(fh, delimiter="\t", fieldnames=dataset_columns)
            writer.writeheader()
            writer.writerows(datasets_by_id[i] for i in ids if i in datasets_by_id)
        print(f"{lineage:<14}{len(rows):>10}{len(ids):>10}  {folder}")

    shared = [i for i, n in seen_datasets.items() if n > 1]
    print(f"\n{sum(seen_datasets.values())} dataset slots over the shards, "
          f"{len(seen_datasets)} distinct, {len(shared)} used by more than one lineage")
    if unassigned:
        print(f"{unassigned} branches had no lineage for their first A species and were skipped")


if __name__ == "__main__":
    main()
