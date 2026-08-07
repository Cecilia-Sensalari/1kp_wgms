#!/usr/bin/env python3
"""Annotate the rooted 1KP species tree with WGM/WGD IDs using ETE3.

[Generated via AI, tested by Cecilia]

The 1KP/Barker release provides the WGM/WGD summary table and the figure PDF,
but not a machine-readable Newick tree with WGM labels already attached to
branches. This script reconstructs a practical branch annotation from the
supplementary tables:

1. Supplementary Table 2 defines the official WGM/WGD IDs.
2. Supplementary Table 3 lists which 1KP taxa carry each WGM/WGD in their WGD
   history columns.
3. The rooted ASTRAL species tree provides the topology.
4. Each WGM/WGD is placed on the MRCA branch of its supporting taxa.

Outputs:

- an NHX-annotated Newick tree containing ``WGM1=...``, ``WGM2=...`` branch features
- a TSV audit table explaining every placement and its support taxa

ETE3 is used for all tree parsing, MRCA calculation, and tree writing. The
supplementary source tables are read from tab-separated text files.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

try:
    from ete3 import Tree
except ImportError as exc:
    raise SystemExit(
        "Missing dependency: ete3. Install it in this environment with `pip install ete3`."
    ) from exc


# Project defaults. Override these on the command line if running from another
# checkout, another platform, or a test directory.
ROOT = Path(r"/group/esb/cesen/1kp")
DEFAULT_TABLE2 = ROOT / "source_data/1.species_dataset/1kp_paper_2019_suptab2_wgm_EDITED.tsv"
DEFAULT_TABLE3 = ROOT / "source_data/1.species_dataset/1kp_paper_2019_suptab3_ks_EDITED.tsv"
DEFAULT_TREE = (
    ROOT
    / "source_data/3.phylogenetic_tree/1kp_trees/"
    "astral_trees_33_percent-FAA_estimated_species_tree.rooted.tree"
)


# ---------------------------------------------------------------------------
# TSV reading helper
# ---------------------------------------------------------------------------


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    """Read a UTF-8 TSV table as a list of row dictionaries.

    ``utf-8-sig`` also accepts files with a byte-order mark. Empty rows are
    omitted, matching the behavior expected by the table interpretation code.
    """
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            {key: value or "" for key, value in row.items() if key is not None}
            for row in reader
            if any(value for value in row.values())
        ]


# ---------------------------------------------------------------------------
# 1KP supplementary table interpretation
# ---------------------------------------------------------------------------


def split_wgm_ids(value: str) -> list[str]:
    """Split a WGM/WGD cell into one or more event IDs.

    Most Table 3 cells contain a single event ID, but the parser also accepts
    comma-, semicolon-, or whitespace-separated lists. ``NA`` and blank cells are
    ignored.
    """
    if not value:
        return []

    ids = []
    for item in re.split(r"[,;]\s*|\s+", value.strip()):
        item = item.strip()
        if item and item.upper() != "NA":
            ids.append(item)
    return ids


def wgm_ids_from_table2(path: Path) -> set[str]:
    """Read the official WGM/WGD IDs from Supplementary Table 2.

    Table 2 is used as a whitelist. By default, IDs seen in Table 3 but absent
    from Table 2 are skipped, because Table 2 is the paper's summary list of
    inferred WGMs/WGDs.
    """
    return {row["WGD ID"] for row in read_tsv_rows(path) if row.get("WGD ID")}


def wgm_taxa_from_table3(path: Path, skipped_list_not_valid_wgms: list[str]) -> dict[str, set[str]]:
    """Build a mapping from WGM/WGD ID to supporting 1KP species codes.

    Supplementary Table 3 is organized by taxon. The columns ``WGD 1``, ``WGD 2``,
    and ``WGD 3`` describe the WGM/WGD history inferred for that taxon. Inverting
    those columns gives the taxon set used to place each event on the tree.
    NOTE: Ignores ID "B2(possible WGD)" because it is not a Ks-derived WGM and is not listed in Table 2.
    """
    out: dict[str, set[str]] = defaultdict(set)
    for row in read_tsv_rows(path):
        code = row.get("1KP Code", "").strip()
        if not code:
            continue
        for column in ("WGD 1", "WGD 2", "WGD 3"):
            for wgm_id in split_wgm_ids(row.get(column, "")):
                # Catch edge cases that have to be manually curated
                # 1. In Table 3, the "WGD2" column has cells with "B2(possible WGD)", which is actually 
                # a MAPS-derived ID and not a Ks plots-derived ID. Therefore Table 2 (listing Ks-derived WGMs)
                # does not list that WGM ID. Since we are only dealing with Ks-derived WGMs, we will skip that WGM ID.
                if wgm_id == "B2(possible" or wgm_id == "WGD)":
                    if wgm_id not in skipped_list_not_valid_wgms:
                        skipped_list_not_valid_wgms.append(wgm_id)
                    continue
                out[wgm_id].add(code)
    return out


# ---------------------------------------------------------------------------
# Tree placement helpers
# ---------------------------------------------------------------------------


def nhx_escape(value: str) -> str:
    """Escape separators that would confuse NHX feature syntax.

    WGM IDs are usually compact strings such as ``AMBOalpha``. This helper makes
    the output more robust if a future table contains spaces or punctuation.
    """
    return (
        value.replace("\\", "_")
        .replace(":", "_")
        .replace(",", "_")
        .replace(";", "_")
        .replace("|", "_")
        .replace("]", "_")
        .replace("[", "_")
        .replace(" ", "_")
    )


def assign_node_ids(tree: Tree) -> None:
    """Attach stable numeric IDs to ETE nodes for the audit table.

    ETE nodes do not come with persistent IDs. The IDs written here are not
    biological identifiers; they are a convenient way to connect a TSV row to a
    branch in this specific run.
    """
    for node_id, node in enumerate(tree.traverse("postorder")):
        node.add_feature("placement_node_id", node_id)


def find_mrca(tree: Tree, taxa: Iterable[str], tree_tip_names: set[str]):
    """Find the MRCA branch for a set of support taxa.

    Returns ``(node, present, missing)``. ``present`` contains support taxa found
    in the tree; ``missing`` contains support taxa listed in Table 3 but absent
    from this tree. If only one support taxon is present, the event is placed on
    that terminal branch.
    """
    taxa = set(taxa)
    present = sorted(taxa.intersection(tree_tip_names))
    missing = sorted(taxa.difference(tree_tip_names))
    if not present:
        return None, present, missing
    if len(present) == 1:
        return tree & present[0], present, missing
    return tree.get_common_ancestor(present), present, missing


def add_wgm_annotation(node, wgm_id: str) -> None:
    """Append one WGM/WGD ID to an ETE node as NHX-compatible features.

    NHX stores annotations as colon-separated ``key=value`` fields. Some tools,
    including iTOL, can misread a custom value containing ``|``. To keep
    multi-WGM branches parser-friendly, each event gets its own key:
    ``WGM1=A:WGM2=B``. ``WGM_label`` is a comma-separated convenience label for
    scripts that want one display string.
    """
    current = []
    if hasattr(node, "WGM_label") and node.WGM_label:
        current = [item.strip() for item in str(node.WGM_label).split(",") if item.strip()]
    elif hasattr(node, "WGM") and node.WGM:
        current = [item.strip() for item in re.split(r"[|,]", str(node.WGM)) if item.strip()]

    current.append(nhx_escape(wgm_id))
    current = sorted(set(current))

    for index, event_id in enumerate(current, start=1):
        node.add_feature(f"WGM{index}", event_id)
    node.add_feature("WGM_count", len(current))
    node.add_feature("WGM_label", ",".join(current))

def placement_note(present_count: int, mrca_tip_count: int) -> tuple[bool, str]:
    """Describe how precise the MRCA placement is.

    Exact clade matches are the easiest to trust. Single-taxon and subset-of-MRCA
    placements are still useful annotations, but should be reviewed before being
    treated as final curated placements.
    """
    exact = mrca_tip_count == present_count
    if present_count == 1:
        return exact, "single sampled support taxon; placed on terminal branch"
    if exact:
        return exact, "support taxa exactly match MRCA clade"
    return exact, "support taxa are a subset of MRCA clade"


# ---------------------------------------------------------------------------
# Outputs and command-line interface
# ---------------------------------------------------------------------------


def default_output_paths(tree_path: Path) -> tuple[Path, Path]:
    """Derive default output filenames beside the input tree."""
    stem = tree_path.name
    if stem.endswith(".tree"):
        stem = stem[: -len(".tree")]
    out_tree = tree_path.with_name(f"{stem}.wgm_suptab3_mrca.nhx.tree")
    out_tsv = tree_path.with_name(f"{stem}.wgm_suptab3_mrca.tsv")
    return out_tree, out_tsv


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    """Write the audit table that explains every tree annotation."""
    fieldnames = [
        "WGM",
        "node_id",
        "support_taxa_count",
        "present_taxa_count",
        "missing_taxa_count",
        "mrca_tip_count",
        "exact_support_clade",
        "placement_note",
        "support_taxa",
        "missing_taxa",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    """Create the command-line interface and document all override points."""
    parser = argparse.ArgumentParser(
        description=(
            "Place 1KP WGM/WGD IDs onto a rooted Newick tree using the MRCA "
            "of supporting taxa listed in Supplementary Table 3. Tree parsing "
            "and MRCA calculations are handled by ETE3."
        )
    )
    parser.add_argument("--tree", type=Path, default=DEFAULT_TREE, help="Rooted 1KP species tree in Newick format.")
    parser.add_argument("--wgm-table", type=Path, default=DEFAULT_TABLE2, help="Supplementary Table 2 WGM/WGD list TSV.")
    parser.add_argument("--ks-table", type=Path, default=DEFAULT_TABLE3, help="Supplementary Table 3 Ks/WGD history TSV.")
    parser.add_argument("--out-tree", type=Path, default=None, help="Output NHX-annotated Newick tree.")
    parser.add_argument("--out-tsv", type=Path, default=None, help="Output placement audit TSV.")
    parser.add_argument(
        "--ete-format",
        type=int,
        default=0,
        help=(
            "ETE3 Newick format code for input/output. Default 0 is intended "
            "for the ASTRAL-style tree with internal support values and branch lengths."
        ),
    )
    parser.add_argument(
        "--include-unlisted",
        action="store_true",
        help="Also annotate event IDs found in Table 3 but absent from the Table 2 WGM list.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Coordinate table parsing, MRCA placement, and output writing."""
    args = build_arg_parser().parse_args(argv)
    out_tree, out_tsv = default_output_paths(args.tree)
    if args.out_tree is not None:
        out_tree = args.out_tree
    if args.out_tsv is not None:
        out_tsv = args.out_tsv

    # ETE3 handles the real tree work: parsing Newick, finding leaves, and
    # calculating MRCAs. Node IDs are only for human review in the TSV.
    tree = Tree(str(args.tree), format=args.ete_format)
    assign_node_ids(tree)
    tree_tip_names = set(tree.get_leaf_names())

    # Inizialize the audit table and bookkeeping for skipped events when listing WGM IDs
    # that are not in Table 2 or have no MRCA in the tree.
    audit_rows: list[dict[str, object]] = []
    skipped = 0
    skipped_list_not_valid_wgms = []
    skipped_list_no_mrca = []
    annotated_node_ids: set[int] = set()

    # Read the two evidence tables: Table 2 defines the valid event IDs, while
    # Table 3 tells us which species codes support each event.
    valid_wgms = wgm_ids_from_table2(args.wgm_table)
    wgm_to_taxa = wgm_taxa_from_table3(args.ks_table, skipped_list_not_valid_wgms)

    # Place every WGM independently. If several WGMs map to one branch, each event
    # gets its own NHX key: WGM1=..., WGM2=..., etc.
    for wgm_id in sorted(wgm_to_taxa):
        if not args.include_unlisted and wgm_id not in valid_wgms:
            skipped += 1
            skipped_list_not_valid_wgms.append(wgm_id)
            continue

        mrca, present, missing = find_mrca(tree, wgm_to_taxa[wgm_id], tree_tip_names)
        if mrca is None:
            skipped += 1
            skipped_list_no_mrca.append(wgm_id)
            continue

        add_wgm_annotation(mrca, wgm_id)
        annotated_node_ids.add(mrca.placement_node_id)

        mrca_tip_count = len(mrca.get_leaf_names())
        exact, note = placement_note(len(present), mrca_tip_count)
        audit_rows.append(
            {
                "WGM": wgm_id,
                "node_id": mrca.placement_node_id,
                "support_taxa_count": len(wgm_to_taxa[wgm_id]),
                "present_taxa_count": len(present),
                "missing_taxa_count": len(missing),
                "mrca_tip_count": mrca_tip_count,
                "exact_support_clade": str(exact).lower(),
                "placement_note": note,
                "support_taxa": ",".join(present),
                "missing_taxa": ",".join(missing),
            }
        )

    # Tell ETE3 which custom node features should be written as NHX tags. We list
    # many WGM slots so branches with multiple events retain separate keys.
    out_tree.parent.mkdir(parents=True, exist_ok=True)
    max_wgms_per_node = max((int(node.WGM_count) for node in tree.traverse() if hasattr(node, "WGM_count")), default=0)
    nhx_features = ["WGM_label", "WGM_count"] + [f"WGM{i}" for i in range(1, max_wgms_per_node + 1)]
    tree.write(outfile=str(out_tree), format=args.ete_format, features=nhx_features)
    write_tsv(out_tsv, sorted(audit_rows, key=lambda row: str(row["WGM"])))

    print(f"Wrote annotated tree: {out_tree}")
    print(f"Wrote audit table:    {out_tsv}")
    print(f"Table 2 WGM IDs:      {len(valid_wgms)}")
    print(f"Table 3 WGM IDs:      {len(wgm_to_taxa)}")

    missing_valid_wgms = set(valid_wgms) - set(wgm_to_taxa)
    if missing_valid_wgms != set():
        print(f"Missing valid WGMs (Tab2 - Tab3): {missing_valid_wgms}")

    print(f"Annotated WGM IDs:    {len(audit_rows)}")
    print(f"Annotated nodes:      {len(annotated_node_ids)}")

    if skipped:
        print(f"Skipped WGM IDs:  {skipped}", file=sys.stderr)
        print(f"  Not in Table 2 (WGM list):     {skipped_list_not_valid_wgms}", file=sys.stderr)
        print(f"  No MRCAs found           :     {skipped_list_no_mrca}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
