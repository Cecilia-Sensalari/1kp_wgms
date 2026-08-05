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

ETE3 is used for all tree parsing, MRCA calculation, and tree writing. XLSX
files are read directly with Python's standard library so that only the tree
library needs to be installed.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable
from xml.etree import ElementTree as ET

try:
    from ete3 import Tree
except ImportError as exc:
    raise SystemExit(
        "Missing dependency: ete3. Install it in this environment with `pip install ete3`."
    ) from exc


# Project defaults. Override these on the command line if running from another
# checkout, another platform, or a test directory.
ROOT = Path(r"/group/esb/cesen/1kp")
DEFAULT_TABLE2 = ROOT / "source_data/1.species_dataset/1kp_paper_2019_suptab2_wgm.xlsx"
DEFAULT_TABLE3 = ROOT / "source_data/1.species_dataset/1kp_paper_2019_suptab3_ks.xlsx"
DEFAULT_TREE = (
    ROOT
    / "source_data/3.phylogenetic_tree/1kp_trees/"
    "astral_trees_33_percent-FAA_estimated_species_tree.rooted.tree"
)


# Namespace abbreviations used inside the XLSX XML files.
XML_NS = {
    "a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}


# ---------------------------------------------------------------------------
# XLSX reading helpers
# ---------------------------------------------------------------------------
#
# XLSX files are ZIP archives containing XML. The 1KP supplementary workbooks
# are simple enough that we can read the first worksheet directly. This avoids
# requiring pandas/openpyxl on top of ETE3.


def column_index(cell_ref: str) -> int:
    """Convert an Excel cell reference such as ``C42`` to a 1-based column.

    Spreadsheet rows are sparse in the XML: blank cells may simply be omitted.
    Converting ``A``, ``B``, ..., ``AA`` to numeric column indices lets us put
    each value under the correct header even when intermediate cells are blank.
    """
    match = re.match(r"[A-Z]+", cell_ref)
    if match is None:
        raise ValueError(f"Cannot parse Excel cell reference: {cell_ref}")

    index = 0
    for char in match.group(0):
        index = index * 26 + (ord(char) - ord("A") + 1)
    return index


def cell_text(cell: ET.Element, shared_strings: list[str]) -> str:
    """Return the displayed text for one XLSX worksheet cell.

    Excel stores repeated strings in a shared string table and stores numbers
    directly in the worksheet XML. This helper hides that difference so the rest
    of the script can treat all cells as strings.
    """
    cell_type = cell.attrib.get("t")
    if cell_type == "s":
        value = cell.findtext("a:v", default="", namespaces=XML_NS)
        return shared_strings[int(value)] if value != "" else ""
    if cell_type == "inlineStr":
        return "".join(cell.itertext()).strip()
    return cell.findtext("a:v", default="", namespaces=XML_NS).strip()


def read_xlsx_rows(path: Path, sheet_xml: str = "xl/worksheets/sheet1.xml") -> list[dict[str, str]]:
    """Read a simple XLSX worksheet as a list of row dictionaries.

    The first row is interpreted as the header. Every following non-empty row is
    returned as ``{header: value}``, with blank cells represented by empty
    strings. This is enough for the 1KP supplementary tables used here.
    """
    with zipfile.ZipFile(path) as zf:
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            shared_root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for item in shared_root.findall("a:si", XML_NS):
                shared_strings.append("".join(item.itertext()))

        sheet = ET.fromstring(zf.read(sheet_xml))
        raw_rows: list[dict[int, str]] = []
        for row in sheet.findall(".//a:sheetData/a:row", XML_NS):
            values: dict[int, str] = {}
            for cell in row.findall("a:c", XML_NS):
                ref = cell.attrib.get("r", "")
                if ref:
                    values[column_index(ref)] = cell_text(cell, shared_strings)
            raw_rows.append(values)

    if not raw_rows:
        return []

    headers = {index: value for index, value in raw_rows[0].items() if value}
    rows: list[dict[str, str]] = []
    for raw in raw_rows[1:]:
        row = {header: raw.get(index, "") for index, header in headers.items()}
        if any(value != "" for value in row.values()):
            rows.append(row)
    return rows


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
    return {row["WGD ID"] for row in read_xlsx_rows(path) if row.get("WGD ID")}


def wgm_taxa_from_table3(path: Path) -> dict[str, set[str]]:
    """Build a mapping from WGM/WGD ID to supporting 1KP species codes.

    Supplementary Table 3 is organized by taxon. The columns ``WGD 1``, ``WGD 2``,
    and ``WGD 3`` describe the WGM/WGD history inferred for that taxon. Inverting
    those columns gives the taxon set used to place each event on the tree.
    """
    out: dict[str, set[str]] = defaultdict(set)
    for row in read_xlsx_rows(path):
        code = row.get("1KP Code", "").strip()
        if not code:
            continue
        for column in ("WGD 1", "WGD 2", "WGD 3"):
            for wgm_id in split_wgm_ids(row.get(column, "")):
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
    parser.add_argument("--wgm-table", type=Path, default=DEFAULT_TABLE2, help="Supplementary Table 2 WGM/WGD list XLSX.")
    parser.add_argument("--ks-table", type=Path, default=DEFAULT_TABLE3, help="Supplementary Table 3 Ks/WGD history XLSX.")
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

    # Read the two evidence tables: Table 2 defines the valid event IDs, while
    # Table 3 tells us which species codes support each event.
    valid_wgms = wgm_ids_from_table2(args.wgm_table)
    wgm_to_taxa = wgm_taxa_from_table3(args.ks_table)

    # ETE3 handles the real tree work: parsing Newick, finding leaves, and
    # calculating MRCAs. Node IDs are only for human review in the TSV.
    tree = Tree(str(args.tree), format=args.ete_format)
    assign_node_ids(tree)
    tree_tip_names = set(tree.get_leaf_names())

    audit_rows: list[dict[str, object]] = []
    skipped = 0
    annotated_node_ids: set[int] = set()

    # Place every WGM independently. If several WGMs map to one branch, each event
    # gets its own NHX key: WGM1=..., WGM2=..., etc.
    for wgm_id in sorted(wgm_to_taxa):
        if not args.include_unlisted and wgm_id not in valid_wgms:
            skipped += 1
            continue

        mrca, present, missing = find_mrca(tree, wgm_to_taxa[wgm_id], tree_tip_names)
        if mrca is None:
            skipped += 1
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
    print(f"Annotated WGM IDs:    {len(audit_rows)}")
    print(f"Annotated nodes:      {len(annotated_node_ids)}")
    if skipped:
        print(f"Skipped WGM IDs:      {skipped}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
