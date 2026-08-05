#!/usr/bin/env python3
"""Annotate a 1KP species tree with WGM/WGD IDs from supplementary tables.

[Generated via AI, tested by Cecilia]

The script places each event on the MRCA branch of the 1KP species codes that
support that event in Supplementary Table 3. It writes:

1. an NHX-commented Newick tree with ``WGM=...`` on annotated branches
2. a TSV audit file describing every placement

It uses only the Python standard library, including direct XLSX parsing through
zip/xml, so no openpyxl/pandas installation is required.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import zipfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable
from xml.etree import ElementTree as ET


#ROOT = Path(r"\\psb.ugent.be\shares\biocomp\groups\group_esb\cesen\1kp")
ROOT = Path(r"/group/esb/cesen/1kp")
DEFAULT_TABLE2 = ROOT / "source_data/1.species_dataset/1kp_paper_2019_suptab2_wgm.xlsx"
DEFAULT_TABLE3 = ROOT / "source_data/1.species_dataset/1kp_paper_2019_suptab3_ks.xlsx"
DEFAULT_TREE = (
    ROOT
    / "source_data/3.phylogenetic_tree/1kp_trees/"
    "astral_trees_33_percent-FAA_estimated_species_tree.rooted.tree"
)


XML_NS = {
    "a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}


# ---------------------------------------------------------------------------
# Minimal tree representation and Newick parser/writer
# ---------------------------------------------------------------------------
#
# The script avoids BioPython/ETE dependencies, so it keeps just enough tree
# structure to parse Newick, walk parent/child relationships, calculate MRCAs,
# and write the tree back with NHX comments.

@dataclass(eq=False)
class Node:
    """Small Newick tree node used for MRCA calculations."""

    name: str = ""
    length: str = ""
    comments: list[str] = field(default_factory=list)
    children: list["Node"] = field(default_factory=list)
    parent: "Node | None" = None
    node_id: int = -1
    tip_names: set[str] = field(default_factory=set)

    @property
    def is_tip(self) -> bool:
        return not self.children


class NewickParser:
    """Recursive-descent parser for the simple Newick produced by ASTRAL."""

    def __init__(self, text: str) -> None:
        self.text = text.strip()
        self.i = 0
        self.next_node_id = 0

    def parse(self) -> Node:
        """Parse the full Newick string and return the root node."""
        root = self._subtree()
        self._skip_ws()
        if self.i < len(self.text) and self.text[self.i] == ";":
            self.i += 1
        self._skip_ws()
        if self.i != len(self.text):
            raise ValueError(f"Unexpected trailing Newick text at offset {self.i}")
        return root

    def _new_node(self, name: str = "") -> Node:
        node = Node(name=name, node_id=self.next_node_id)
        self.next_node_id += 1
        return node

    def _subtree(self) -> Node:
        """Parse either an internal subtree ``(...)label:length`` or a tip."""
        self._skip_ws()
        if self._peek() == "(":
            self.i += 1
            node = self._new_node()
            while True:
                child = self._subtree()
                child.parent = node
                node.children.append(child)
                self._skip_ws()
                char = self._peek()
                if char == ",":
                    self.i += 1
                    continue
                if char == ")":
                    self.i += 1
                    break
                raise ValueError(f"Expected ',' or ')' at offset {self.i}")
            node.name = self._read_label()
            node.comments.extend(self._read_comments())
            node.length = self._read_length()
            return node

        node = self._new_node(self._read_label())
        if not node.name:
            raise ValueError(f"Expected tip label at offset {self.i}")
        node.comments.extend(self._read_comments())
        node.length = self._read_length()
        return node

    def _peek(self) -> str:
        return self.text[self.i] if self.i < len(self.text) else ""

    def _skip_ws(self) -> None:
        while self.i < len(self.text) and self.text[self.i].isspace():
            self.i += 1

    def _read_label(self) -> str:
        self._skip_ws()
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in ",():;[]":
            self.i += 1
        return self.text[start : self.i].strip()

    def _read_comments(self) -> list[str]:
        """Preserve any existing Newick comments instead of discarding them."""
        comments: list[str] = []
        self._skip_ws()
        while self._peek() == "[":
            start = self.i
            depth = 0
            while self.i < len(self.text):
                if self.text[self.i] == "[":
                    depth += 1
                elif self.text[self.i] == "]":
                    depth -= 1
                    if depth == 0:
                        self.i += 1
                        comments.append(self.text[start : self.i])
                        break
                self.i += 1
            else:
                raise ValueError(f"Unclosed Newick comment starting at offset {start}")
            self._skip_ws()
        return comments

    def _read_length(self) -> str:
        self._skip_ws()
        if self._peek() != ":":
            return ""
        self.i += 1
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in ",();":
            self.i += 1
        return self.text[start : self.i].strip()


def quote_newick_label(label: str) -> str:
    """Quote labels only when required by Newick syntax."""
    if not label:
        return ""
    if re.search(r"[\s,:;()\[\]']", label):
        return "'" + label.replace("'", "''") + "'"
    return label


def nhx_escape(value: str) -> str:
    """Escape separators that would confuse NHX consumers."""
    return (
        value.replace("\\", "_")
        .replace(":", "_")
        .replace(",", "_")
        .replace(";", "_")
        .replace("]", "_")
        .replace("[", "_")
        .replace(" ", "_")
    )


def to_newick(node: Node, annotations: dict[int, list[str]]) -> str:
    """Serialize the tree, adding ``[&&NHX:WGM=...]`` to annotated branches."""
    comments = list(node.comments)
    if node.node_id in annotations:
        wgms = "|".join(nhx_escape(wgm) for wgm in sorted(annotations[node.node_id]))
        comments.append(f"[&&NHX:WGM={wgms}]")

    comment_text = "".join(comments)
    length_text = f":{node.length}" if node.length else ""
    if node.is_tip:
        return f"{quote_newick_label(node.name)}{comment_text}{length_text}"

    child_text = ",".join(to_newick(child, annotations) for child in node.children)
    return f"({child_text}){quote_newick_label(node.name)}{comment_text}{length_text}"


def iter_nodes_postorder(node: Node) -> Iterable[Node]:
    """Yield children before parents, which makes clade-tip indexing easy."""
    for child in node.children:
        yield from iter_nodes_postorder(child)
    yield node


def index_tree(root: Node) -> dict[str, Node]:
    """Map tip names to nodes and cache the descendant tips below each node."""
    tip_to_node: dict[str, Node] = {}
    for node in iter_nodes_postorder(root):
        if node.is_tip:
            node.tip_names = {node.name}
            if node.name in tip_to_node:
                raise ValueError(f"Duplicate tip label in tree: {node.name}")
            tip_to_node[node.name] = node
        else:
            node.tip_names = set()
            for child in node.children:
                node.tip_names.update(child.tip_names)
    return tip_to_node


def ancestors(node: Node) -> list[Node]:
    """Return the path from a node up to the root."""
    out = []
    while node is not None:
        out.append(node)
        node = node.parent
    return out


def find_mrca(taxa: Iterable[str], tip_to_node: dict[str, Node]) -> tuple[Node | None, list[str], list[str]]:
    """Find the most recent common ancestor of the input taxon labels.

    The function also reports which requested taxa were present or missing in
    the tree, so the audit table can flag placements based on incomplete data.
    """
    present = sorted(taxon for taxon in set(taxa) if taxon in tip_to_node)
    missing = sorted(taxon for taxon in set(taxa) if taxon not in tip_to_node)
    if not present:
        return None, present, missing

    ancestor_lists = [ancestors(tip_to_node[taxon]) for taxon in present]
    common_ids = {node.node_id for node in ancestor_lists[0]}
    for lineage in ancestor_lists[1:]:
        common_ids.intersection_update(node.node_id for node in lineage)

    for node in ancestor_lists[0]:
        if node.node_id in common_ids:
            return node, present, missing
    raise RuntimeError("Could not find MRCA despite present taxa")


# ---------------------------------------------------------------------------
# XLSX helpers
# ---------------------------------------------------------------------------
#
# XLSX files are ZIP archives containing XML. The 1KP supplementary workbooks
# are simple enough that we can read sheet1 directly and avoid optional Python
# packages on HPC/login nodes.

def column_index(cell_ref: str) -> int:
    """Convert an Excel cell reference such as ``C42`` to a 1-based column."""
    letters = re.match(r"[A-Z]+", cell_ref).group(0)
    index = 0
    for char in letters:
        index = index * 26 + (ord(char) - ord("A") + 1)
    return index


def cell_text(cell: ET.Element, shared_strings: list[str]) -> str:
    """Return the displayed text for a worksheet cell."""
    cell_type = cell.attrib.get("t")
    if cell_type == "s":
        value = cell.findtext("a:v", default="", namespaces=XML_NS)
        return shared_strings[int(value)] if value != "" else ""
    if cell_type == "inlineStr":
        return "".join(cell.itertext()).strip()
    return cell.findtext("a:v", default="", namespaces=XML_NS).strip()


def read_xlsx_rows(path: Path, sheet_xml: str = "xl/worksheets/sheet1.xml") -> list[dict[str, str]]:
    """Read the first worksheet of a simple XLSX file as dictionaries.

    Row 1 is interpreted as the header. Later rows are returned as
    ``{header: value}`` dictionaries, preserving empty strings for blank cells.
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
                if not ref:
                    continue
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
    """Split WGM cells that may contain one or more comma/space separated IDs."""
    if not value:
        return []
    ids = []
    for item in re.split(r"[,;]\s*|\s+", value.strip()):
        item = item.strip()
        if item and item.upper() != "NA":
            ids.append(item)
    return ids


def wgm_ids_from_table2(path: Path) -> set[str]:
    """Read the official WGM/WGD IDs from Supplementary Table 2."""
    return {row["WGD ID"] for row in read_xlsx_rows(path) if row.get("WGD ID")}


def wgm_taxa_from_table3(path: Path) -> dict[str, set[str]]:
    """Build ``WGM ID -> supporting 1KP species codes`` from Table 3."""
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
# Outputs and command-line interface
# ---------------------------------------------------------------------------

def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    """Write a placement audit table that explains every tree annotation."""
    path.parent.mkdir(parents=True, exist_ok=True)
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
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def default_output_paths(tree_path: Path) -> tuple[Path, Path]:
    """Derive default output filenames beside the input tree."""
    stem = tree_path.name
    if stem.endswith(".tree"):
        stem = stem[: -len(".tree")]
    out_tree = tree_path.with_name(f"{stem}.wgm_suptab3_mrca.nhx.tree")
    out_tsv = tree_path.with_name(f"{stem}.wgm_suptab3_mrca.tsv")
    return out_tree, out_tsv


def build_arg_parser() -> argparse.ArgumentParser:
    """Create the CLI, keeping project-specific paths as overridable defaults."""
    parser = argparse.ArgumentParser(
        description=(
            "Place 1KP WGM/WGD IDs onto a rooted Newick tree using the MRCA of "
            "supporting taxa listed in Supplementary Table 3."
        )
    )
    parser.add_argument("--tree", type=Path, default=DEFAULT_TREE, help="Rooted 1KP species tree in Newick format.")
    parser.add_argument("--wgm-table", type=Path, default=DEFAULT_TABLE2, help="Supplementary Table 2 WGM/WGD list XLSX.")
    parser.add_argument("--ks-table", type=Path, default=DEFAULT_TABLE3, help="Supplementary Table 3 Ks/WGD history XLSX.")
    parser.add_argument("--out-tree", type=Path, default=None, help="Output NHX-annotated Newick tree.")
    parser.add_argument("--out-tsv", type=Path, default=None, help="Output placement audit TSV.")
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

    valid_wgms = wgm_ids_from_table2(args.wgm_table)
    wgm_to_taxa = wgm_taxa_from_table3(args.ks_table)

    root = NewickParser(args.tree.read_text(encoding="utf-8")).parse()
    tip_to_node = index_tree(root)

    # Each event is placed independently. Multiple WGM/WGD IDs may map to the
    # same branch, so annotations are stored as node_id -> list of IDs.
    annotations: dict[int, list[str]] = defaultdict(list)
    audit_rows: list[dict[str, object]] = []
    skipped = 0
    for wgm_id in sorted(wgm_to_taxa):
        if not args.include_unlisted and wgm_id not in valid_wgms:
            skipped += 1
            continue

        mrca, present, missing = find_mrca(wgm_to_taxa[wgm_id], tip_to_node)
        if mrca is None:
            skipped += 1
            continue

        annotations[mrca.node_id].append(wgm_id)

        # These notes make it easy to review placements that are less precise
        # than a clean "all and only these taxa" clade match.
        exact = len(mrca.tip_names) == len(present)
        if len(present) == 1:
            note = "single sampled support taxon; placed on terminal branch"
        elif exact:
            note = "support taxa exactly match MRCA clade"
        else:
            note = "support taxa are a subset of MRCA clade"

        audit_rows.append(
            {
                "WGM": wgm_id,
                "node_id": mrca.node_id,
                "support_taxa_count": len(wgm_to_taxa[wgm_id]),
                "present_taxa_count": len(present),
                "missing_taxa_count": len(missing),
                "mrca_tip_count": len(mrca.tip_names),
                "exact_support_clade": str(exact).lower(),
                "placement_note": note,
                "support_taxa": ",".join(present),
                "missing_taxa": ",".join(missing),
            }
        )

    out_tree.parent.mkdir(parents=True, exist_ok=True)
    out_tree.write_text(to_newick(root, annotations) + ";\n", encoding="utf-8")
    write_tsv(out_tsv, sorted(audit_rows, key=lambda row: str(row["WGM"])))

    print(f"Wrote annotated tree: {out_tree}")
    print(f"Wrote audit table:    {out_tsv}")
    print(f"Table 2 WGM IDs:      {len(valid_wgms)}")
    print(f"Table 3 WGM IDs:      {len(wgm_to_taxa)}")
    print(f"Annotated WGM IDs:    {len(audit_rows)}")
    print(f"Annotated nodes:      {len(annotations)}")
    if skipped:
        print(f"Skipped WGM IDs:      {skipped}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

