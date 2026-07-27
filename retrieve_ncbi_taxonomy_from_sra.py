#!/usr/bin/env python3
"""Retrieve rank-labelled NCBI taxonomy for the 1KP species workbook.

[Generated via AI, tested by Cecilia]

The script uses the SRA accessions in ``NCBI SRA Accession IDs`` rather than
matching the reported species names.  It first resolves each accession in the
NCBI SRA database, extracts the organism TaxID from the SRA record, and then
retrieves the corresponding ranked lineage from the NCBI Taxonomy database.

Only the Python standard library is required, including for reading XLSX files.
Results are cached in JSON, so an interrupted run can be restarted without
repeating completed NCBI requests.

Some workbook rows contain more than one SRA accession.  Normally all
accessions for one 1KP sample resolve to the same organism TaxID.  If they
resolve to two or more distinct TaxIDs, however, the output status is
``conflicting_taxids``.  In that case the script reports every resolved TaxID
in the semicolon-separated ``taxid`` column but deliberately leaves the
ranked taxonomy columns empty: choosing one lineage automatically could hide
a sample mix-up, contamination, or inconsistent NCBI records.  The cache's
``accessions`` section records the accession-to-TaxID mapping for diagnosis.
The user has therefore to manually choose the NCBI SRA Accession ID that matches
the species name from the 1KP paper.
E.g.
KJAA,Mollugo pentaphylla,ERS3670308;ERS1829300
Because the two accessions that resolve to different TaxIDs, the output
will have status ``conflicting_taxids`` and the user will have to choose the correct accession
to use for the second row:
KJAA,Mollugo pentaphylla,ERS1829300

Example::

    python code/retrieve_ncbi_taxonomy_from_sra.py \
        --email your.name@example.org

NCBI asks API users to identify themselves with an email address.  Supplying an
NCBI API key through ``--api-key`` or the ``NCBI_API_KEY`` environment variable
raises the E-utilities request limit from 3 to 10 requests per second.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Sequence
from xml.etree import ElementTree


PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    PROJECT / "source_data/1.species_dataset/1kp_paper_2019_suptab1_species.xlsx"
)
DEFAULT_OUTPUT = (
    PROJECT / "source_data/1.species_dataset/1kp_paper_2019_suptab1_species_ncbi_taxonomy.csv"
)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL_NAME = "onekp_taxonomy"
ACCESSION_PATTERN = re.compile(r"\b[A-Z]{3}\d+\b")
SEARCH_BATCH_SIZE = 40
SUMMARY_BATCH_SIZE = 200
TAXONOMY_BATCH_SIZE = 200
STANDARD_RANKS = (
    "superkingdom",
    "kingdom",
    "subkingdom",
    "phylum",
    "subphylum",
    "class",
    "subclass",
    "order",
    "suborder",
    "family",
    "subfamily",
    "tribe",
    "genus",
    "species",
    "subspecies",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Resolve 1KP SRA accessions to rank-labelled NCBI taxonomy."
    )
    parser.add_argument("-i", "--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("-o", "--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--cache",
        type=Path,
        help="JSON cache path (default: OUTPUT with .cache.json appended)",
    )
    parser.add_argument(
        "--email",
        default=os.environ.get("NCBI_EMAIL"),
        help="Contact email sent to NCBI (or set NCBI_EMAIL)",
    )
    parser.add_argument(
        "--api-key",
        default=os.environ.get("NCBI_API_KEY"),
        help="Optional NCBI API key (or set NCBI_API_KEY)",
    )
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Ignore cached accession and taxonomy results",
    )
    args = parser.parse_args()
    if not args.email:
        parser.error("--email is required (or set NCBI_EMAIL)")
    if args.cache is None:
        args.cache = Path(str(args.output) + ".cache.json")
    return args


def excel_column_number(cell_reference: str) -> int:
    """Convert the letters in an XLSX cell reference (e.g. AB7) to 1-based index."""
    match = re.match(r"[A-Z]+", cell_reference)
    if not match:
        raise ValueError(f"Invalid XLSX cell reference: {cell_reference}")
    number = 0
    for character in match.group():
        number = number * 26 + ord(character) - ord("A") + 1
    return number


def read_xlsx_rows(path: Path) -> List[Dict[str, str]]:
    """Read named columns from the first worksheet using the XLSX XML files."""
    namespace = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(path) as archive:
        shared_root = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
        shared_strings = [
            "".join(part.text or "" for part in item.findall(".//x:t", namespace))
            for item in shared_root.findall("x:si", namespace)
        ]
        sheet_root = ElementTree.fromstring(
            archive.read("xl/worksheets/sheet1.xml")
        )

    matrix: List[Dict[int, str]] = []
    for xml_row in sheet_root.findall(".//x:sheetData/x:row", namespace):
        row: Dict[int, str] = {}
        for cell in xml_row.findall("x:c", namespace):
            column = excel_column_number(cell.attrib["r"])
            value_node = cell.find("x:v", namespace)
            if value_node is None:
                inline = cell.find("x:is/x:t", namespace)
                if inline is not None:
                    row[column] = inline.text or ""
                continue
            value = value_node.text or ""
            if cell.attrib.get("t") == "s":
                value = shared_strings[int(value)]
            row[column] = value
        matrix.append(row)

    if not matrix:
        return []
    header = {column: value.strip() for column, value in matrix[0].items() if value}
    required = {"1KP Index ID", "Species", "NCBI SRA Accession IDs"}
    if not required.issubset(header.values()):
        missing = sorted(required - set(header.values()))
        raise ValueError(f"Missing expected XLSX column(s): {', '.join(missing)}")

    return [
        {name: row.get(column, "").strip() for column, name in header.items()}
        for row in matrix[1:]
    ]


def batches(values: Sequence[str], size: int) -> Iterator[Sequence[str]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


class EntrezClient:
    """Small E-utilities client with NCBI-compliant throttling and retries."""

    def __init__(self, email: str, api_key: str | None = None) -> None:
        self.common = {"tool": TOOL_NAME, "email": email}
        if api_key:
            self.common["api_key"] = api_key
        self.minimum_interval = 0.11 if api_key else 0.34
        self.last_request = 0.0

    def post(self, endpoint: str, parameters: Mapping[str, str]) -> bytes:
        payload = dict(self.common)
        payload.update(parameters)
        request = urllib.request.Request(
            f"{EUTILS}/{endpoint}",
            data=urllib.parse.urlencode(payload).encode("ascii"),
            headers={"User-Agent": f"{TOOL_NAME}/1.0 ({self.common['email']})"},
        )
        for attempt in range(6):
            delay = self.minimum_interval - (time.monotonic() - self.last_request)
            if delay > 0:
                time.sleep(delay)
            try:
                with urllib.request.urlopen(request, timeout=60) as response:
                    body = response.read()
                self.last_request = time.monotonic()
                return body
            except urllib.error.HTTPError as error:
                self.last_request = time.monotonic()
                if error.code not in {429, 500, 502, 503, 504} or attempt == 5:
                    raise
            except urllib.error.URLError:
                self.last_request = time.monotonic()
                if attempt == 5:
                    raise
            time.sleep(2**attempt)
        raise RuntimeError("Unreachable retry state")


def load_cache(path: Path, refresh: bool) -> MutableMapping[str, object]:
    if not refresh and path.exists():
        with path.open(encoding="utf-8") as handle:
            cache = json.load(handle)
        if cache.get("version") == 1:
            return cache
    return {"version": 1, "accessions": {}, "taxonomy": {}}


def save_cache(path: Path, cache: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(cache, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(path)


def resolve_sra_accessions(
    client: EntrezClient,
    accessions: Sequence[str],
    cache: MutableMapping[str, object],
    cache_path: Path,
) -> None:
    """Populate cache['accessions'] with SRA UID and organism TaxID values."""
    accession_cache = cache["accessions"]
    assert isinstance(accession_cache, dict)
    unresolved = sorted(set(accessions) - set(accession_cache))

    for number, group in enumerate(batches(unresolved, SEARCH_BATCH_SIZE), 1):
        query = " OR ".join(f"{accession}[ACCN]" for accession in group)
        search = json.loads(
            client.post(
                "esearch.fcgi",
                {"db": "sra", "term": query, "retmode": "json", "retmax": "10000"},
            )
        )
        uids = search.get("esearchresult", {}).get("idlist", [])
        summaries: List[Mapping[str, object]] = []
        for uid_group in batches(uids, SUMMARY_BATCH_SIZE):
            response = json.loads(
                client.post(
                    "esummary.fcgi",
                    {"db": "sra", "id": ",".join(uid_group), "retmode": "json"},
                )
            )
            result = response.get("result", {})
            summaries.extend(result[uid] for uid in uid_group if uid in result)

        for accession in group:
            matches = []
            for summary in summaries:
                searchable = " ".join(str(value) for value in summary.values())
                if re.search(rf"(?<![A-Z0-9]){re.escape(accession)}(?![A-Z0-9])", searchable, re.I):
                    organism = re.search(
                        r"<Organism\b[^>]*\btaxid=\"(\d+)\"", searchable, re.I
                    )
                    if organism:
                        matches.append(
                            {"uid": str(summary.get("uid", "")), "taxid": organism.group(1)}
                        )
            unique = {(match["uid"], match["taxid"]) for match in matches}
            accession_cache[accession] = [
                {"uid": uid, "taxid": taxid} for uid, taxid in sorted(unique)
            ]

        save_cache(cache_path, cache)
        print(
            f"Resolved SRA batch {number}: {min(number * SEARCH_BATCH_SIZE, len(unresolved))}"
            f"/{len(unresolved)} accessions",
            file=sys.stderr,
        )


def parse_taxon(taxon: ElementTree.Element) -> Dict[str, object]:
    taxid = taxon.findtext("TaxId", "")
    scientific_name = taxon.findtext("ScientificName", "")
    rank = taxon.findtext("Rank", "no rank")
    lineage_nodes = list(taxon.findall("./LineageEx/Taxon"))
    lineage_nodes.append(taxon)

    ranked: Dict[str, str] = {}
    lineage_names: List[str] = []
    lineage_taxids: List[str] = []
    lineage_ranks: List[str] = []
    for node in lineage_nodes:
        node_name = node.findtext("ScientificName", "")
        node_taxid = node.findtext("TaxId", "")
        node_rank = node.findtext("Rank", "no rank")
        lineage_names.append(node_name)
        lineage_taxids.append(node_taxid)
        lineage_ranks.append(node_rank)
        if node_rank in STANDARD_RANKS:
            ranked[node_rank] = node_name

    return {
        "taxid": taxid,
        "scientific_name": scientific_name,
        "rank": rank,
        "ranked": ranked,
        "lineage_names": lineage_names,
        "lineage_taxids": lineage_taxids,
        "lineage_ranks": lineage_ranks,
    }


def retrieve_taxonomy(
    client: EntrezClient,
    taxids: Sequence[str],
    cache: MutableMapping[str, object],
    cache_path: Path,
) -> None:
    taxonomy_cache = cache["taxonomy"]
    assert isinstance(taxonomy_cache, dict)
    unresolved = sorted(set(taxids) - set(taxonomy_cache), key=int)
    for number, group in enumerate(batches(unresolved, TAXONOMY_BATCH_SIZE), 1):
        xml = client.post(
            "efetch.fcgi",
            {"db": "taxonomy", "id": ",".join(group), "retmode": "xml"},
        )
        root = ElementTree.fromstring(xml)
        returned = set()
        for taxon in root.findall("./Taxon"):
            parsed = parse_taxon(taxon)
            returned.add(str(parsed["taxid"]))
            taxonomy_cache[str(parsed["taxid"])] = parsed
        for missing in set(group) - returned:
            taxonomy_cache[missing] = None
        save_cache(cache_path, cache)
        print(
            f"Retrieved taxonomy batch {number}: "
            f"{min(number * TAXONOMY_BATCH_SIZE, len(unresolved))}/{len(unresolved)} TaxIDs",
            file=sys.stderr,
        )


def sample_result(
    source: Mapping[str, str], cache: Mapping[str, object]
) -> Dict[str, str]:
    """Combine all SRA results for one workbook row into one output row.

    A 1KP sample can have multiple SRA accessions, for example when its reads
    were deposited in more than one experiment.  Repeated accessions and
    repeated SRA records do not constitute a conflict: the set operation below
    collapses them, and the sample is accepted when exactly one distinct TaxID
    remains.

    ``conflicting_taxids`` means that the accessions collectively resolve to
    more than one *distinct* organism TaxID.  The semicolon-separated ``taxid``
    and ``sra_uids`` fields are retained, but no lineage is selected and all
    rank columns remain empty.  To see which input accession produced each
    TaxID, inspect ``accessions`` in the JSON cache.  This conservative behavior
    prevents an arbitrary first match from masking inconsistent source data.
    """
    accession_cache = cache["accessions"]
    taxonomy_cache = cache["taxonomy"]
    assert isinstance(accession_cache, dict) and isinstance(taxonomy_cache, dict)
    accession_text = source.get("NCBI SRA Accession IDs", "")
    accessions = ACCESSION_PATTERN.findall(accession_text.upper())
    result = {
        "sample_id": source.get("1KP Index ID", ""),
        "reported_species": source.get("Species", ""),
        "sra_accessions": ";".join(accessions),
        "sra_uids": "",
        "taxid": "",
        "ncbi_scientific_name": "",
        "ncbi_rank": "",
        **{rank: "" for rank in STANDARD_RANKS},
        "lineage_names": "",
        "lineage_taxids": "",
        "lineage_ranks": "",
        "status": "",
    }
    if not accessions:
        result["status"] = "no_sra_accession"
        return result

    matches = [match for accession in accessions for match in accession_cache.get(accession, [])]
    if not matches:
        result["status"] = "sra_not_found"
        return result
    result["sra_uids"] = ";".join(sorted({match["uid"] for match in matches}))

    # Reduce all records from all of this sample's accessions to distinct
    # organism TaxIDs. Multiple records are harmless when their TaxID agrees.
    taxids = sorted({match["taxid"] for match in matches}, key=int)
    result["taxid"] = ";".join(taxids)
    if len(taxids) > 1:
        # Do not guess which organism is correct.  Keeping the competing TaxIDs
        # while leaving lineage columns blank makes the conflict visible and
        # allows it to be reviewed against the accession-level JSON cache.
        result["status"] = "conflicting_taxids"
        return result

    taxonomy = taxonomy_cache.get(taxids[0])
    if not taxonomy:
        result["status"] = "taxonomy_not_found"
        return result
    result["ncbi_scientific_name"] = str(taxonomy["scientific_name"])
    result["ncbi_rank"] = str(taxonomy["rank"])
    result.update({key: str(value) for key, value in taxonomy["ranked"].items()})
    # Rebuild ranked fields from the cached full lineage as well.  Taxonomy
    # records cached by older versions of this script have no ``subspecies``
    # entry in their ``ranked`` mapping, although their lineage arrays already
    # contain it.  This makes old caches compatible without another NCBI call.
    result.update(
        {
            str(rank): str(name)
            for rank, name in zip(
                taxonomy["lineage_ranks"], taxonomy["lineage_names"]
            )
            if rank in STANDARD_RANKS
        }
    )
    result["lineage_names"] = ";".join(taxonomy["lineage_names"])
    result["lineage_taxids"] = ";".join(taxonomy["lineage_taxids"])
    result["lineage_ranks"] = ";".join(taxonomy["lineage_ranks"])
    result["status"] = "ok"
    return result


def main() -> None:
    args = parse_args()
    if not args.input.exists():
        raise SystemExit(f"Input workbook not found: {args.input}")
    source_rows = read_xlsx_rows(args.input)
    all_accessions = sorted(
        {
            accession
            for row in source_rows
            for accession in ACCESSION_PATTERN.findall(
                row.get("NCBI SRA Accession IDs", "").upper()
            )
        }
    )
    cache = load_cache(args.cache, args.refresh)
    client = EntrezClient(args.email, args.api_key)
    resolve_sra_accessions(client, all_accessions, cache, args.cache)

    accession_cache = cache["accessions"]
    assert isinstance(accession_cache, dict)
    taxids = sorted(
        {
            match["taxid"]
            for matches in accession_cache.values()
            for match in matches
        },
        key=int,
    )
    retrieve_taxonomy(client, taxids, cache, args.cache)

    output_rows = [sample_result(row, cache) for row in source_rows]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(output_rows[0]) if output_rows else []
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)
    counts: Dict[str, int] = {}
    for row in output_rows:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    print(f"Wrote {len(output_rows)} rows to {args.output}")
    print("Statuses: " + ", ".join(f"{key}={value}" for key, value in sorted(counts.items())))


if __name__ == "__main__":
    main()
