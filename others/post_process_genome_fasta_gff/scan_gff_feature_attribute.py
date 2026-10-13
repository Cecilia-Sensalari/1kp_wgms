#!/usr/bin/env python3
"""Find, for each genome species, which GFF feature and attribute identify the genes of its CDS FASTA.

ksrates' collinearity analysis needs two settings per focal species, gff_feature (the feature type in
column 3 of the GFF, e.g. mRNA) and gff_attribute (the attribute in column 9, e.g. Name). They are used by
wgd's gff_parser (wgd_ksrates/colinearity.py), which makes one gene-list entry for every GFF line of that
feature, with the whole value of that attribute as the gene ID and the coordinates of that same line.
So a valid pair has to identify every gene exactly once:
  - the attribute values of that feature are the sequence IDs of the species' CDS FASTA file, and
  - there is only one such GFF line per ID (univocal). For example the exon lines of a gene all carry the same
    Parent, so exon/Parent would put every gene in the list once per exon: it matches the IDs but is wrong.
Pairs that are not univocal are discarded, they are never reported.
The goal of the pair is to give the coordinates of the whole transcript. Only features whose line spans the
whole transcript are therefore considered (WHOLE_TRANSCRIPT_FEATURES: mRNA, gene, transcript); features that
cover only a part of it, such as exon or CDS segments or UTRs, are never candidates.

This script follows the parser's own rules: attribute names are compared case-sensitively and are reported with
the spelling used in the GFF (give exactly that spelling to ksrates for attributes other than ID, Name and Parent);
the value of an attribute is taken whole (it is not split on commas); lines whose strand is not + or - are
ignored; a line that does not have 9 columns makes the parser fail (reported as BAD_LINES).

For every species of the genome manifest it tests all (feature, attribute) combinations present in the GFF and
reports the three best univocal ones. The "coverage" of a candidate is the fraction of the FASTA IDs found among
its attribute values (1.000 = all of them). The FASTA files hold only the primary transcripts, so GFF entries
without a FASTA ID (secondary transcripts) are expected and ignored. Equally good candidates are ordered mRNA,
gene, transcript, then by name. Always look at the species with something in the "note" column:
  NO_MATCH          no attribute value of any feature is a FASTA ID
  NOT_UNIVOCAL      FASTA IDs are found, but only in pairs with several GFF lines per ID: no valid pair
  NO_TRANSCRIPT_FEATURE(feature/attribute)  univocal pairs exist only for features that do not span the whole
                    transcript, e.g. CDS: the best of them is shown in the brackets, but it is not a candidate
  LOW_COVERAGE      the best candidate finds less than 90% of the FASTA IDs
  TIE               the 2nd candidate is as good as the best one (both valid: choose one)
  KSRATES_CHANGES_SPELLING  ksrates would rewrite the best candidate's spelling (e.g. id -> ID, mrna -> mRNA),
                    which then no longer matches the GFF: it cannot be used through the config as it is
  BAD_LINES(n)      n GFF lines do not have 9 columns: gff_parser would stop with an error
  CDS_FILE_MISSING / GFF_FILE_MISSING

Run with a Python 3 that has only the standard library, e.g.:
    python3 scan_gff_feature_attribute.py
    python3 scan_gff_feature_attribute.py --manifest MANIFEST.tsv --out RESULT.tsv

Takes a few minutes for ~30 genomes (it reads every GFF once).
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

DEFAULT_MANIFEST = Path("/group/esb/cesen/1kp/source_data/2.genomes/1kp/genome_species_manifest.tsv")
DEFAULT_OUT = Path("/group/esb/cesen/1kp/source_data/2.genomes/1kp/gff_feature_attribute_scan.tsv")

# Coverage below this is flagged for a manual check
LOW_COVERAGE = 0.9

# Features whose GFF line spans the whole transcript, in order of preference among equally good candidates
WHOLE_TRANSCRIPT_FEATURES = ["mRNA", "gene", "transcript"]

# ksrates rewrites the spelling of these common names before giving them to gff_parser (see get_feature and
# get_attribute in ksrates/fc_configfile.py); anything else is passed unchanged
KSRATES_FEATURE_SPELLING = {
    "gene": "gene", "mrna": "mRNA", "cds": "CDS", "cdna_match": "cDNA_match",
    "transcript": "transcript", "region": "region", "exon": "exon",
}
KSRATES_ATTRIBUTE_SPELLING = {"id": "ID", "name": "Name", "parent": "Parent"}

N_CANDIDATES = 3
CANDIDATE_FIELDS = ["feature", "attribute", "coverage"]
COLUMNS = (
    ["species_id", "cds_file", "gff_file", "n_fasta_ids"]
    + [f"{prefix}{field}" for prefix in ("best_", "2nd_", "3rd_") for field in CANDIDATE_FIELDS]
    + ["note"]
)

Pair = Tuple[str, str]


def read_fasta_ids(path: Path) -> Set[str]:
    """First word of every header line of a FASTA file."""
    ids = set()
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith(">") and len(line.strip()) > 1:
                ids.add(line[1:].split(None, 1)[0])
    return ids


def parse_attributes(text: str) -> Dict[str, str]:
    """Column 9 as gff_parser reads it: 'key=value' fields separated by ';'; fields with anything else are dropped."""
    attributes = {}
    for field in text.split(";"):
        parts = field.split("=")
        if len(parts) == 2:
            attributes[parts[0]] = parts[1]
    return attributes


def scan_gff(path: Path, fasta_ids: Set[str]) -> Tuple[Dict[Pair, Set[str]], Counter, int]:
    """
    For every (feature, attribute) pair of the GFF, collects the FASTA IDs found among its values.

    :return: (FASTA IDs found per pair, number of GFF lines of the pair that carry a FASTA ID,
              number of lines without 9 columns)
    """
    matched: Dict[Pair, Set[str]] = defaultdict(set)
    lines_with_id: Counter = Counter()
    bad_lines = 0
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or line == "\n":
                continue
            columns = line.strip().split("\t")
            if len(columns) != 9:
                bad_lines += 1
                continue
            if columns[6] not in ("+", "-"):  # gff_parser skips lines without strand
                continue
            feature = columns[2]
            for key, value in parse_attributes(columns[8]).items():
                if value in fasta_ids:
                    matched[(feature, key)].add(value)
                    lines_with_id[(feature, key)] += 1
    return matched, lines_with_id, bad_lines


def feature_rank(feature: str) -> int:
    return WHOLE_TRANSCRIPT_FEATURES.index(feature)


def ksrates_changes_spelling(feature: str, attribute: str) -> bool:
    """True if ksrates would rewrite this feature or attribute spelling, so the GFF's own spelling can never match."""
    return (KSRATES_FEATURE_SPELLING.get(feature.lower(), feature) != feature
            or KSRATES_ATTRIBUTE_SPELLING.get(attribute.lower(), attribute) != attribute)


def scan_species(species_id: str, cds_path: str, gff_path: str) -> List[str]:
    """One output row for a species."""
    # File names (not paths), to find the files in the genome directory
    first_columns = [species_id, Path(cds_path).name if cds_path else "", Path(gff_path).name if gff_path else ""]
    if not cds_path or not Path(cds_path).is_file():
        return first_columns + [""] * (len(COLUMNS) - 4) + ["CDS_FILE_MISSING"]
    if not gff_path or not Path(gff_path).is_file():
        return first_columns + [""] * (len(COLUMNS) - 4) + ["GFF_FILE_MISSING"]

    fasta_ids = read_fasta_ids(Path(cds_path))
    matched, lines_with_id, bad_lines = scan_gff(Path(gff_path), fasta_ids)
    n_ids = len(fasta_ids)

    # Univocal pairs only: each FASTA ID found on exactly one GFF line of the feature
    univocal = [pair for pair in matched if lines_with_id[pair] == len(matched[pair])]
    # Candidates: univocal pairs of a feature that spans the whole transcript. Most IDs found first, then the
    # feature preference, then the name (same result on every run)
    candidates = [pair for pair in univocal if pair[0] in WHOLE_TRANSCRIPT_FEATURES]
    ranking = sorted(candidates, key=lambda pair: (-len(matched[pair]), feature_rank(pair[0]), pair[1]))

    row = first_columns + [str(n_ids)]
    for position in range(N_CANDIDATES):
        if position < len(ranking):
            pair = ranking[position]
            row.extend([pair[0], pair[1], f"{len(matched[pair]) / n_ids:.3f}"])
        else:
            row.extend(["-", "-", "0"])

    notes = []
    if not matched:
        notes.append("NO_MATCH")
    elif not univocal:
        notes.append("NOT_UNIVOCAL")
    elif not ranking:
        # Univocal pairs exist, but only of features covering part of the transcript (e.g. CDS)
        other = min(univocal, key=lambda pair: (-len(matched[pair]), pair))
        notes.append(f"NO_TRANSCRIPT_FEATURE({other[0]}/{other[1]})")
    else:
        best = ranking[0]
        if len(matched[best]) / n_ids < LOW_COVERAGE:
            notes.append("LOW_COVERAGE")
        if ksrates_changes_spelling(best[0], best[1]):
            notes.append("KSRATES_CHANGES_SPELLING")
        if len(ranking) > 1 and len(matched[ranking[1]]) == len(matched[best]):
            notes.append("TIE")
    if bad_lines:
        notes.append(f"BAD_LINES({bad_lines})")
    row.append(";".join(notes))
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST,
                        help="genome manifest TSV with columns species_id, cds_filepath, gff_filepath")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="output TSV")
    args = parser.parse_args()

    with args.manifest.open(newline="") as handle:
        manifest_rows = list(csv.DictReader(handle, delimiter="\t"))

    with args.out.open("w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        for manifest_row in manifest_rows:
            row = scan_species(
                manifest_row["species_id"],
                manifest_row.get("cds_filepath", ""),
                manifest_row.get("gff_filepath", ""),
            )
            assert len(row) == len(COLUMNS), row
            writer.writerow(row)
            out_handle.flush()
            column = {name: row[position] for position, name in enumerate(COLUMNS)}
            best = (f"{column['best_feature']}/{column['best_attribute']} (coverage {column['best_coverage']})"
                    if column["best_feature"] not in ("", "-") else "-")
            print(f"{column['species_id']}: {best} {column['note']}", flush=True)
    print(f"Written to {args.out}")


if __name__ == "__main__":
    main()
