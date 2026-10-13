#!/usr/bin/env python3
"""Collect the BUSCO v3 results of run_busco_genomes.sh into one table.

Reads OUT_DIR/output_busco/<species>/run_<species>/short_summary_<species>.txt for every species and writes a TSV with
  species_id, busco_complete_fragmented (%), complete (%), fragmented (%), missing (%), n_buscos
where busco_complete_fragmented = 100 - missing is the same quantity as the 1KP table column
"% BUSCOs (complete+fragmented)". Species without a result are listed on stderr.

Usage:   collect_busco_scores.py [--busco-dir DIR] [--out FILE]
"""

import argparse
import re
import sys
from pathlib import Path

DEFAULT_DIR = Path("/group/esb/cesen/1kp/source_data/2.genomes/1kp/genome_busco/output_busco")
DEFAULT_MANIFEST = Path("/group/esb/cesen/1kp/source_data/2.genomes/1kp/genome_species_manifest.tsv")
SUMMARY = re.compile(
    r"C:(?P<c>[\d.]+)%\[S:[\d.]+%,D:[\d.]+%\],F:(?P<f>[\d.]+)%,M:(?P<m>[\d.]+)%,n:(?P<n>\d+)"
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--busco-dir", type=Path, default=DEFAULT_DIR)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST, help="genome manifest TSV, species_id in column 1")
    parser.add_argument("--out", type=Path, default=None, help="default: BUSCO_DIR/genome_busco_scores.tsv")
    args = parser.parse_args()
    out = args.out or args.busco_dir / "genome_busco_scores.tsv"
    if out.exists():
        sys.exit(f"ERROR: [{out}] already exists, not overwriting it")

    rows, missing = [], []
    with args.manifest.open() as handle:
        species_ids = [line.split("\t")[0] for line in handle.read().splitlines()[1:] if line]
    for species in species_ids:
        species_dir = args.busco_dir / species
        summary_file = species_dir / f"run_{species}" / f"short_summary_{species}.txt"
        match = SUMMARY.search(summary_file.read_text()) if summary_file.is_file() else None
        if not match:
            missing.append(species)
            continue
        complete, fragmented, absent = (float(match[key]) for key in ("c", "f", "m"))
        rows.append([species, f"{100 - absent:.1f}", f"{complete:.1f}", f"{fragmented:.1f}", f"{absent:.1f}", match["n"]])

    with out.open("w") as handle:
        handle.write("species_id\tbusco_complete_fragmented\tcomplete\tfragmented\tmissing\tn_buscos\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")
    print(f"{len(rows)} species written to {out}")
    if missing:
        print(f"WARNING: no BUSCO result for: {', '.join(missing)}", file=sys.stderr)


if __name__ == "__main__":
    main()
