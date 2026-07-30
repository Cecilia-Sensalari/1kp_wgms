#!/usr/bin/env python3
"""

List direct children or descendants of a chosen rank in 1KP taxonomy,
using as input the taxonomy file based on the 1KP species list.

[Generate via AI, tested by Cecilia]

python list_direct_subclades.py species "Oenothera elata" subspecies
Oenothera elata subsp. elata
Oenothera elata subsp. hookeri

"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


DEFAULT_TAXONOMY = (
    Path(__file__).resolve().parent.parent.parent
    / "source_data"
    / "1.species_dataset"
    / "1kp_paper_2019_suptab1_species_ncbi_taxonomy.csv"
)


def lineage(row: dict[str, str]) -> list[tuple[str, str]]:
    """Return (rank, name) pairs, ignoring malformed lineage fields."""
    names = [value.strip() for value in row.get("lineage_names", "").split(";")]
    ranks = [value.strip() for value in row.get("lineage_ranks", "").split(";")]
    if len(names) != len(ranks):
        return []
    return list(zip(ranks, names))


def subclades(
    csv_path: Path,
    clade_type: str,
    clade_name: str,
    target_type: str | None = None,
) -> list[tuple[str, str]]:
    """Find direct children, or all descendants having ``target_type``."""
    wanted_rank = clade_type.strip().casefold()
    wanted_name = clade_name.strip().casefold()
    target_rank = target_type.strip().casefold() if target_type else None
    children: set[tuple[str, str]] = set()

    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"lineage_names", "lineage_ranks"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"missing CSV column(s): {', '.join(sorted(missing))}")

        for row in reader:
            nodes = lineage(row)
            for index, (rank, name) in enumerate(nodes):
                if rank.casefold() == wanted_rank and name.casefold() == wanted_name:
                    descendants = nodes[index + 1 :]
                    if target_rank is None:
                        descendants = descendants[:1]
                    for child_rank, child_name in descendants:
                        if target_rank is None or child_rank.casefold() == target_rank:
                            children.add((child_rank, child_name))

    return sorted(children, key=lambda child: (child[1].casefold(), child[0].casefold()))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "List descendants of a clade in the 1KP dataset. If TARGET_TYPE is "
            "omitted, list only its direct children. Matching is case-insensitive."
        )
    )
    parser.add_argument("clade_type", help='rank of the clade, for example "genus"')
    parser.add_argument("clade_name", help='name of the clade, for example "Draba"')
    parser.add_argument(
        "target_type",
        nargs="?",
        help=(
            'rank to return, for example "genus"; omit it to return direct children'
        ),
    )
    parser.add_argument(
        "--taxonomy-file",
        type=Path,
        default=DEFAULT_TAXONOMY,
        help=f"taxonomy CSV (default: {DEFAULT_TAXONOMY})",
    )
    parser.add_argument(
        "--show-ranks",
        action="store_true",
        help="print each result's rank after its name",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        children = subclades(
            args.taxonomy_file,
            args.clade_type,
            args.clade_name,
            args.target_type,
        )
    except (OSError, csv.Error, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    if not children:
        target = f" of rank {args.target_type!r}" if args.target_type else ""
        print(f"No subclades{target} found for "
              f"{args.clade_type} {args.clade_name!r}.", file=sys.stderr)
        return 1

    for rank, name in children:
        print(f"{name}\t{rank}" if args.show_ranks else name)
    print(f"Retrieved subclades: {len(children)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
