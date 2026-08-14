#!/usr/bin/env python3
"""Score 1KP transcriptomes using BUSCO and TransRate percentile ranks.

The script ranks each transcriptome relative to the full 1KP species table,
then combines the two ranks into one weighted quality score:

    quality_score = 0.8 * busco_percentile + 0.2 * transrate_percentile

Examples:
    python3 score_transcriptome_quality.py --query VYLQ
    python3 score_transcriptome_quality.py --busco-weight 0.7 --transrate-weight 0.3 --query VYLQ
    python3 score_transcriptome_quality.py --top 50
"""

# [Generated via AI, tested by Cecilia]

import argparse
from pathlib import Path
import pandas as pd

BUSCO_COL = "% BUSCOs (complete+fragmented)"
TRANSRATE_COL = "TransRate Score"

ROOT = Path(r"/group/esb/cesen/1kp")
SCORE_TABLE = ROOT / "source_data/1.species_dataset/1kp_paper_2019_suptab1_species.tsv"

def main():
    parser = argparse.ArgumentParser(
        description="Score 1KP transcriptomes using BUSCO and TransRate percentile ranks."
    )
    parser.add_argument(
        "--top",
        type=int,
        default=20,
        help="Number of top-scoring transcriptomes to print",
    )
    parser.add_argument(
        "--query",
        help="Optional case-insensitive text filter for Species or 1KP Index ID",
    )
    parser.add_argument("--busco-weight", type=float, default=0.8)
    parser.add_argument("--transrate-weight", type=float, default=0.2)
    args = parser.parse_args()

    df = pd.read_csv(SCORE_TABLE, sep="\t")

    df[BUSCO_COL] = pd.to_numeric(df[BUSCO_COL], errors="coerce")
    df[TRANSRATE_COL] = pd.to_numeric(df[TRANSRATE_COL], errors="coerce")

    df["busco_percentile"] = df[BUSCO_COL].rank(pct=True)
    df["transrate_percentile"] = df[TRANSRATE_COL].rank(pct=True)
    df["quality_score"] = (
        args.busco_weight * df["busco_percentile"]
        + args.transrate_weight * df["transrate_percentile"]
    )

    out = df
    if args.query:
        query = args.query.lower()
        index_id = df["1KP Index ID"].astype(str).str.lower()
        species = df["Species"].astype(str).str.lower()
        out = df[index_id.str.contains(query, na=False) | species.str.contains(query, na=False)]

    cols = [
        "1KP Index ID",
        "Species",
        BUSCO_COL,
        TRANSRATE_COL,
        "busco_percentile",
        "transrate_percentile",
        "quality_score",
    ]

    out = out.sort_values("quality_score", ascending=False)
    print(out[cols].head(args.top).to_string(index=False))


if __name__ == "__main__":
    main()
