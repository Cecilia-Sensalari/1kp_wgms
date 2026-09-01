#!/bin/bash
#
#SBATCH -p all # partition (queue)
#SBATCH -c 1 # number of cores
#SBATCH --mem 7G # memory pool for all cores
#SBATCH -t 0-72:00 # time (D-HH:MM)
#SBATCH -o cngb.%N.%j.out # STDOUT
#SBATCH -e cngb.%N.%j.err # STDERR

set -euo pipefail

BASE="https://ftp.cngb.org/pub/SciRAID/onekp/assemblies"
OUTDIR="/group/esb/cesen/1kp/source_data/2.transcriptomes/1kp/cngb_transcriptomes"

INDEX_HTML="$OUTDIR/onekp_assemblies_index.html"
ASSEMBLY_LISTING="$OUTDIR/onekp_assemblies_listing.txt"
ASSEMBLY_DIRS="$OUTDIR/onekp_assembly_dirs.txt"
URL_LIST="$OUTDIR/translated_nucleotides_urls.txt"
MISSING_LOG="$OUTDIR/missing_transcriptomes_cngb.tsv"

mkdir -p "$OUTDIR"

echo "1. Downloading top-level CNGB assemblies listing..."
wget -qO "$INDEX_HTML" "$BASE/"

echo "2. Extracting and HTML-decoding linked entries..."

python3 - "$INDEX_HTML" > "$ASSEMBLY_LISTING" <<'PY'
import sys
import html
from html.parser import HTMLParser

class LinkParser(HTMLParser):
    def handle_starttag(self, tag, attrs):
        if tag.lower() != "a":
            return
        for key, value in attrs:
            if key.lower() == "href" and value:
                print(html.unescape(value))

with open(sys.argv[1], encoding="utf-8", errors="replace") as handle:
    parser = LinkParser()
    parser.feed(handle.read())
PY

grep -v '^\.\./\?$' "$ASSEMBLY_LISTING" \
    | grep -v '^$' \
    | sort -u \
    > "$ASSEMBLY_LISTING.tmp"

mv "$ASSEMBLY_LISTING.tmp" "$ASSEMBLY_LISTING"

echo "Saved full decoded listing to:"
echo "  $ASSEMBLY_LISTING"

echo "3. Extracting assembly directory names..."

grep '/$' "$ASSEMBLY_LISTING" \
    | sed 's#/$##' \
    | grep -v '^\.\.$' \
    | grep -E '^[A-Za-z0-9]+-.+$' \
    | sort -u \
    > "$ASSEMBLY_DIRS"

N_DIRS=$(wc -l < "$ASSEMBLY_DIRS")
echo "Found $N_DIRS assembly directories."

echo "4. Generating expected translated-nucleotides FASTA URLs..."

: > "$URL_LIST"

while read -r dir; do
    id="${dir%%-*}"
    echo "$BASE/$dir/${id}-translated-nucleotides.fa.gz"
done < "$ASSEMBLY_DIRS" > "$URL_LIST"

echo "Saved URL list to:"
echo "  $URL_LIST"

echo "5. Initialising missing transcriptome log..."

printf "id\tdirectory\turl\treason\n" > "$MISSING_LOG"

echo "6. Downloading and gunzipping files sequentially..."

while read -r dir; do
    id="${dir%%-*}"

    gz_file="$OUTDIR/${id}-translated-nucleotides.fa.gz"
    fa_file="$OUTDIR/${id}-translated-nucleotides.fa"
    url="$BASE/$dir/${id}-translated-nucleotides.fa.gz"

    if [[ -s "$fa_file" ]]; then
        echo "Skipping $id: already decompressed."
        continue
    fi

    echo "Downloading $id from directory:"
    echo "  $dir"

    if wget \
        --continue \
        --tries=3 \
        --timeout=60 \
        --waitretry=10 \
        --directory-prefix="$OUTDIR" \
        "$url"
    then
        if [[ ! -s "$gz_file" ]]; then
            echo "Missing or empty downloaded file for $id."
            printf "%s\t%s\t%s\t%s\n" "$id" "$dir" "$url" "downloaded_gz_missing_or_empty" >> "$MISSING_LOG"
            continue
        fi

        echo "Checking gzip integrity for $id..."

        if gunzip -t "$gz_file"; then
            echo "Gunzipping $id..."
            gunzip -f "$gz_file"
        else
            echo "Invalid gzip file for $id; removing corrupted file."
            rm -f "$gz_file"
            printf "%s\t%s\t%s\t%s\n" "$id" "$dir" "$url" "invalid_gzip" >> "$MISSING_LOG"
            continue
        fi

    else
        wget_exit_code=$?
        echo "Download failed for $id; continuing with next transcriptome."
        printf "%s\t%s\t%s\twget_failed_exit_%s\n" "$id" "$dir" "$url" "$wget_exit_code" >> "$MISSING_LOG"
        continue
    fi

done < "$ASSEMBLY_DIRS"

echo "Done."
echo "Downloaded and decompressed FASTA files are in:"
echo "  $OUTDIR/"

echo "Missing/failed transcriptomes were logged in:"
echo "  $MISSING_LOG"

N_MISSING=$(($(wc -l < "$MISSING_LOG") - 1))
echo "Number of missing/failed transcriptomes: $N_MISSING"