#!/usr/bin/env bash
# Rename the sequences of an NCBI-style CDS FASTA file to their locus_tag, e.g.
#   >lcl|BFEA01000001.1_cds_GBG58599.1_1 [locus_tag=CBR_g1] [protein=hypothetical protein] ...
# becomes
#   >CBR_g1
# so that the sequence names are the values that the GFF (e.g. gene Name=CBR_g1) uses to identify the genes.
#
# A locus tag can have several CDS records (alternative isoforms, different protein_id and length), but a
# sequence name must be unique. Only one record is kept per locus tag, the LONGEST (the first one if equally
# long), i.e. the primary transcript. Everything that is dropped is listed in a log file.
#
# Usage:   rename_fasta_to_locus_tag.sh INPUT.fna [OUTPUT.fna]
# Default OUTPUT: INPUT with ".locus_tag" before the extension, in the same directory.
# The input file is never changed, and an existing output file is never overwritten. Besides OUTPUT it writes
#   OUTPUT.mapping.tsv   new name, original sequence name, length (the records kept)
#   OUTPUT.dropped.tsv   locus tag, original sequence name, length, kept record (the records left out)
# It stops with an error, without writing OUTPUT, if a header has no [locus_tag=...].

set -euo pipefail

usage() {
    sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//' >&2
    exit 1
}

[ $# -ge 1 ] && [ $# -le 2 ] || usage
input=$1
if [ $# -eq 2 ]; then
    output=$2
else
    case "$input" in
        *.*) output="${input%.*}.locus_tag.${input##*.}" ;;
        *)   output="${input}.locus_tag" ;;
    esac
fi
mapping="${output}.mapping.tsv"
dropped="${output}.dropped.tsv"

[ -r "$input" ] || { echo "ERROR: cannot read input file [$input]" >&2; exit 1; }
for file in "$output" "$mapping" "$dropped"; do
    [ ! -e "$file" ] || { echo "ERROR: [$file] already exists, not overwriting it" >&2; exit 1; }
done

# Written under temporary names and moved to their final names only if everything went well
tmp_output=$(mktemp "${output}.tmp.XXXXXX")
tmp_mapping=$(mktemp "${mapping}.tmp.XXXXXX")
tmp_dropped=$(mktemp "${dropped}.tmp.XXXXXX")
trap 'rm -f "$tmp_output" "$tmp_mapping" "$tmp_dropped"' EXIT

# The input is read twice: the first pass measures every record and chooses the one kept per locus tag,
# the second pass writes the records kept.
awk -v mapping="$tmp_mapping" -v dropped="$tmp_dropped" '
function locus_tag(header) {
    # value of [locus_tag=...]
    if (match(header, /\[locus_tag=[^]]*\]/)) {
        return substr(header, RSTART + 11, RLENGTH - 12)
    }
    return ""
}

FNR == NR {                                   # first pass
    if (substr($0, 1, 1) == ">") {
        n++
        id[n] = substr($1, 2)
        tag[n] = locus_tag($0)
        if (tag[n] == "") {
            missing++
            print "ERROR: no [locus_tag=...] in header: " $0 > "/dev/stderr"
        }
    } else if (n > 0) {
        len[n] += length($0)
    }
    next
}

FNR == 1 {                                    # start of the second pass
    if (missing > 0) exit 1
    for (i = 1; i <= n; i++) {                # longest record per locus tag, the first one if equally long
        if (!(tag[i] in best) || len[i] > len[best[tag[i]]]) best[tag[i]] = i
    }
}

{                                             # second pass
    if (substr($0, 1, 1) == ">") {
        m++
        keep = (best[tag[m]] == m)
        if (keep) {
            print ">" tag[m]
            print tag[m] "\t" id[m] "\t" len[m] > mapping
        } else {
            print tag[m] "\t" id[m] "\t" len[m] "\t" id[best[tag[m]]] > dropped
        }
    } else if (keep) {
        print
    }
}
' "$input" "$input" > "$tmp_output"

# The names of the output must be unique
repeated=$(grep '^>' "$tmp_output" | sort | uniq -d | wc -l)
if [ "$repeated" -ne 0 ]; then
    echo "ERROR: $repeated sequence names are not unique in the output" >&2
    exit 1
fi

mv "$tmp_output" "$output"
mv "$tmp_mapping" "$mapping"
mv "$tmp_dropped" "$dropped"
trap - EXIT

echo "Input:   $input ($(grep -c '^>' "$input") sequences)"
echo "Output:  $output ($(grep -c '^>' "$output") sequences, one per locus tag)"
echo "Dropped: $(wc -l < "$dropped") sequences of locus tags with several records, listed in $dropped"
echo "Mapping: $mapping"
