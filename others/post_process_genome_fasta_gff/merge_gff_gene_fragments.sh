#!/usr/bin/env bash
# Merge the "gene" lines of a GFF that share the same locus_tag into a single gene line.
#
# Some annotations (e.g. NCBI "is_ordered=true" genes) split one gene into several gene lines with the same
# ID, Name and locus_tag but different coordinates, e.g. for CBR_g48825:
#   gene  BFEA01000716.1  662009..704099  -  ID=gene-CBR_g48825;Name=CBR_g48825;locus_tag=CBR_g48825
#   gene  BFEA01000716.1  656210..657741  -  ID=gene-CBR_g48825;Name=CBR_g48825;locus_tag=CBR_g48825
#   gene  BFEA01000716.1  654661..656150  -  ID=gene-CBR_g48825;Name=CBR_g48825;locus_tag=CBR_g48825
# ksrates' collinearity analysis makes one gene-list entry per GFF line, so such a gene would be listed once per
# fragment. This script replaces the fragments by ONE line, from the lowest start to the highest end of the
# fragments, placed where the first fragment was. It keeps the columns and the attributes of that first fragment
# and appends ";merged_fragments=N" to its attributes. All other lines are copied unchanged.
#
# Fragments are merged only if they are on the same contig and strand; otherwise they are left as they are and
# reported. Lines without a locus_tag attribute are never touched.
#
# For each merged gene the span is compared with the span of the mRNA line(s) having the same locus_tag, to
# spot merges that look wrong.
#
# Usage:   merge_gff_gene_fragments.sh INPUT.gff [OUTPUT.gff]
# Default OUTPUT: INPUT with ".merged_genes" before the extension, in the same directory.
# The input file is never changed, and an existing output file is never overwritten. Besides OUTPUT it writes
#   OUTPUT.merge_log.tsv   one row per locus tag that has several gene lines
# Columns of the log: locus_tag, contig, strand, n_fragments, fragments, merged_start, merged_end,
#   merged_length, mrna_span, status. status is one of
#   MERGED_MATCHES_MRNA      the merged span equals the mRNA span
#   MERGED_DIFFERS_FROM_MRNA the merged span is different from the mRNA span: CHECK these genes
#   MERGED_NO_MRNA           merged, but there is no mRNA line to compare with
#   NOT_MERGED_CONTIG_OR_STRAND  fragments on different contigs or strands: left unchanged

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
        *.*) output="${input%.*}.merged_genes.${input##*.}" ;;
        *)   output="${input}.merged_genes" ;;
    esac
fi
log="${output}.merge_log.tsv"

[ -r "$input" ] || { echo "ERROR: cannot read input file [$input]" >&2; exit 1; }
for file in "$output" "$log"; do
    [ ! -e "$file" ] || { echo "ERROR: [$file] already exists, not overwriting it" >&2; exit 1; }
done

# Written under temporary names and moved to their final names only if everything went well
tmp_output=$(mktemp "${output}.tmp.XXXXXX")
tmp_log=$(mktemp "${log}.tmp.XXXXXX")
trap 'rm -f "$tmp_output" "$tmp_log" "${tmp_log}.unsorted"' EXIT
: > "${tmp_log}.unsorted"   # exists even when no gene is split (the awk below then writes nothing to it)

# The input is read twice: the first pass collects the gene fragments and the mRNA spans, the second pass
# writes the output.
awk -F'\t' -v OFS='\t' -v logfile="${tmp_log}.unsorted" '
function attr(text, key,    n, parts, i) {
    # value of "key=value" in the attributes column, "" if absent
    n = split(text, parts, ";")
    for (i = 1; i <= n; i++) {
        if (index(parts[i], key "=") == 1) return substr(parts[i], length(key) + 2)
    }
    return ""
}
function usable(    ) { return substr($0, 1, 1) != "#" && NF == 9 }

FNR == NR {                                   # first pass
    if (!usable()) next
    if ($3 == "gene") {
        tag = attr($9, "locus_tag")
        if (tag == "") next
        n[tag]++
        if (n[tag] == 1) {
            contig[tag] = $1; strand[tag] = $7
            lo[tag] = $4 + 0; hi[tag] = $5 + 0
            fragments[tag] = $4 ".." $5
        } else {
            if ($1 != contig[tag] || $7 != strand[tag]) different[tag] = 1
            if ($4 + 0 < lo[tag]) lo[tag] = $4 + 0
            if ($5 + 0 > hi[tag]) hi[tag] = $5 + 0
            fragments[tag] = fragments[tag] "," $4 ".." $5
        }
    } else if ($3 == "mRNA") {
        tag = attr($9, "locus_tag")
        if (tag == "") next
        if (!(tag in mlo) || $4 + 0 < mlo[tag]) mlo[tag] = $4 + 0
        if (!(tag in mhi) || $5 + 0 > mhi[tag]) mhi[tag] = $5 + 0
    }
    next
}

{                                             # second pass
    if (usable() && $3 == "gene") {
        tag = attr($9, "locus_tag")
        if (tag != "" && n[tag] > 1 && !(tag in different)) {
            if (!(tag in written)) {
                written[tag] = 1
                $4 = lo[tag]; $5 = hi[tag]
                sub(/;$/, "", $9)
                $9 = $9 ";merged_fragments=" n[tag]
                print
            }
            next
        }
    }
    print
}

END {
    for (tag in n) {
        if (n[tag] < 2) continue
        mrna = (tag in mlo) ? mlo[tag] ".." mhi[tag] : "-"
        if (tag in different) status = "NOT_MERGED_CONTIG_OR_STRAND"
        else if (!(tag in mlo)) status = "MERGED_NO_MRNA"
        else if (mlo[tag] == lo[tag] && mhi[tag] == hi[tag]) status = "MERGED_MATCHES_MRNA"
        else status = "MERGED_DIFFERS_FROM_MRNA"
        print tag, contig[tag], strand[tag], n[tag], fragments[tag], lo[tag], hi[tag], hi[tag] - lo[tag] + 1, mrna, status > logfile
    }
}
' "$input" "$input" > "$tmp_output"

{
    printf 'locus_tag\tcontig\tstrand\tn_fragments\tfragments\tmerged_start\tmerged_end\tmerged_length\tmrna_span\tstatus\n'
    sort -t$'\t' -k1,1 "${tmp_log}.unsorted"
} > "$tmp_log"
rm -f "${tmp_log}.unsorted"

mv "$tmp_output" "$output"
mv "$tmp_log" "$log"
trap - EXIT

genes_before=$(awk -F'\t' '!/^#/ && NF == 9 && $3 == "gene"' "$input" | wc -l)
genes_after=$(awk -F'\t' '!/^#/ && NF == 9 && $3 == "gene"' "$output" | wc -l)
echo "Input:   $input ($genes_before gene lines)"
echo "Output:  $output ($genes_after gene lines)"
echo "Log:     $log"
echo "Locus tags with several gene lines, by status:"
tail -n +2 "$log" | cut -f10 | sort | uniq -c | sed 's/^/  /'
