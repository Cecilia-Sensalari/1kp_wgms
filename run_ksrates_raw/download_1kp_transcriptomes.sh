# Download, extract, merge, and deduplicate 1KP transcriptomes.
#
# [Generated via AI, tested by Cecilia]
#
# Usage:
#   bash download_cyverse_transcriptomes.sh SPECIES_LIST [both|unfiltered|filtered] [OUTPUT_DIR] [gdrive|cyverse]
#
# SPECIES_LIST contains one 1KP sample ID followed by a name per line. For example:
#   ACSA-Species_name_rest
#
# The data type defaults to "unfiltered". OUTPUT_DIR defaults to the 1KP
# transcriptome source-data directory. SOURCE defaults to "gdrive".
# Deduplication requires seqkit. Google Drive downloads require the Python
# package gdown.

# Stop when a command fails, an undefined variable is used, or a pipeline fails.
set -euo pipefail

# Check that the required species list and no more than three optional arguments
# were supplied.
if (( $# < 1 || $# > 4 )); then
    echo "Usage: $0 SPECIES_LIST [both|unfiltered|filtered] [OUTPUT_DIR] [gdrive|cyverse]" >&2
    exit 1
fi

# Command-line settings. The :- syntax supplies a default when an optional
# argument was not provided.
species_list=$1
selection=${2:-unfiltered}
output_dir=${3:-/group/esb/cesen/1kp/source_data/2.transcriptomes/1kp}
download_source=${4:-gdrive}
cyverse_base_url=https://de.cyverse.org/anon-files/iplant/home/shared/commons_repo/curated/oneKP_capstone_2019/transcript_assemblies
cyverse_resolver_url=https://datacommons.cyverse.org/api/list/iplant/home/shared/commons_repo/curated/oneKP_capstone_2019/transcript_assemblies
gdrive_folder_id=18AOvneP_1l5uzE7tVWPKVR9MAkhrtA2N
gdrive_folder_url="https://drive.google.com/drive/folders/$gdrive_folder_id"
remote_index_file=$(mktemp "${TMPDIR:-/tmp}/onekp_cyverse_directory_index.XXXXXX")
keep_dup_reports=${KEEP_DUP_REPORTS:-false}
trap 'rm -f "$remote_index_file"' EXIT

if [[ $output_dir == gdrive || $output_dir == cyverse ]]; then
    download_source=$output_dir
    output_dir=/group/esb/cesen/1kp/source_data/2.transcriptomes/1kp
fi

missing_log_file="$output_dir/missing_transcriptomes.tsv"

case $download_source in
    gdrive|cyverse) ;;
    *)
        echo "Invalid download source '$download_source'; use gdrive or cyverse." >&2
        exit 1
        ;;
esac

if [[ ! -f $species_list ]]; then
    echo "Species list does not exist: $species_list" >&2
    exit 1
fi

# Convert the requested selection into an array. The same processing loop can
# then handle one data type or both data types without duplicated code.
case $selection in
    both)       data_types=(unfiltered filtered) ;;
    unfiltered) data_types=(unfiltered) ;;
    filtered)   data_types=(filtered) ;;
    *)
        echo "Invalid data type '$selection'; use both, unfiltered, or filtered." >&2
        exit 1
        ;;
esac

# Load seqkit on the cluster only if it is not already on PATH. Keeping this
# check here allows the download/extraction stages to remain usable on systems
# where seqkit is installed normally rather than through environment modules.
ensure_seqkit() {
    if command -v seqkit >/dev/null 2>&1; then
        return 0
    fi

    if command -v module >/dev/null 2>&1; then
        echo "  Loading seqkit module"
        if ! module load seqkit/x86_64/0.7.1; then
            echo "Could not load the seqkit module (seqkit/x86_64/0.7.1)." >&2
            return 1
        fi
    fi

    if ! command -v seqkit >/dev/null 2>&1; then
        echo "seqkit is required for deduplication but was not found on PATH." >&2
        return 1
    fi
}

# Remove only the intermediate files belonging to one species and data type.
remove_intermediates() {
    local archive=$1
    local extract_dir=$2
    local merged_file=$3

    rm -f "$archive" "$merged_file"
    rm -rf "$extract_dir"
}

remove_duplicate_reports() {
    local duplicate_sequences_file=$1
    local duplicate_ids_file=$2

    if [[ $keep_dup_reports == true ]]; then
        return 0
    fi

    rm -f "$duplicate_sequences_file" "$duplicate_ids_file"
}

log_missing_transcriptome() {
    local local_name=$1
    local remote_name_dir=$2
    local code=$3
    local data_type=$4
    local reason=$5
    local attempted_url=$6

    mkdir -p "$(dirname "$missing_log_file")"
    if [[ ! -f $missing_log_file ]]; then
        printf "timestamp\tcode\tlocal_name\tremote_directory\tdata_type\treason\tattempted_url\n" > "$missing_log_file"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$(date -Is)" "$code" "$local_name" "$remote_name_dir" "$data_type" "$reason" "$attempted_url" \
        >> "$missing_log_file"
}

# Check that a downloaded tar.bz2 archive is present, non-empty, and readable.
# This prevents a previous partial download from being reused on reruns.
archive_is_valid() {
    local archive=$1

    [[ -s $archive ]] || return 1
    tar -tjf "$archive" >/dev/null 2>&1
}

# Download through a temporary file and move it into place only after the
# download has completed and the archive passes a tar/bzip2 integrity check.
download_archive() {
    local url=$1
    local archive=$2
    local tmp_archive="${archive}.tmp"

    rm -f "$tmp_archive"
    if ! wget -O "$tmp_archive" "$url"; then
        rm -f "$tmp_archive"
        echo "Archive download failed: $url" >&2
        return 1
    fi

    if archive_is_valid "$tmp_archive"; then
        mv "$tmp_archive" "$archive"
    else
        rm -f "$tmp_archive"
        echo "Downloaded archive failed integrity check: $url" >&2
        return 1
    fi
}

download_gdrive_archive() {
    local remote_name_dir=$1
    local remote_name=$2
    local archive=$3
    local tmp_archive="${archive}.tmp"
    local remote_dir_id

    remote_dir_id=$(awk -F '\t' -v dir="$remote_name_dir" '$1 == dir {print $2; exit}' "$remote_index_file")
    if [[ -z $remote_dir_id ]]; then
        echo "Google Drive directory ID not found for $remote_name_dir" >&2
        return 1
    fi

    rm -f "$tmp_archive"
    if ! python3 - "$remote_name_dir" "$remote_name" "$remote_dir_id" "$tmp_archive" <<'PY'; then
from html.parser import HTMLParser
import html
import re
import sys
import urllib.parse

import requests

remote_name_dir = sys.argv[1]
remote_name = sys.argv[2]
remote_dir_id = sys.argv[3]
tmp_archive = sys.argv[4]

try:
    import gdown
except Exception as exc:
    raise SystemExit(
        "Python package gdown is required for Google Drive downloads. "
        "Install it or rerun with source 'cyverse'. "
        f"Import error: {exc}"
    )

class LinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self._href = None
        self._text = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            self._href = dict(attrs).get("href", "")
            self._text = []

    def handle_data(self, data):
        if self._href is not None:
            self._text.append(data)

    def handle_endtag(self, tag):
        if tag == "a" and self._href is not None:
            name = html.unescape("".join(self._text)).strip()
            self.links.append((self._href, name))
            self._href = None
            self._text = []

def embedded_folder_links(folder_id):
    query = urllib.parse.urlencode({"id": folder_id})
    response = requests.get(
        f"https://drive.google.com/embeddedfolderview?{query}",
        headers={"User-Agent": "Mozilla/5.0"},
        timeout=60,
    )
    response.raise_for_status()
    text = response.text
    parser = LinkParser()
    parser.feed(text)
    return parser.links

matches = []
for href, name in embedded_folder_links(remote_dir_id):
    file_match = re.match(r"https://drive\.google\.com/file/d/([-\w]{25,})/view", href)
    if file_match and name == remote_name:
        matches.append(file_match.group(1))

if len(matches) != 1:
    raise SystemExit(f"Expected one Google Drive file named {remote_name} in {remote_name_dir}, found {len(matches)}")

gdown.download(id=matches[0], output=tmp_archive, quiet=False, use_cookies=False)
PY
        rm -f "$tmp_archive"
        echo "Google Drive archive download failed: $remote_name_dir/$remote_name" >&2
        return 1
    fi

    if archive_is_valid "$tmp_archive"; then
        mv "$tmp_archive" "$archive"
    else
        rm -f "$tmp_archive"
        echo "Downloaded Google Drive archive failed integrity check: $remote_name_dir/$remote_name" >&2
        return 1
    fi
}

download_remote_archive() {
    local archive_locator=$1
    local archive=$2
    local remote_name_dir=$3
    local remote_name=$4

    case $download_source in
        gdrive) download_gdrive_archive "$remote_name_dir" "$remote_name" "$archive" ;;
        cyverse) download_archive "$archive_locator" "$archive" ;;
    esac
}

# Build a temporary list of public OneKP directory names from the selected
# source. The remote directory names are not always exactly CODE + species name,
# for example:
#   QFND-Cyanophora_paradoxa-CCAC_0074
# This lets us resolve the exact remote directory from the 1KP sample ID while
# keeping predictable local output filenames based on the species-list entry.
ensure_cyverse_remote_index() {
    if [[ -s $remote_index_file ]]; then
        return 0
    fi

    echo "Resolving CyVerse directory names from $cyverse_resolver_url" >&2
    rm -f "$remote_index_file"

    if ! python3 - "$cyverse_resolver_url" "$remote_index_file" <<'PY'; then
import json
import re
import sys
import urllib.parse
import urllib.request

resolver_url = sys.argv[1]
remote_index_file = sys.argv[2]
labels = []
page = 0

while True:
    query = urllib.parse.urlencode({"page": page, "sort-col": "NAME", "sort-dir": "asc"})
    with urllib.request.urlopen(f"{resolver_url}?{query}") as response:
        payload = json.load(response)

    folders = payload.get("folders", [])
    labels.extend(
        item["label"]
        for item in folders
        if re.match(r"^[A-Za-z0-9_.]+-", item.get("label", ""))
    )

    total = int(payload.get("total", 0))
    page += 1
    if not folders or page * 100 >= total:
        break

with open(remote_index_file, "w", encoding="utf-8") as handle:
    handle.write("\n".join(f"{label}\t" for label in labels))
    handle.write("\n")
PY
        rm -f "$remote_index_file"
        echo "Could not download or parse OneKP directory index from $cyverse_resolver_url" >&2
        return 1
    fi

    if [[ ! -s $remote_index_file ]]; then
        rm -f "$remote_index_file"
        echo "OneKP directory index from $cyverse_resolver_url was empty" >&2
        return 1
    fi
}

ensure_gdrive_remote_index() {
    if [[ -s $remote_index_file ]]; then
        return 0
    fi

    echo "Resolving Google Drive directory names from $gdrive_folder_url" >&2
    rm -f "$remote_index_file"

    if ! python3 - "$gdrive_folder_id" "$remote_index_file" <<'PY'; then
from html.parser import HTMLParser
import html
import re
import sys
import urllib.parse

import requests

folder_id = sys.argv[1]
remote_index_file = sys.argv[2]

class LinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self._href = None
        self._text = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            self._href = dict(attrs).get("href", "")
            self._text = []

    def handle_data(self, data):
        if self._href is not None:
            self._text.append(data)

    def handle_endtag(self, tag):
        if tag == "a" and self._href is not None:
            name = html.unescape("".join(self._text)).strip()
            self.links.append((self._href, name))
            self._href = None
            self._text = []

query = urllib.parse.urlencode({"id": folder_id})
response = requests.get(
    f"https://drive.google.com/embeddedfolderview?{query}",
    headers={"User-Agent": "Mozilla/5.0"},
    timeout=60,
)
response.raise_for_status()
text = response.text

parser = LinkParser()
parser.feed(text)

rows = []
for href, name in parser.links:
    folder_match = re.match(r"https://drive\.google\.com/drive/folders/([-\w]{25,})", href)
    if folder_match and name and re.match(r"^[A-Za-z0-9_.]+-", name):
        rows.append((name, folder_match.group(1)))

with open(remote_index_file, "w", encoding="utf-8") as handle:
    for directory_name, directory_id in rows:
        handle.write(f"{directory_name}\t{directory_id}\n")
PY
        rm -f "$remote_index_file"
        echo "Could not download or parse Google Drive directory index from $gdrive_folder_url" >&2
        return 1
    fi

    if [[ ! -s $remote_index_file ]]; then
        rm -f "$remote_index_file"
        echo "Google Drive directory index from $gdrive_folder_url was empty" >&2
        return 1
    fi
}

ensure_remote_index() {
    case $download_source in
        gdrive) ensure_gdrive_remote_index ;;
        cyverse) ensure_cyverse_remote_index ;;
    esac
}

resolve_remote_name() {
    local requested_name=$1
    local code=$2
    local -a matches

    if ! ensure_remote_index; then
        echo "  Could not resolve remote directory from [$code-]; skipping [$requested_name]" >&2
        return 1
    fi

    mapfile -t matches < <(awk -F '\t' -v code="$code" 'index($1, code "-") == 1 && !seen[$1]++ {print $1}' "$remote_index_file")

    if (( ${#matches[@]} == 0 )); then
        echo "  No remote directory starts with [$code-]; skipping [$requested_name]" >&2
        return 1
    fi

    if (( ${#matches[@]} == 1 )); then
        echo "${matches[0]}"
        return 0
    fi

    echo "  Multiple remote directories match [$code]; skipping [$requested_name]" >&2
    printf '    %s\n' "${matches[@]}" >&2
    return 1
}

# Process one species and one data type. Keeping these operations in a function
# lets the main loop call exactly the same code for filtered and unfiltered data.
process_transcriptome() {
    local local_name=$1
    local remote_name_dir=$2
    local code=$3
    local data_type=$4
    local remote_name archive extract_dir fna_dir merged_file deduplicated_file
    local duplicate_output_dir processed_output_dir
    local duplicate_sequences_file duplicate_ids_file
    local archive_locator
    local -a fna_files

    # CyVerse uses different archive names for the two data types. Local files
    # are also separated into filtered and unfiltered directories.
    if [[ $data_type == filtered ]]; then
        remote_name="${code}.filtered.tar.bz2"
        archive="$output_dir/1.downloaded_compressed/filtered/${local_name}.filtered.tar.bz2"
        extract_dir="$output_dir/2.merged/filtered/${local_name}_filtered_FNA"
        fna_dir="$extract_dir/FILTERED/FNA"
        merged_file="$output_dir/2.merged/filtered/${local_name}_filtered.duplic.FNA"
        duplicate_output_dir="$output_dir/3.rm_duplicates_by_id/filtered"
        duplicate_sequences_file="$duplicate_output_dir/${local_name}_filtered_dup_seqs.fasta"
        duplicate_ids_file="$duplicate_output_dir/${local_name}_filtered_dup_num_id.txt"
        processed_output_dir="$output_dir/4.processed_transcriptomes/filtered"
        deduplicated_file="$processed_output_dir/${local_name}_filtered.FNA"
    else
        remote_name="${code}.fna.tar.bz2"
        archive="$output_dir/1.downloaded_compressed/unfiltered/${local_name}.fna.tar.bz2"
        extract_dir="$output_dir/2.merged/unfiltered/${local_name}_unfiltered_FNA"
        fna_dir="$extract_dir/FNA"
        merged_file="$output_dir/2.merged/unfiltered/${local_name}_unfiltered.duplic.FNA"
        duplicate_output_dir="$output_dir/3.rm_duplicates_by_id/unfiltered"
        duplicate_sequences_file="$duplicate_output_dir/${local_name}_unfiltered_dup_seqs.fasta"
        duplicate_ids_file="$duplicate_output_dir/${local_name}_unfiltered_dup_num_id.txt"
        processed_output_dir="$output_dir/4.processed_transcriptomes/unfiltered"
        deduplicated_file="$processed_output_dir/${local_name}_unfiltered.FNA"
    fi

    # mkdir -p is harmless when these parent directories already exist.
    mkdir -p \
        "$(dirname "$archive")" \
        "$(dirname "$merged_file")" \
        "$duplicate_output_dir" \
        "$processed_output_dir"

    # A completed deduplication is the final result. On a rerun, clean up any
    # stale intermediates and skip downloading and processing this species.
    if [[ -f $deduplicated_file ]]; then
        echo "  [$data_type] Deduplicated file exists; skipping processing."
        echo "  [$data_type] Removing reproducible intermediates, if still present"
        remove_intermediates "$archive" "$extract_dir" "$merged_file"
        remove_duplicate_reports "$duplicate_sequences_file" "$duplicate_ids_file"
        return 0
    fi

    # Download only when the renamed local archive is not already available.
    # If a previous partial/corrupt archive is found, remove and redownload it.
    if [[ $download_source == gdrive ]]; then
        archive_locator="gdrive://$remote_name_dir/$remote_name"
    else
        archive_locator="$cyverse_base_url/$remote_name_dir/$remote_name"
    fi

    downloaded_archive=false
    if [[ -f $archive ]]; then
        if archive_is_valid "$archive"; then
            echo "  [$data_type] Archive exists and passed integrity check; skipping download."
        else
            echo "  [$data_type] Existing archive is incomplete or corrupt; redownloading."
            rm -f "$archive"
            rm -rf "$extract_dir"
            rm -f "$merged_file"
            if ! download_remote_archive "$archive_locator" "$archive" "$remote_name_dir" "$remote_name"; then
                log_missing_transcriptome "$local_name" "$remote_name_dir" "$code" "$data_type" "archive_download_or_integrity_failed" "$archive_locator"
                echo "  [$data_type] Logged missing/invalid archive and skipped." >&2
                return 0
            fi
            downloaded_archive=true
        fi
    else
        echo "  [$data_type] Downloading $remote_name"
        if ! download_remote_archive "$archive_locator" "$archive" "$remote_name_dir" "$remote_name"; then
            log_missing_transcriptome "$local_name" "$remote_name_dir" "$code" "$data_type" "archive_download_or_integrity_failed" "$archive_locator"
            echo "  [$data_type] Logged missing/invalid archive and skipped." >&2
            return 0
        fi
        downloaded_archive=true
    fi

    if [[ $downloaded_archive == true ]]; then
        rm -rf "$extract_dir"
        rm -f "$merged_file"
    fi

    # nullglob makes an unmatched *.FNA pattern expand to nothing instead of
    # being stored as the literal text "*.FNA".
    shopt -s nullglob
    fna_files=("$fna_dir"/*.FNA)
    shopt -u nullglob

    # Existing FNA files indicate that this archive was already extracted.
    if (( ${#fna_files[@]} > 0 )); then
        echo "  [$data_type] Extracted FNA files exist; skipping extraction."
    else
        echo "  [$data_type] Extracting archive"
        mkdir -p "$extract_dir"
        if ! tar -xjf "$archive" -C "$extract_dir"; then
            log_missing_transcriptome "$local_name" "$remote_name_dir" "$code" "$data_type" "archive_extraction_failed" "$archive_locator"
            echo "  [$data_type] Logged extraction failure and skipped." >&2
            rm -rf "$extract_dir"
            rm -f "$merged_file"
            return 0
        fi
        # Refresh the file list after extraction so it can be used for merging.
        shopt -s nullglob
        fna_files=("$fna_dir"/*.FNA)
        shopt -u nullglob
    fi

    # Do not overwrite an existing merged result. If extraction produced no
    # FNA files, report the problem instead of creating an empty result.
    if [[ -f $merged_file ]]; then
        echo "  [$data_type] Merged file exists; skipping merge."
    elif (( ${#fna_files[@]} == 0 )); then
        echo "  [$data_type] No extracted .FNA files found in $fna_dir" >&2
        log_missing_transcriptome "$local_name" "$remote_name_dir" "$code" "$data_type" "no_extracted_fna_files" "$archive_locator"
        echo "  [$data_type] Logged missing FNA files and skipped." >&2
        return 0
    else
        echo "  [$data_type] Merging ${#fna_files[@]} FNA file(s)"
        cat "${fna_files[@]}" > "$merged_file"
    fi

    # Remove duplicate sequences with seqkit. Skip this step when its primary
    # output already exists, preserving the same resumable behavior as the
    # download, extraction, and merge stages.
    if [[ -f $deduplicated_file ]]; then
        echo "  [$data_type] Deduplicated file exists; skipping deduplication."
    else
        ensure_seqkit
        echo "  [$data_type] Removing duplicate sequences"
        seqkit rmdup \
            -d "$duplicate_sequences_file" \
            -D "$duplicate_ids_file" \
            "$merged_file" > "$deduplicated_file"
        echo "  [$data_type] Retained $(awk '/^>/{count++} END{print count + 0}' "$deduplicated_file") sequence(s)"
    fi

    # Once the final FNA has been created, the downloaded archive, extracted
    # species directory, merged duplicate-containing input, and duplicate
    # reports are reproducible and are removed by default.
    echo "  [$data_type] Removing reproducible intermediates"
    remove_intermediates "$archive" "$extract_dir" "$merged_file"
    remove_duplicate_reports "$duplicate_sequences_file" "$duplicate_ids_file"
}

# Read every entry in the species list. The condition after || also processes a
# final line when the input file does not end with a newline character.
while IFS= read -r onekp_cyverse_name || [[ -n $onekp_cyverse_name ]]; do
    # Accept Windows line endings and ignore empty lines and comments.
    onekp_cyverse_name=${onekp_cyverse_name%$'\r'}
    [[ -z $onekp_cyverse_name || $onekp_cyverse_name == \#* ]] && continue

    if [[ $onekp_cyverse_name == *-* ]]; then
        # Split "CODE-Species_name" at the first hyphen. Only CODE is used for
        # remote CyVerse lookup; the full input line is kept for local filenames.
        onekp_code=${onekp_cyverse_name%%-*}
        species_name=${onekp_cyverse_name#*-}
    else
        onekp_code=$onekp_cyverse_name
        species_name=
    fi

    echo "$onekp_cyverse_name"
    echo "  - source: $download_source"
    echo "  - code: $onekp_code"
    if [[ -n $species_name ]]; then
        echo "  - species: $species_name"
    fi

    if ! resolved_cyverse_name=$(resolve_remote_name "$onekp_cyverse_name" "$onekp_code"); then
        echo "  - remote directory: unresolved"
        for data_type in "${data_types[@]}"; do
            log_missing_transcriptome "$onekp_cyverse_name" "" "$onekp_code" "$data_type" "no_unique_remote_directory" "$download_source"
        done
        continue
    fi

    echo "  - remote directory: $resolved_cyverse_name"

    # This loop runs once for a single selection or twice when "both" was used.
    for data_type in "${data_types[@]}"; do
        process_transcriptome "$onekp_cyverse_name" "$resolved_cyverse_name" "$onekp_code" "$data_type"
    done
done < "$species_list"
