#!/usr/bin/env bash

# Download, extract, merge, and deduplicate 1KP transcriptomes from CyVerse.
#
# [Generated via AI, tested by Cecilia]
#
# Usage:
#   bash download_cyverse_transcriptomes.sh SPECIES_LIST [both|unfiltered|filtered] [OUTPUT_DIR]
#
# SPECIES_LIST contains one CyVerse directory name per line, for example:
#   ACSA-Species_name_rest
#
# The data type defaults to "unfiltered". OUTPUT_DIR defaults to the 1KP
# transcriptome source-data directory. Deduplication requires seqkit.

# Stop when a command fails, an undefined variable is used, or a pipeline fails.
set -euo pipefail

# Check that the required species list and no more than two optional arguments
# were supplied.
if (( $# < 1 || $# > 3 )); then
    echo "Usage: $0 SPECIES_LIST [both|unfiltered|filtered] [OUTPUT_DIR]" >&2
    exit 1
fi

# Command-line settings. The :- syntax supplies a default when an optional
# argument was not provided.
species_list=$1
selection=${2:-unfiltered}
output_dir=${3:-/group/esb/cesen/1kp/source_data/2.transcriptomes/1kp}
base_url=https://de.cyverse.org/anon-files/iplant/home/shared/commons_repo/curated/oneKP_capstone_2019/transcript_assemblies
resolver_url=https://web.corral.tacc.utexas.edu/OneKP/
remote_index_file="/group/esb/cesen/1kp/code/1kp_wgms/run_ksrates_raw/onekp_cyverse_directory_index.txt"
missing_log_file="$output_dir/missing_transcriptomes.tsv"
keep_dup_reports=${KEEP_DUP_REPORTS:-false}

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

# Build a local cache of public OneKP directory names. The CyVerse directory
# names are not always exactly CODE + species name, for example:
#   QFND-Cyanophora_paradoxa-CCAC_0074
# This index lets us resolve a CODE or CODE-species hint to the exact remote
# directory while keeping predictable local output filenames.
ensure_remote_index() {
    local tmp_index="${remote_index_file}.tmp"

    if [[ -s $remote_index_file ]] && ! grep -q '"' "$remote_index_file"; then
        return 0
    fi

    echo "Resolving CyVerse directory names from $resolver_url" >&2
    mkdir -p "$(dirname "$remote_index_file")"
    rm -f "$tmp_index"

    if ! wget -q -O - "$resolver_url" \
        | grep -o 'href="[^"]*/"' \
        | sed 's/^href="//; s#/"$##' \
        | grep -E '^[A-Za-z0-9_.]+-' > "$tmp_index"; then
        rm -f "$tmp_index"
        echo "Could not download or parse OneKP directory index from $resolver_url" >&2
        return 1
    fi

    if [[ ! -s $tmp_index ]]; then
        rm -f "$tmp_index"
        echo "OneKP directory index from $resolver_url was empty" >&2
        return 1
    fi

    mv "$tmp_index" "$remote_index_file"
}

resolve_remote_cyverse_name() {
    local requested_name=$1
    local code=$2
    local -a matches prefix_matches

    if ! ensure_remote_index; then
        echo "  Could not resolve exact remote directory; trying requested name [$requested_name]" >&2
        echo "$requested_name"
        return 0
    fi

    mapfile -t matches < <(awk -v code="$code" 'index($0, code "-") == 1 {print}' "$remote_index_file")

    if (( ${#matches[@]} == 0 )); then
        echo "  No remote directory starts with [$code-]; trying requested name [$requested_name]" >&2
        echo "$requested_name"
        return 0
    fi

    for match in "${matches[@]}"; do
        if [[ $match == "$requested_name" ]]; then
            echo "$match"
            return 0
        fi
    done

    mapfile -t prefix_matches < <(printf '%s\n' "${matches[@]}" | awk -v requested="$requested_name" '$0 == requested || index($0, requested "-") == 1 {print}')
    if (( ${#prefix_matches[@]} == 1 )); then
        echo "${prefix_matches[0]}"
        return 0
    fi

    if (( ${#matches[@]} == 1 )); then
        echo "${matches[0]}"
        return 0
    fi

    echo "  Multiple remote directories match [$code]; choosing first: ${matches[0]}" >&2
    printf '    %s\n' "${matches[@]}" >&2
    echo "${matches[0]}"
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
    local archive_url
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
    archive_url="$base_url/$remote_name_dir/$remote_name"
    downloaded_archive=false
    if [[ -f $archive ]]; then
        if archive_is_valid "$archive"; then
            echo "  [$data_type] Archive exists and passed integrity check; skipping download."
        else
            echo "  [$data_type] Existing archive is incomplete or corrupt; redownloading."
            rm -f "$archive"
            rm -rf "$extract_dir"
            rm -f "$merged_file"
            if ! download_archive "$archive_url" "$archive"; then
                log_missing_transcriptome "$local_name" "$remote_name_dir" "$code" "$data_type" "archive_download_or_integrity_failed" "$archive_url"
                echo "  [$data_type] Logged missing/invalid archive and skipped." >&2
                return 0
            fi
            downloaded_archive=true
        fi
    else
        echo "  [$data_type] Downloading $remote_name"
        if ! download_archive "$archive_url" "$archive"; then
            log_missing_transcriptome "$local_name" "$remote_name_dir" "$code" "$data_type" "archive_download_or_integrity_failed" "$archive_url"
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
            log_missing_transcriptome "$local_name" "$remote_name_dir" "$code" "$data_type" "archive_extraction_failed" "$archive_url"
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
        log_missing_transcriptome "$local_name" "$remote_name_dir" "$code" "$data_type" "no_extracted_fna_files" "$archive_url"
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

    if [[ $onekp_cyverse_name != *-* ]]; then
        echo "Invalid species entry (expected CODE-Species_name): $onekp_cyverse_name" >&2
        exit 1
    fi

    # Split "CODE-Species_name" at the first hyphen.
    onekp_code=${onekp_cyverse_name%%-*}
    species_name=${onekp_cyverse_name#*-}
    resolved_cyverse_name=$(resolve_remote_cyverse_name "$onekp_cyverse_name" "$onekp_code")

    echo "$onekp_cyverse_name"
    echo "  - code: $onekp_code"
    echo "  - species: $species_name"
    echo "  - remote directory: $resolved_cyverse_name"

    # This loop runs once for a single selection or twice when "both" was used.
    for data_type in "${data_types[@]}"; do
        process_transcriptome "$onekp_cyverse_name" "$resolved_cyverse_name" "$onekp_code" "$data_type"
    done
done < "$species_list"
