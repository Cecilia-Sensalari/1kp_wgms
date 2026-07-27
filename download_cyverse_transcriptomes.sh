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
# The data type defaults to "both". OUTPUT_DIR defaults to the project directory
# used by the original script. Deduplication requires seqkit.

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
selection=${2:-both}
output_dir=${3:-source_data/2.transcriptomes/gymno_li_barker}
base_url=https://de.cyverse.org/anon-files/iplant/home/shared/commons_repo/curated/oneKP_capstone_2019/transcript_assemblies

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
# Other species and the retained seqkit reports are left untouched.
remove_intermediates() {
    local archive=$1
    local extract_dir=$2
    local merged_file=$3

    rm -f "$archive" "$merged_file"
    rm -rf "$extract_dir"
}

# Process one species and one data type. Keeping these operations in a function
# lets the main loop call exactly the same code for filtered and unfiltered data.
process_transcriptome() {
    local cyverse_name=$1
    local code=$2
    local data_type=$3
    local remote_name archive extract_dir fna_dir merged_file deduplicated_file
    local duplicate_output_dir processed_output_dir
    local duplicate_sequences_file duplicate_ids_file
    local -a fna_files

    # CyVerse uses different archive names for the two data types. Local files
    # are also separated into filtered and unfiltered directories.
    if [[ $data_type == filtered ]]; then
        remote_name="${code}.filtered.tar.bz2"
        archive="$output_dir/1.downloaded_compressed/filtered/${cyverse_name}.filtered.tar.bz2"
        extract_dir="$output_dir/2.merged/filtered/${cyverse_name}_filtered_FNA"
        fna_dir="$extract_dir/FILTERED/FNA"
        merged_file="$output_dir/2.merged/filtered/${cyverse_name}_filtered.duplic.FNA"
        duplicate_output_dir="$output_dir/3.rm_duplicates_by_id/filtered"
        duplicate_sequences_file="$duplicate_output_dir/${cyverse_name}_filtered_dup_seqs.fasta"
        duplicate_ids_file="$duplicate_output_dir/${cyverse_name}_filtered_dup_num_id.txt"
        processed_output_dir="$output_dir/4.processed_transcriptomes/filtered"
        deduplicated_file="$processed_output_dir/${cyverse_name}_filtered.FNA"
    else
        remote_name="${code}.fna.tar.bz2"
        archive="$output_dir/1.downloaded_compressed/unfiltered/${cyverse_name}.fna.tar.bz2"
        extract_dir="$output_dir/2.merged/unfiltered/${cyverse_name}_unfiltered_FNA"
        fna_dir="$extract_dir/FNA"
        merged_file="$output_dir/2.merged/unfiltered/${cyverse_name}_unfiltered.duplic.FNA"
        duplicate_output_dir="$output_dir/3.rm_duplicates_by_id/unfiltered"
        duplicate_sequences_file="$duplicate_output_dir/${cyverse_name}_unfiltered_dup_seqs.fasta"
        duplicate_ids_file="$duplicate_output_dir/${cyverse_name}_unfiltered_dup_num_id.txt"
        processed_output_dir="$output_dir/4.processed_transcriptomes/unfiltered"
        deduplicated_file="$processed_output_dir/${cyverse_name}_unfiltered.FNA"
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
        echo "  [$data_type] Removing download and merge intermediates, if still llexists"
        remove_intermediates "$archive" "$extract_dir" "$merged_file"
        return 0
    fi

    # Download only when the renamed local archive is not already available.
    if [[ -f $archive ]]; then
        echo "  [$data_type] Archive exists; skipping download."
    else
        echo "  [$data_type] Downloading $remote_name"
        wget -O "$archive" "$base_url/$cyverse_name/$remote_name"
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
        tar -xjf "$archive" -C "$extract_dir"
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
        return 1
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
    # species directory, and merged duplicate-containing input are no longer
    # needed. Keep the seqkit duplicate reports for later inspection.
    echo "  [$data_type] Removing download and merge intermediates"
    remove_intermediates "$archive" "$extract_dir" "$merged_file"
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

    echo "$onekp_cyverse_name"
    echo "  - code: $onekp_code"
    echo "  - species: $species_name"

    # This loop runs once for a single selection or twice when "both" was used.
    for data_type in "${data_types[@]}"; do
        process_transcriptome "$onekp_cyverse_name" "$onekp_code" "$data_type"
    done
done < "$species_list"
