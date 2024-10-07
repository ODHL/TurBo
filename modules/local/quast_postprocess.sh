#!/bin/bash

# Check if input argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <report_file>"
    exit 1
fi

TSV_FILE="$1"

# Check if the TSV file exists
if [ ! -f "$TSV_FILE" ]; then
    echo "Error: File '$TSV_FILE' not found!"
    exit 1
fi

# Read the TSV file and extract required values
while IFS=$'\t' read -r line; do
    case "$line" in
        *"Total length"*)
            echo "${line##*$'\t'}" > "GENOME_LENGTH"
            ;;
        *"# contigs"*)
            echo "${line##*$'\t'}" > "NUMBER_CONTIGS"
            ;;
        *"N50"*)
            echo "${line##*$'\t'}" > "N50_VALUE"
            ;;
        *"GC"*)
            echo "${line##*$'\t'}" > "GC_PERCENT"
            ;;
        *"Largest contig"*)
            echo "${line##*$'\t'}" > "LARGEST_CONTIG"
            ;;
        *"# N's per 100 kbp"*)
            echo "${line##*$'\t'}" > "UNCALLED_BASES"
            ;;
    esac
done < "$TSV_FILE"