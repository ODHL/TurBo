#!/bin/bash

input=$1
###########################################################################################################################
## COVERAGE
###########################################################################################################################
# Output file
for type in raw cleaned; do
    # set output
    output_file="$input/cov_${type}_output.tsv"

    # Create the output header
    echo -e "filename\trname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq" > "$output_file"

    # Loop through all the files
    for file in $input/flagstat/*$type*cov*; do
        # Check if the file is not empty and is a valid file
        if [[ -f "$file" && -s "$file" ]]; then
            # Process the file and append to output
            awk -v filename="$file" 'NR > 1 {printf "%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.5f\t%d\t%.1f\t\n", filename, $1, $2, $3, $4, $5, $6, $7, $8, $9}' "$file" >> "$output_file"
        fi
    done

    echo "Combined output written to $output_file"
done

###########################################################################################################################
## stats
###########################################################################################################################
# Output file
output_file="$input/clean_output.tsv"

# Create or clear the output file
echo -e "filename\tgenome\tReads_processed\tUnmapped\tOne_hit_one_genome" > "tmp.txt"

# Process each FASTQ screen file
for file in $input/fastq_screen/*.txt; do  # Adjust the path as needed
    filename=$(basename "$file" .txt)  # Get the filename without extension

    # Extract relevant data using awk
    awk -v fname="$filename" 'NR > 2 {  # Skip the first three lines (header)
        printf "%s\t%s\t%s\t%s\t%s\n", fname, $1, $2, $3, $5
    }' "$file" >> "tmp.txt"
done
cat tmp.txt | grep -v "Hit_no_genomes" > "$output_file"
rm tmp.txt
echo "Combined output written to $output_file"

## bash cov.sh /home/ubuntu/output/OH-VH00648-240920_TB