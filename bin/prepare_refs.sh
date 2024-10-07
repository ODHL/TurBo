#!/usr/bin bash
# https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly

################################################################################################################
## Usage
################################################################################################################
if [ "$#" -ne 2 ]; then
    echo "usage: $0 <ntm_accessions_file> <contamination_file"
    exit 1
fi
# Example: bash prepare_refs.sh ntm_test.tsv contam.csv

################################################################################################################
# Inputs
ntm_accessions_file=$1
species_file=$2

# set path
eutilsPath="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text"

# set files
merged_fasta=merged.fasta
if [[ -f $merged_fasta ]]; then rm $merged_fasta; fi
touch $merged_fasta

merged_meta=metadata_tmp.tsv
merged_meta_final=metadata.tsv
if [[ -f $merged_meta ]]; then rm $merged_meta; fi
touch $merged_meta
################################################################################################################
# Main
declare -a accessions
while IFS="\n" read -r ntm; do
    accessions+=("$ntm")
    echo "$ntm"

    if [[ -f ${ntm}.fa ]]; then
        echo "Skipping $ntm"
    else
        if [[ $ntm =~ ^[A-Z]{4}[0-9]{8}$ ]]; then
            wget -O ${ntm}.fa "${eutilsPath}&id=${ntm}"
        else
            wget -O ${ntm}.fa "${eutilsPath}&id=${ntm}"
        fi
    fi
    # grab headerline
    header=`cat ${ntm}.fa | grep ">" | cut -f1 -d"," | sed "s/ /|/" | sed "s/ /_/g"`

    # add information to fasta
    cat ${ntm}.fa >> $merged_fasta

    # add information to metadata
    echo -e "contam_group\t1\t$header" >> $merged_meta

    # cleanup
    rm ${ntm}.fa
done < $ntm_accessions_file


# Read the CSV file into the arrays
while IFS=, read -r sampleID category url; do
    # download
    filename=$(basename "$url")
    wget -O $filename "$url"

    # check if it's a bz file
    if [[ $filename =~ ".bz2" ]]; then
        scaffolds=($(grep "scaffolds.fa.bz2" <(tar --list -f "$filename")))
        tar -O -x -f "$filename" "${scaffolds[0]}" | bunzip2 -c > $sampleID.fa
        rm $filename
        filename=$sampleID.fa
    fi

    # add information to fasta
    cat ${filename} >> $merged_fasta

    # add information to metadata
    if [[ $category == "reference" ]]; then
        awk '{print "reference\t0\t" $0}' <(grep ">" "${filename}" | cut -f1 -d"," | sed "s/ /|/" | sed "s/ /_/g") >> $merged_meta
    else
        awk '{print "is_contam\t1\t" $0}' <(grep ">" "${filename}" | cut -f1 -d"," | sed "s/ /|/" | sed "s/ /_/g") >> $merged_meta
    fi

    # cleanup
    rm ${filename}

done < $species_file

# cleanup metadata file
cat $merged_meta | grep ">" | sed "s/>//g" | sort | uniq > $merged_meta_final

# zip fasta
gzip $merged_fasta