#!/bin/bash

# Get the current timestamp for logging purposes
timestamp=$(date +'%Y%m%d%H%M%S')

# Create output directories and log files
mkdir Output_$timestamp
mkdir Output_$timestamp/QC
mkdir Output_$timestamp/tmp
mkdir Output_$timestamp/trimmomatic
mkdir Output_$timestamp/clockwork

# Copy config.yml to the current directory and update paths
cp $(dirname $0)/config.yml config.yml
sed -i "s#tools:#tools: $(dirname $0)#g" config.yml
sed -i "s#scripts:#scripts: $(dirname $0)#g" config.yml
sed -i "s#other:#other: $(dirname $0)#g" config.yml

# Read configuration from config.yml
nextflow=$(yaml get tools.nextflow < config.yml)
remove_contam=$(yaml get tools.remove_contam < config.yml)
ref_fasta=$(yaml get tools.ref_fasta < config.yml)
ref_metadata=$(yaml get tools.ref_metadata < config.yml)
trimmomatic=$(yaml get tools.trimmomatic < config.yml)
bwa=$(yaml get tools.bwa < config.yml)
samtools=$(yaml get tools.samtools < config.yml)
picard=$(yaml get tools.picard < config.yml)
gatk=$(yaml get tools.gatk < config.yml)
bcftools=$(yaml get tools.bcftools < config.yml)
bedtools=$(yaml get tools.bedtools < config.yml)
annotator=$(yaml get tools.annotator < config.yml)
parser=$(yaml get scripts.parser < config.yml)
creater=$(yaml get scripts.creater < config.yml)
interpreter=$(yaml get scripts.interpreter < config.yml)
included=$(yaml get scripts.included < config.yml)
reported=$(yaml get scripts.reported < config.yml)
stats_estimator=$(yaml get scripts.stats_estimator < config.yml)
genome_stats_estimator=$(yaml get scripts.genome_stats_estimator < config.yml)
bedlist_amp=$(yaml get scripts.bedlist_amp < config.yml)
bedlist_one=$(yaml get scripts.bedlist_one < config.yml)
bedlist_two=$(yaml get scripts.bedlist_two < config.yml)
bedstruct=$(yaml get scripts.bedstruct < config.yml)
target_estimator=$(yaml get scripts.target_cov_estimator < config.yml)
genome_cov_estimator=$(yaml get scripts.genome_cov_estimator < config.yml)
mutationloci=$(yaml get scripts.mutationloci < config.yml)
lineage_parser=$(yaml get scripts.lineage_parser < config.yml)
lineages=$(yaml get scripts.lineages < config.yml)
structparser=$(yaml get scripts.struct_parser < config.yml)
create_report=$(yaml get scripts.create_report < config.yml)

# Create log files
log=Output_$timestamp/$1.log
qlog=Output_$timestamp/QC/QC.log
touch $log
touch $qlog

# Write arguments to the log file
echo "$@" >> $log
echo "" >> $log

# Clockwork decontamination
echo "Performing clockwork decontamination."
if [ ! -z "$6" ]; then
    $nextflow run $remove_contam --ref_fasta $ref_fasta --ref_metadata_tsv $ref_metadata --reads_in1 $1 --reads_in2 $6 --outprefix Output_$timestamp/clockwork/$4 --mapping_threads $8
    input=Output_$timestamp/clockwork/$4.remove_contam.1.fq.gz
    input2=Output_$timestamp/clockwork/$4.remove_contam.2.fq.gz
fi

# QC Trimmomatic
echo "Performing trimmomatic trimming."
if [ ! -z "$6" ]; then
    java -jar $trimmomatic PE -threads $8 -trimlog Output_$timestamp/trimmomatic/trimLog.txt $1 $6 Output_$timestamp/trimmomatic/$4_paired_1.fastq.gz Output_$timestamp/trimmomatic/$4_unpaired_1.fastq.gz Output_$timestamp/trimmomatic/$4_paired_2.fastq.gz Output_$timestamp/trimmomatic/$4_unpaired_2.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40
    rm Output_$timestamp/trimmomatic/$4_unpaired_1.fastq.gz
    rm Output_$timestamp/trimmomatic/$4_unpaired_2.fastq.gz
    input=Output_$timestamp/trimmomatic/$4_paired_1.fastq.gz
    input2=Output_$timestamp/trimmomatic/$4_paired_2.fastq.gz
else
    java -jar $trimmomatic SE -threads $8 -trimlog Output_$timestamp/trimmomatic/trimLog.txt $1 Output_$timestamp/trimmomatic/$4_paired.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40
    input=Output_$timestamp/trimmomatic/$4_paired.fastq.gz
fi

# Align reads using BWA
echo "Running BWA."
mkdir -p Output_$timestamp/tmp/bwa/index
cp $2 Output_$timestamp/tmp/bwa/index/ref.fa
bwa index Output_$timestamp/tmp/bwa/index/ref.fa
java -jar $picard CreateSequenceDictionary R=Output_$timestamp/tmp/bwa/index/ref.fa O=Output_$timestamp/tmp/bwa/index/ref.dict
samtools faidx Output_$timestamp/tmp/bwa/index/ref.fa
if [ ! -z "$6" ]; then
    bwa mem -t $8 -R "@RG\\tID:$4\\tSM:$4\\tPL:ILLUMINA" Output_$timestamp/tmp/bwa/index/ref.fa $input $input2 > Output_$timestamp/tmp/bwa/bwa.sam
else
    bwa mem -t $8 -R "@RG\\tID:$4\\tSM:$4\\tPL:ILLUMINA" Output_$timestamp/tmp/bwa/index/ref.fa $input > Output_$timestamp/tmp/bwa/bwa.sam
fi

# Filter alignment using GATK and Picard-Tools
echo "Filtering alignment with GATK and Picard-Tools."
mkdir -p Output_$timestamp/GATK
mkdir -p Output_$timestamp/SamTools
java -jar $picard SamFormatConverter INPUT=Output_$timestamp/tmp/bwa/bwa.sam VALIDATION_STRINGENCY=LENIENT OUTPUT=Output_$timestamp/GATK/GATK.bam
java -Xmx8g -Djava.io.tmpdir=Output_$timestamp/tmp -jar $picard SortSam INPUT=Output_$timestamp/GATK/GATK.bam SORT_ORDER=coordinate OUTPUT=Output_$timestamp/GATK/GATK_s.bam VALIDATION_STRINGENCY=LENIENT TMP_DIR=Output_$timestamp/tmp
java -Xmx8g -jar $picard MarkDuplicates INPUT=Output_$timestamp/GATK/GATK_s.bam OUTPUT=Output_$timestamp/GATK/GATK_sdr.bam METRICS_FILE=Output_$timestamp/GATK/MarkDupes.metrics ASSUME_SORTED=true REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT
java -Xmx8g -jar $picard BuildBamIndex INPUT=Output_$timestamp/GATK/GATK_sdr.bam VALIDATION_STRINGENCY=LENIENT
samtools view -bhF 4 -o Output_$timestamp/$4_sdrcsm.bam Output_$timestamp/GATK/GATK_sdr.bam
java -Xmx8g -jar $picard BuildBamIndex INPUT=Output_$timestamp/$4_sdrcsm.bam VALIDATION_STRINGENCY=LENIENT
rm -r Output_$timestamp/tmp

# Run genome coverage statistics
echo "Running genome coverage statistics."
samtools depth -a Output_$timestamp/$4_sdrcsm.bam > Output_$timestamp/SamTools/coverage.txt
bedtools coverage -abam Output_$timestamp/$4_sdrcsm.bam -b $bedlist_one > Output_$timestamp/SamTools/bed_1_coverage.txt
bedtools coverage -abam Output_$timestamp/$4_sdrcsm.bam -b $bedlist_two > Output_$timestamp/SamTools/bed_2_coverage.txt
sort -nk 6 Output_$timestamp/SamTools/bed_1_coverage.txt > Output_$timestamp/SamTools/bed_1_sorted_coverage.txt
sort -nk 6 Output_$timestamp/SamTools/bed_2_coverage.txt > Output_$timestamp/SamTools/bed_2_sorted_coverage.txt
python $target_estimator $bedlist_amp Output_$timestamp/SamTools/coverage.txt $4 > Output_$timestamp/$4_target_region_coverage.txt
sort -nk 3 Output_$timestamp/$4_target_region_coverage.txt > Output_$timestamp/$4_target_region_coverage_sorted.txt
python $genome_stats_estimator Output_$timestamp/SamTools/$4_genome_stats.txt Output_$timestamp/SamTools/coverage.txt $4
python $genome_cov_estimator Output_$timestamp/SamTools/bed_1_sorted_coverage.txt Output_$timestamp/SamTools/coverage.txt $4 > Output_$timestamp/SamTools/genome_region_coverage_1.txt
python $genome_cov_estimator Output_$timestamp/SamTools/bed_2_sorted_coverage.txt Output_$timestamp/SamTools/coverage.txt $4 > Output_$timestamp/SamTools/genome_region_coverage_2.txt
cat Output_$timestamp/SamTools/genome_region_coverage_1.txt Output_$timestamp/SamTools/genome_region_coverage_2.txt > Output_$timestamp/SamTools/genome_region_coverage.txt
sort -nk 3 Output_$timestamp/SamTools/genome_region_coverage.txt > Output_$timestamp/$4_genome_region_coverage_sorted.txt
sed -i '/^#/d' Output_$timestamp/$4_genome_region_coverage_sorted.txt

# Variant calling and annotation
echo "Variant calling and annotation."
samtools mpileup -f Output_$timestamp/tmp/bwa/index/ref.fa Output_$timestamp/$4_sdrcsm.bam > Output_$timestamp/$4.pileup
bcftools call -c -O v -o Output_$timestamp/$4.vcf Output_$timestamp/$4.pileup
perl $annotator -v Output_$timestamp/$4.vcf -b Output_$timestamp/$4.bam -f Output_$timestamp/tmp/bwa/index/ref.fa -d $included -r $reported -s $stats_estimator > Output_$timestamp/$4_report.txt

# Parse results
echo "Parsing results."
perl $parser Output_$timestamp/$4_report.txt Output_$timestamp/tmp/bwa/index/ref.fa > Output_$timestamp/$4_mutations.txt
perl $creater Output_$timestamp/$4_mutations.txt > Output_$timestamp/$4_mutation_summary.txt
perl $interpreter Output_$timestamp/$4_mutation_summary.txt Output_$timestamp/$4.bed > Output_$timestamp/$4_interpreted_mutations.txt

# Analyze lineages
echo "Analyzing lineages."
perl $lineage_parser $lineages > Output_$timestamp/$4_parsed_lineages.txt

# Clean up
rm config.yml
