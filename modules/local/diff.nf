process make_mask_and_diff {
    
    // Define input parameters
    input:
    path bam
    boolean force_diff = false
    boolean histograms = false
    float max_ratio_low_coverage_sites_per_sample
    int min_coverage_per_site
    path tbmf
    path vcf

    // Runtime attributes
    int addldisk = 10
    int cpu = 8
    int retries = 1
    int memory = 16 GB
    int preempt = 1

    // Determine the output paths and sizes
    String basename_bam = file(bam).baseName
    String basename_vcf = file(vcf).baseName
    int finalDiskSize = Math.ceil(bam.size() / (1024 * 1024 * 1024)) * 2 + Math.ceil(vcf.size() / (1024 * 1024 * 1024)) * 2 + addldisk

    // Command to execute
    script:
    """
    set -euxo pipefail
    start=\$(date +%s)

    # Determine the mask file
    if [[ ! -e "${tbmf}" ]]; then
        mask="/mask/R00000039_repregions.bed"
    else
        mask="${tbmf}"
    fi

    echo "Copying bam..."
    cp "${bam}" .

    echo "Sorting bam..."
    samtools sort -u "${basename_bam}.bam" > sorted_u_"${basename_bam}.bam"

    echo "Calculating coverage..."
    bedtools genomecov -ibam sorted_u_"${basename_bam}.bam" -bga | \
        awk '\$4 < ${min_coverage_per_site}' > \
        "${basename_bam}_below_${min_coverage_per_site}x_coverage.bedgraph"

    if [[ "${histograms}" == "true" ]]; then
        echo "Generating histograms..."
        bedtools genomecov -ibam sorted_u_"${basename_bam}.bam" > histogram.txt
    fi

    echo "Pulling diff script..."
    wget https://raw.githubusercontent.com/aofarrel/parsevcf/1.1.8/vcf_to_diff_script.py
    echo "Running script..."
    python3 vcf_to_diff_script.py -v "${vcf}" \
        -d . \
        -tbmf "${mask}" \
        -bed "${basename_bam}_below_${min_coverage_per_site}x_coverage.bedgraph" \
        -cd "${min_coverage_per_site}"

    # Check for low coverage sites
    this_files_info=\$(awk -v file_to_check="${basename_vcf}.diff" '\$1 == file_to_check' "${basename_vcf}.report")
    if [[ -n "\$this_files_info" ]]; then
        echo "\$this_files_info" > temp
        amount_low_coverage=\$(cut -f2 temp)
        percent_low_coverage=\$(echo "\$amount_low_coverage * 100" | bc)
        echo "\$percent_low_coverage percent of ${basename_vcf} is below ${min_coverage_per_site}x coverage."

        is_bigger=\$(echo "\$amount_low_coverage > ${max_ratio_low_coverage_sites_per_sample}" | bc)
        if [[ \$is_bigger == 0 ]]; then
            echo "PASS" >> ERROR
        else
            if [[ "${force_diff}" == "false" ]]; then
                rm "${basename_vcf}.diff"
            fi
            pretty_percent=\$(printf "%0.2f" "\$percent_low_coverage")
            echo "FAILURE - \$pretty_percent% is above ${max_ratio_low_coverage_sites_per_sample}% cutoff"
            echo "VCF2DIFF_\$pretty_percent_PCT_BELOW_${min_coverage_per_site}x_COVERAGE" >> ERROR
        fi
    fi

    end=\$(date +%s)
    seconds=\$(echo "\$end - \$start" | bc)
    minutes=\$(echo "\$seconds / 60" | bc)
    echo "Finished in about \$minutes minutes (\$seconds sec)"
    ls -lha
    """

    // Define output
    output:
    file("${basename_bam}_below_${min_coverage_per_site}x_coverage.bedgraph") into mask_file
    file("${basename_vcf}.diff") optional into diff
    file("${basename_vcf}.report") optional into report
    file("histogram.txt") optional into histogram
    string(read_string("ERROR")) into errorcode

    // Define resource requirements
    resources:
    cpus cpu
    memory memory
    disk "${finalDiskSize} GB"
    retries retries
    preemptible preempt

    // Metadata
    meta:
    author "Lily Karim (WDLization by Ash O'Farrell)"
}
