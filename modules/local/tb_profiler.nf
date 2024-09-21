process TB_PROFILER {
    tag "tb_profiler"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'staphb/tbprofiler:4.4.2'}"

    input:
        tuple val(meta), path(reads)
        path(report)
        val(mapper)
        val(caller)
        val(min_depth)
        val(min_af)
        val(min_af_pred)
        val(cov_frac_threshold)
        val(threads)

    output:
        tuple val(meta), path("*/*.results.json"),            emit: results_json
        tuple val(meta), path("*/*.results.csv"),             emit: results_csv
        tuple val(meta), path("*/*.results.txt"),             emit: results_txt
        tuple val(meta), path("*/*.results.docx"),            emit: results_docx
        tuple val(meta), path("*/*.bam"),                     emit: bam
        tuple val(meta), path("*/*.bam.bai"),                 emit: bai
        tuple val(meta), path("*/*.targets.csq.vcf"),         emit: vcf

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tb-profiler profile \
        --read1 ${reads[0]} \
        --read2 ${reads[1]} \
        --threads $task.cpus \
        --prefix ${prefix} \
        --mapper ${mapper} \
        --caller ${caller} \
        --min_depth ${min_depth} \
        --af ${min_af} \
        --reporting_af ${min_af_pred} \
        --coverage_fraction_threshold ${cov_frac_threshold} \
        --csv \
        --txt \
        --docx $report

    gunzip vcf/${prefix}.targets.csq.vcf.gz
    """

    stub:
    """
    for EXT in html txt png; do
        touch ${meta.id}.trimmed_screen.\${EXT}
    done
    """
}