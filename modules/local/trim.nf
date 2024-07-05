process TRIMM {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      '':
      'staphb/trimmomatic:0.39'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*trimmed*.fastq.gz"),       emit: reads
    path("*log*"),                              emit: log

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    trimmomatic_minlen="40"
    trimmomatic_window_size="4"
    trimmomatic_quality_trim_score="15"
    threads="4"

    trimmomatic PE \
    -threads \${threads} \
    ${reads[0]} ${reads[1]} \
    -baseout ${prefix}.fastq.gz \
    -trimlog ${prefix}.trim.log \
    -summary ${prefix}.trim.stats.txt \
    LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:\${trimmomatic_window_size}:\${trimmomatic_quality_trim_score} \
    MINLEN:\${trimmomatic_minlen} 2> ${prefix}.trim.err

    # remove length=#
    zcat ${prefix}_1P.fastq.gz | sed 's/ length=[0-9]*//' | gzip > ${prefix}.trimmed.R1.fastq.gz
    zcat ${prefix}_2P.fastq.gz | sed 's/ length=[0-9]*//' | gzip > ${prefix}.trimmed.R2.fastq.gz

    # zip log for space
    gzip ${prefix}.trim.log
    """

    stub:
    """
    touch ${prefix}.R1.fastq.gz
    touch ${prefix}.R2.fastq.gz
    """
}