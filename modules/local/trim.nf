process TRIMM {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      '':
      'staphb/trimmomatic:0.39'}"

    input:
    tuple val(meta), path(reads)
    val(minlen)
    val(window_size)
    val(quality_trim_score)


    output:
    tuple val(meta), path("*trimmed*.fastq.gz"),       emit: reads
    path("*log*"),                              emit: log

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    trimmomatic PE \
      -threads $task.cpus \
      ${reads[0]} ${reads[1]} \
      -baseout ${prefix}.fastq.gz \
      -trimlog ${prefix}.trim.log \
      -summary ${prefix}.trim.stats.txt \
      LEADING:3 TRAILING:3 \
      SLIDINGWINDOW:${window_size}:${quality_trim_score} \
      MINLEN:${minlen} 2> ${prefix}.trim.err

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