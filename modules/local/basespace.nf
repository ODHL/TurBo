process BASESPACE {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*bs.*fastq.gz"),       emit: reads

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ~/tools/basespace download biosample -n ${prefix}
    
    # move final files
    mv ${prefix}*ds*/*R1* ${prefix}.bs.R1.fastq.gz
    mv ${prefix}*ds*/*R2* ${prefix}.bs.R2.fastq.gz
    """

    stub:
    """
    touch ${prefix}.R1.fastq.gz
    touch ${prefix}.R2.fastq.gz
    """
}