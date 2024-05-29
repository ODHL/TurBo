process PREPSRR {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    path("*.gz"),                   emit: fastq

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "----SRA"
    prefetch $prefix
                
    echo "----FASTQ"
    cd $prefix
    fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip ${prefix}.sra

    echo "----cleaning"
    mv fastq/${prefix}_pass_1.fastq.gz ../${prefix}.R1.fastq.gz
    mv fastq/${prefix}_pass_2.fastq.gz ../${prefix}.R2.fastq.gz    
    """

    stub:
    """
    touch ${prefix}.R1.fastq.gz
    touch ${prefix}.R2.fastq.gz
    """
}