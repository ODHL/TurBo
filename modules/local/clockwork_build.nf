process BUILD_SIF {
    tag "buildSIF"
    label 'process_single'

    input:
    path(clockworkDEF)

    output:
    path("clockwork.sif")                 , emit: sif

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    """
    singularity build --fakeroot clockwork.sif $clockworkDEF
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    touch merged_fasta.gz
    touch merged_metadata.tsv    
    """
}