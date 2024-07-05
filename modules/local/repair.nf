process REPAIR {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'staphb/bbtools:39.01'}"
    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("*repaired*gz"),      emit: repaired_reads
        path("*repair.log"),                         emit: log

    script:
    def maxmem = task.memory.toGiga()-(task.attempt*12) // keep heap mem low so and rest of mem is for java expansion.
    """
    maxmem=\$(echo \"$maxmem GB\"| sed 's/ GB/g/g' | sed 's/-//g')
    set -ex
    repair.sh \\
        -Xmx\$maxmem \\
        in1=${reads[0]} \
        in2=${reads[1]} \
        out1=${meta.id}_repaired.R1.fastq.gz \
        out2=${meta.id}_repaired.R2.fastq.gz \
        outs=${meta.id}_singletons.fastq.gz

    cp .command.log ${meta.id}.repair.log
    """

    stub:
    """
        touch ${meta.id}_repaired.R1.fastq.gz
        touch ${meta.id}_repaired.R2.fastq.gz
        touch ${meta.id}_singletons.fastq.gz
    """
}