process CLOCKWORK {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'dbest/clockwork:v1.0.0'}"

    input:
        tuple val(meta), path(reads)
        path(clockContam)
        path(clockMeta)

    output:
        tuple val(meta), path("*decontamination.bam"),       emit: deconBam
        tuple val(meta), path("*_clockwork_cleaned*"),       emit: reads
        path("*stats*"),                                     emit: stats

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    clockwork map_reads \
        --threads ${task.cpus} \
        --unsorted_sam ${prefix} \
        ${clockContam} \
        ${prefix}_clockwork_decontamination.bam \
        ${reads[0]} \
        ${reads[1]}

    clockwork remove_contam \
        ${clockMeta} \
        ${prefix}_clockwork_decontamination.bam \
        ${prefix}_clockwork_decontamination_stats.txt \
        ${prefix}_clockwork_cleaned_1.fastq.gz \
        ${prefix}_clockwork_cleaned_2.fastq.gz
 
    """

    stub:
    """
    for EXT in html txt png; do
        touch ${meta.id}.trimmed_screen.\${EXT}
    done
    """
}