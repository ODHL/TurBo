process CLOCKWORK_REMOVE {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'ghcr.io/iqbal-lab-org/clockwork:latest'}"
        //'dbest/clockwork:v1.0.0'}"

    input:
        tuple val(meta), path(sam)
        path(clockMeta)

    output:
        tuple val(meta), path("*fq.gz"),                    emit: deconFq
        tuple val(meta), path("*stats*"),                   emit: stats

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    clockwork remove_contam \
        ${clockMeta} \
        ${sam} \
        ${prefix}_clockwork_decontamination_stats.txt \
        ${prefix}_decontam_1.fq.gz\
        ${prefix}_decontam_2.fq.gz
 
    """

    stub:
    """
    for EXT in html txt png; do
        touch ${meta.id}.trimmed_screen.\${EXT}
    done
    """
}