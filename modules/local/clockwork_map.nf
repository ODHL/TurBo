process CLOCKWORK_MAP {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'ghcr.io/iqbal-lab-org/clockwork:latest'}"
        //'dbest/clockwork:v1.0.0'}"

    input:
        tuple val(meta), path(reads)
        path(clockContam)

    output:
        tuple val(meta), path("*decontamination.sam"),      emit: deconSam

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # unzip file
    gzFile=\$(readlink -f ${clockContam})
    gunzip -c \${gzFile} > merged.fasta

    clockwork map_reads \
        --threads ${task.cpus} \
        --unsorted_sam ${prefix} \
        merged.fasta \
        ${prefix}_decontamination.sam \
        ${reads[0]} \
        ${reads[1]} 
    """

    stub:
    """
    for EXT in html txt png; do
        touch ${meta.id}.trimmed_screen.\${EXT}
    done
    """
}