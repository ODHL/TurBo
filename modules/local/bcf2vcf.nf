process BCF2VCF {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'staphb/bcftools:1.17'}"

    input:
        tuple val(meta), path(bcf)

    output:
        tuple val(meta), path("*vcf"),       emit: vcf

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools view ${bcf} > ${meta.id}.vcf
    """

    stub:
    """
    for EXT in html txt png; do
        touch ${meta.id}.vcf
    done
    """
}
