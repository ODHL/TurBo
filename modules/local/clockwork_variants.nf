process CLOCKWORK_VARIANTS {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'ghcr.io/iqbal-lab-org/clockwork:latest'}"
        //'dbest/clockwork:v1.0.0'}"

    input:
        tuple val(meta), path(fq)
        path(indexedRefs)

    when:
        params.clockworkCallVariants

    output:
        tuple val(meta), path("*/final.vcf"),                               emit: vcf
        tuple val(meta), path("*/*bam"),                                    emit: sortedBam
        tuple val(meta), path("*/*bai"),                                    emit: sortedBai

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir H37Rv
    mv ref* H37Rv
    
    clockwork variant_call_one_sample \
        --sample_name ${prefix} \
        --keep_bam \
        H37Rv Var_call_${prefix} \
        ${fq[0]} ${fq[1]}
    """

    stub:
    """
    for EXT in html txt png; do
        touch ${meta.id}.trimmed_screen.\${EXT}
    done
    """
}