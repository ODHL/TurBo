process CLOCKWORK_PREP {
    tag "prep"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'ghcr.io/iqbal-lab-org/clockwork:latest'}"
        //'dbest/clockwork:v1.0.0'}"

    input:
        path(fa)

    when:
        params.clockworkCallVariants

    output:
        path("ref/ref*"),                    emit: indexedRefs

    script:
    """
    clockwork reference_prepare \
      --outdir ref \
      $fa

    """

    stub:
    """
    for EXT in html txt png; do
        touch trimmed_screen.\${EXT}
    done
    """
}