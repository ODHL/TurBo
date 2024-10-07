process PREP_DB {
    tag "prepareDB"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'ghcr.io/iqbal-lab-org/clockwork:latest'}"

    input:
    path(dbIni)

    output:
    path("dbv1.0")                 , emit: db

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    """
    clockwork make_empty_db $dbIni
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    touch db.ini
    """
}