process CENTRIFUGE {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'dbest/centrifuge:v1.0.4'}"

    input:
        tuple val(meta), path(reads)
        path(centrifugeZipped)

    output:
        tuple val(meta), path("*summary.report.tsv"),       emit: summary
        tuple val(meta), path("*sorted.tsv"),               emit: classification

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # unzip index and contam files in new dir
    tarFile=\$(readlink -f ${centrifugeZipped})
    tar -xvf \${tarFile}
    
    # unzip reads
    gunzip -c ${reads[0]} > ${prefix}.R1.fastq
    gunzip -c ${reads[1]} > ${prefix}.R2.fastq

    centrifuge -x index/test \
        --threads $task.cpus \
        -q \
        -1 ${prefix}.R1.fastq -2 ${prefix}.R2.fastq \
        --report-file ${prefix}.centrifuge.summary.report.tsv \
        -S ${prefix}.centrifuge.classification.tsv

    (head -n1 ${prefix}.centrifuge.summary.report.tsv ; tail -n+2 ${prefix}.centrifuge.summary.report.tsv | sort -t \$'\t' -r -g -k7 ) > ${prefix}.centrifuge.summary.report.sorted.tsv
    """

    stub:
    """
    for EXT in html txt png; do
        touch ${meta.id}.trimmed_screen.\${EXT}
    done
    """
}