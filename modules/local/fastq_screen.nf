process FASTQ_SCREEN {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastq-screen:0.15.3--pl5321hdfd78af_0':
        'dbest/fastq_screen:v0.15.3'}"

    input:
        tuple val(meta), path(fastq)
        path(conf)
        path(contamZipped)

    output:
        tuple val(meta), path("*html"),     emit: screen

    script:
    """
    # unzip index and contam files in new dir
    tarFile=\$(readlink -f ${contamZipped})
    tar -xvf \${tarFile}

    fastq_screen \
      --force \
      --threads $task.cpus \
      --outdir . \
      --aligner "bwa" \
      --conf $conf \
      --subset 100000 \
      --nohits $fastq
    """

    stub:
    """
    for EXT in html txt png; do
        touch ${meta.id}.trimmed_screen.\${EXT}
    done
    """
}