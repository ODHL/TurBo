process TBPROFILER {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'staphb/tbprofiler:6.2.1'}"

    input:
    tuple val(meta), path(reads)
    path(template)

    output:
    tuple val(meta), path("results/*results.csv"),          emit: full_report_csv

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tb-profiler profile \
      --threads ${task.cpus} \
      --platform ${params.platform} \
      --mapper ${params.mapper} \
      --caller ${params.caller} \
      --depth ${params.min_depth} \
      --af ${params.min_af_used_for_calling} \
      --read1 ${reads[0]} \
      --read2 ${reads[1]} \
      --prefix ${prefix} \
      --docx $template \
      --csv \
      --call_whole_genome
    """
}
