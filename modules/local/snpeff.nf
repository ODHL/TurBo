process SNPEFF {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'dbest/snpeff:v5.2a'}"

    input:
        tuple val(meta), path(bcf)

    output:
        tuple val(meta), path("*snpeff.vcf"),       emit: outputVcf
        tuple val(meta), path("*snpEff_summary.csv"),       emit: snpEff_summary_csv
        tuple val(meta), path("*snpEff_summary.html"),       emit: snpEff_summary_html

    script:
    def prefix              = task.ext.prefix ?: "${meta.id}"
    def hgvs                = params.hgvs ? "-hgvs" : "-noHgvs"
    def lof                 = params.lof ? "-lof" : "-noLof"
    def noDownstream        = params.noDownstream ? "-no-downstream" : ""
    def noIntergenic        = params.noIntergenic ? "-no-intergenic" : ""
    def noShiftHgvs         = params.noShiftHgvs ? "-noShiftHgvs" : ""
    """
    mkdir -p "$(dirname ~{outputPath})"
    
    unzip ~{dataDir}

    File vcf
      File config
      File dataDir
      String genome
      Int? upDownStreamLen
      String memory = "9G"

    snpEff ann \
        -verbose \
        -noDownload \
        -noLog \
        -stats ${meta.id}_snpEff_summary.html \
        -csvStats ${meta.id}_snpEff_summary.csv \
        -config ~{config} \
        -dataDir "$PWD"/data \
        $hgvs \
        $lof \
        $noDownstream \
        $noIntergenic \
        $noShiftHgvs
        ~{"-upDownStreamLen " + upDownStreamLen} \
        ~{genome} \
        ~{vcf} \
        > ${meta.id}_snpeff.vcf
    
    rm -r "$PWD"/data
    """

    stub:
    """
    for EXT in html txt png; do
        touch ${meta.id}.vcf
    done
    """
}