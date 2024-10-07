process QUAST {
    tag "$meta.id"
    label 'process_low'
    container 'staphb/quast:5.0.2'

    input:
        tuple val(meta), path(assembly)
        min_contig_length

    output:
        path "${prefix}_report.tsv",        emit: quast_report

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
        """
        quast.py ${assembly} \
            -o . \
            --min-contig ${min_contig_length}
        
        mv report.tsv ${prefix}_report.tsv

        bash quast_postprocess.sh ${prefix}_report.tsv
        """
}
