process CLOCKWORK_REFS {
    tag "processRefs"
    label 'process_single'

    input:
    path(contam_list)
    path(ntm_list)

    output:
    path("merged.fasta.gz")          , emit: fasta
    path("metadata.tsv")             , emit: metadata

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    """
    bash prepare_refs.sh $ntm_list $contam_list

    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    touch merged_fasta.gz
    touch merged_metadata.tsv    
    """
}