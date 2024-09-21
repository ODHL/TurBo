process CENTRIFUGE_DB {
    tag "centrifuge_db"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'dbest/centrifuge:v1.0.4':
        'quay.io/biocontainers/ngsutils:0.5.9'}"

    output:
        path("*cf"),               emit: cfs

    script:
    """
    # build index using 16 cores and a small bucket size, which will require less memory
    centrifuge-build \
        -p 16 \
        --bmax 1342177280 \
        --conversion-table gi_taxid_nucl.map \
    --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp nt.fa nt
    """

    stub:
    """
    for EXT in html txt png; do
        touch ${meta.id}.trimmed_screen.\${EXT}
    done
    """
}