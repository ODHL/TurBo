process CENTRIFUGE_DB {
    tag "create_db"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'dbest/centrifuge:v1.0.4'}"

    output:
        path("*sorted.tsv"),               emit: dbs

    script:
    """
    wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz
    gunzip nt.gz && mv -v nt nt.fa

    # Get mapping file
    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gitaxid_nucl.dmp.gz
    gunzip -c gi_taxid_nucl.dmp.gz | sed "s/^/gi|/" > gi_taxid_nucl.map

    # build index using 16 cores and a small bucket size, which will require less memory
    centrifuge-build -p 16 --bmax 1342177280 --conversion-table gi_taxid_nucl.map \
        --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \ 
        nt.fa nt
    """

    stub:
    """
    for EXT in html txt png; do
        touch ${meta.id}.trimmed_screen.\${EXT}
    done
    """
}