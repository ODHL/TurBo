process ANNOTATEVCF {
    
    // Input parameters
    input:
    tuple val(meta),path(vcf)
    string name
    path parser
    path fOut
    boolean verbose = false


    output:
    path("*annotation.txt"),                emit: annoation

    // Command execution
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    java -Xmx4g -jar \
        "${annotator}" \
        -noLog NC_000962 \
        "${vcf}" > "${prefix}_annotation.txt"
    """
}
