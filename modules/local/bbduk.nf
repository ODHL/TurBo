process BBDUK {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'staphb/bbtools:39.01'}"
    input:
        tuple val(meta), path(reads)
        path(contamZipped)

    output:
        tuple val(meta), path("*clean*gz"),      emit: repaired_reads
        path("*repair.log"),                     emit: log
        path("*stats*"),                         emit: stats

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maxmem = task.memory.toGiga()-(task.attempt*12)
    """
    #maxmem=\$(echo \"$maxmem GB\"| sed 's/ GB/g/g' | sed 's/-//g')
    maxmem="7g"
    cpus=$task.cpus
    
    # unzip index and contam files in new dir
    tarFile=\$(readlink -f ${contamZipped})
    tar -xvf \${tarFile}

    repair.sh \\
        -Xmx\$maxmem \\
        in1=${reads[0]} \
        in2=${reads[1]} \
        out1=${prefix}_repaired.R1.fastq.gz \
        out2=${prefix}_repaired.R2.fastq.gz \
        outs=${prefix}_singletons.fastq.gz

    # adapters
    bbduk.sh -Xmx\$maxmem threads=\$cpus \
        in1=${prefix}_repaired.R1.fastq.gz \
        in2=${prefix}_repaired.R2.fastq.gz \
        out1=${prefix}_no_adapter_1.fastq.gz \
        out2=${prefix}_no_adapter_2.fastq.gz \
        ref=/opt/bbmap/resources/adapters.fa \
        stats=${prefix}_adapters.stats.txt \
        ktrim=r k=23 mink=11 hdist=1 tpe tbo

    # phix
    bbduk.sh -Xmx\$maxmem threads=\$cpus \
        in1=${prefix}_no_adapter_1.fastq.gz \
        in2=${prefix}_no_adapter_2.fastq.gz \
        out1=${prefix}_no_phix_1.fastq.gz \
        out2=${prefix}_no_phix_2.fastq.gz \
        outm=${prefix}_matched_phix.fq.gz \
        ref=/opt/bbmap/resources/phix174_ill.ref.fa.gz \
        stats=${prefix}_phix.stats.txt \
        k=31 hdist=1

    # covid
    bbduk.sh -Xmx\$maxmem threads=\$cpus \
        in1=${prefix}_no_phix_1.fastq.gz \
        in2=${prefix}_no_phix_2.fastq.gz \
        out1=${prefix}_no_covid_1.fastq.gz \
        out2=${prefix}_no_covid_2.fastq.gz \
        outm=${prefix}_matched_covid.fq.gz \
        ref=/opt/bbmap/resources/Covid19_ref.fa \
        stats=${prefix}_Covid19.stats.txt \
        k=31 hdist=1

    # polyA
    bbduk.sh -Xmx\$maxmem threads=\$cpus \
        in1=${prefix}_no_covid_1.fastq.gz \
        in2=${prefix}_no_covid_2.fastq.gz \
        out1=${prefix}_no_polyA_1.fastq.gz \
        out2=${prefix}_no_polyA_2.fastq.gz \
        outm=${prefix}_matched_polyA.fq.gz \
        ref=/opt/bbmap/resources/polyA.fa.gz \
        stats=${prefix}_polyA.stats.txt \
        k=31 hdist=1

    # ecoli
    CONTAM=contamination/Escherichia_coli_gca_001606525.ASM160652v1.dna.toplevel.fa.gz
    bbduk.sh -Xmx\$maxmem -Xms4g threads=\$cpus \
        in1=${prefix}_no_polyA_1.fastq.gz \
        in2=${prefix}_no_polyA_2.fastq.gz \
        out1=${prefix}_no_Ecoli_1.fastq.gz \
        out2=${prefix}_no_Ecoli_2.fastq.gz \
        outm=${prefix}_matched_Ecoli.fq.gz \
        ref=\${CONTAM} \
        stats=${prefix}_Ecoli.stats.txt \
        k=10 hdist=1

    # NZ_CP    
    CONTAM=contamination/NZ_CP008889.1.fa.gz
    bbduk.sh -Xmx\$maxmem threads=\$cpus \
        in1=${prefix}_no_Ecoli_1.fastq.gz \
        in2=${prefix}_no_Ecoli_2.fastq.gz \
        out1=${prefix}_no_NZ_CP008889_1.fastq.gz \
        out2=${prefix}_no_NZ_CP008889_2.fastq.gz \
        outm=${prefix}_matched_NZ_CP008889.fq.gz \
        ref=\${CONTAM} \
        stats=${prefix}_NZ_CP008889.stats.txt \
        k=31 hdist=1

    # plasmid    
    CONTAM=contamination/NZ_CP008889.1_plasmid.fa.gz
    bbduk.sh -Xmx\$maxmem threads=\$cpus \
        in1=${prefix}_no_NZ_CP008889_1.fastq.gz \
        in2=${prefix}_no_NZ_CP008889_2.fastq.gz \
        out1=${prefix}_clean_1.fastq.gz \
        out2=${prefix}_clean_2.fastq.gz \
        outm=${prefix}_matched_NZ_CP008889_plasmid.fq.gz \
        ref=\${CONTAM} \
        stats=${prefix}_NZ_CP008889_plasmid.stats.txt \
        k=31 hdist=1
    
    # cleanup
    rm -frv ~{contamination}
    rm -fv ./*no_adapter*.fastq.gz
    rm -fv ./*no_phix*.fastq.gz
    rm -fv ./*no_covid*.fastq.gz
    rm -fv ./*no_polyA*.fastq.gz
    rm -fv ./*no_Ecoli*.fastq.gz
    rm -fv ./*no_NZ_*.fastq.gz

    cp .command.log ${meta.id}.repair.log
    """

    stub:
    """
        touch ${meta.id}_clean_1.fastq.gz
        touch ${meta.id}_clean_2.fastq.gz
        touch ${meta.id}.repair.log
        touch ${meta.id}_NZ_CP008889_plasmid.stats.txt
    """
}