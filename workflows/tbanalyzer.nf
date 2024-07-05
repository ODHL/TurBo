/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/local/fastqc'
include { TRIMM                  } from '../modules/local/trim'
include { REPAIR                 } from '../modules/local/repair'
include { FASTQ_SCREEN           } from '../modules/local/fastq_screen'
include { CENTRIFUGE             } from '../modules/local/centrifuge'
include { CENTRIFUGE_DB          } from '../modules/local/centrifuge_db'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_tbAnalyzer_pipeline'
include { CREATE_INPUT_CHANNEL    } from '../subworkflows/local/create_input_channel'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow tbAnalyzer {

    take:
    ch_reads // channel: 

    main:

    // prep empty channels
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // MODULE: Repair
    TRIMM (
        ch_reads
    )
    ch_trimmed=TRIMM.out.reads

    // MODULE: Repair
    REPAIR (
        ch_trimmed
    )
    ch_trimmed.view()

    //
    // MODULE: FASTQ Screen
    FASTQ_SCREEN (
        ch_trimmed,
        params.fastq_screen_configuration,
        params.contam_dir
    )

    //
    // MODULES: CENTRIFUGE DB
    CENTRIFUGE_DB()

    //
    // MODULES: CENTRIFUGE
    //centrifuge_files = Channel.fromPath(params.centrifuge_dir)
    CENTRIFUGE(
        ch_trimmed,
        CENTRIFUGE_DB.out.dbs
    )


    //   if ( run_centrifuge ) {
    //     call centrifuge.wf_centrifuge {
    //       input:
    //       read1 = task_repair.repaired_out1,
    //       read2 = task_repair.repaired_out2,
    //       samplename = samplename,
    //       indexFiles = indexFiles
    //     }
    //   } 


    // //
    // // Collate and save software versions
    // //
    // softwareVersionsToYAML(ch_versions)
    //     .collectFile(
    //         storeDir: "${params.outdir}/pipeline_info",
    //         name: 'nf_core_pipeline_software_mqc_versions.yml',
    //         sort: true,
    //         newLine: true
    //     ).set { ch_collated_versions }

    // //
    // // MODULE: MultiQC
    // //
    // ch_multiqc_config        = Channel.fromPath(
    //     "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    // ch_multiqc_custom_config = params.multiqc_config ?
    //     Channel.fromPath(params.multiqc_config, checkIfExists: true) :
    //     Channel.empty()
    // ch_multiqc_logo          = params.multiqc_logo ?
    //     Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
    //     Channel.empty()

    // summary_params      = paramsSummaryMap(
    //     workflow, parameters_schema: "nextflow_schema.json")
    // ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    // ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
    //     file(params.multiqc_methods_description, checkIfExists: true) :
    //     file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    // ch_methods_description                = Channel.value(
    //     methodsDescriptionText(ch_multiqc_custom_methods_description))

    // ch_multiqc_files = ch_multiqc_files.mix(
    //     ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    // ch_multiqc_files = ch_multiqc_files.mix(
    //     ch_methods_description.collectFile(
    //         name: 'methods_description_mqc.yaml',
    //         sort: true
    //     )
    // )

    // MULTIQC (
    //     ch_multiqc_files.collect(),
    //     ch_multiqc_config.toList(),
    //     ch_multiqc_custom_config.toList(),
    //     ch_multiqc_logo.toList()
    // )

    // emit:
    // multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    // versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// """ Clockwork Decontamination """
//     def runClockwork(self):
//         self.__ifVerbose("Performing clockwork decontamination.")
//         if self.paired:
//            self.__CallCommand('nextflow remove contamination', [self.__nextflow, 'run', self.__remove_contam, '--ref_fasta', self.__ref_fasta, '--ref_metadata_tsv', 
//                               self.__ref_metadata, '--reads_in1', self.input, '--reads_in2', self.input2, '--outprefix', self.clockwork + "/" + self.name,
//                               '--mapping_threads', self.__threads])
//            self.input   = self.clockwork + "/" + self.name + '.remove_contam.1.fq.gz'
//            self.input2  = self.clockwork + "/" + self.name + '.remove_contam.2.fq.gz'

//     """ QC Trimmomatic """
//     def runTrimmomatic(self):
//         self.__ifVerbose("Performing trimmomatic trimming.")
//         if self.paired:
//            self.__CallCommand('trimmomatic', ['java', '-jar', self.__trimmomatic, 'PE', '-threads', self.__threads, '-trimlog',
//                               self.trimmomatic + "/" + 'trimLog.txt', self.input, self.input2,
//                               self.trimmomatic + "/" + self.name + '_paired_1.fastq.gz', self.trimmomatic + "/" + self.name + '_unpaired_1.fastq.gz',
//                               self.trimmomatic + "/" + self.name + '_paired_2.fastq.gz', self.trimmomatic + "/" + self.name + '_unpaired_2.fastq.gz',
//                               'LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:4:15', 'MINLEN:40'])
//         else:
//            self.__CallCommand('trimmomatic', ['java', '-jar', self.__trimmomatic, 'SE', '-threads', self.__threads, 
//                               '-trimlog', self.trimmomatic + "/" + 'trimLog.txt',
//                               self.input, self.trimmomatic + "/" + self.name + '_paired.fastq.gz',
//                               'LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:4:15', 'MINLEN:40'])
//         if self.paired:
//            self.__CallCommand('rm', ['rm', self.trimmomatic + "/" + self.name + "_unpaired_1.fastq.gz",
//                               self.trimmomatic + "/" + self.name + "_unpaired_2.fastq.gz"])
//            self.input  = self.trimmomatic + "/" + self.name + "_paired_1.fastq.gz"
//            self.input2 = self.trimmomatic + "/" + self.name + "_paired_2.fastq.gz"
//         else:
//            self.input = self.trimmomatic + "/" + self.name + "_paired.fastq.gz"
    
//     """ Aligners """ 
//     def runBWA(self, bwa):
//         """ Align reads against the reference using bwa."""
//         self.__ranBWA = True
//         self.__ifVerbose("Running BWA.")
//         self.__logFH.write("########## Running BWA. ##########\n")
//         bwaOut = os.path.join(self.outdir, "bwa")
//         self.__CallCommand('mkdir', ['mkdir', '-p', bwaOut])
//         self.__ifVerbose("   Building BWA index.")
//         self.__bwaIndex(bwaOut + "/index")
//         self.__alnSam = bwaOut + "/bwa.sam"
//         self.__bwaLongReads(bwaOut)
//         self.__ifVerbose("") 
//         self.__processAlignment()
          
//     def __bwaIndex(self, out):
//         """ Make an index of the given reference genome. """ 
//         self.__CallCommand('mkdir', ['mkdir', '-p', out])
//         self.__CallCommand('cp', ['cp', self.reference, out + "/ref.fa"])
//         self.reference = out + "/ref.fa"
//         self.__CallCommand('bwa index', [self.__bwa, 'index', self.reference])
//         self.__CallCommand('CreateSequenceDictionary', ['java', '-jar', self.__picard, 
//                            'CreateSequenceDictionary', 'R='+self.reference,'O='+ out + "/ref.dict"])
//         self.__CallCommand('samtools faidx', ['samtools', 'faidx', self.reference ])

//     def __bwaLongReads(self, out):
//         """ Make use of bwa mem """
//         if self.paired:
//             self.__ifVerbose("   Running BWA mem on paired end reads.")
//             self.__CallCommand(['bwa mem', self.__alnSam], [self.__bwa, 'mem','-t',self.__threads,'-R', 
//                                "@RG\\tID:" + self.name + "\\tSM:" + self.name + "\\tPL:ILLUMINA", 
//                                 self.reference, self.input, self.input2])
//         else:
//             self.__ifVerbose("   Running BWA mem on single end reads.")
//             self.__CallCommand(['bwa mem', self.__alnSam], [self.__bwa, 'mem','-t', self.__threads, '-R', 
//                                "@RG\\tID:" + self.name + "\\tSM:" + self.name + "\\tPL:ILLUMINA", 
//                                 self.reference, self.input])       

//     def __processAlignment(self):
//         """ Filter alignment using GATK and Picard-Tools """
//         self.__ifVerbose("Filtering alignment with GATK and Picard-Tools.")
//         self.__logFH.write("########## Filtering alignment with GATK and Picard-Tools. ##########\n")
//         GATKdir = os.path.join(self.outdir, "GATK")
//         self.__CallCommand('mkdir', ['mkdir', '-p', GATKdir])
//         samDir  = os.path.join(self.outdir, "SamTools")
//         self.__CallCommand('mkdir', ['mkdir', '-p', samDir])

//         """ Convert SAM to BAM"""
//         if (self.__ranBWA):
//             self.__ifVerbose("   Running SamFormatConverter.")
//             self.__CallCommand('SamFormatConverter', ['java', '-Xmx4g', '-jar', self.__picard, 'SamFormatConverter',  
//                                'INPUT='+ self.__alnSam, 'VALIDATION_STRINGENCY=LENIENT', 
//                                'OUTPUT='+ GATKdir +'/GATK.bam', ])
//         else:
//             self.__CallCommand('cp', ['cp', self.__alnSam, GATKdir +'/GATK.bam'])


//         """ Run mapping Report and Mark duplicates using Picard-Tools"""
//         self.__ifVerbose("   Running SortSam.")
//         self.__CallCommand('SortSam', ['java', '-Xmx8g', '-Djava.io.tmpdir=' + self.tmp, '-jar', self.__picard, 'SortSam',  
//                            'INPUT='+ GATKdir +'/GATK.bam', 'SORT_ORDER=coordinate', 'OUTPUT='+ GATKdir +'/GATK_s.bam', 
//                            'VALIDATION_STRINGENCY=LENIENT', 'TMP_DIR=' + self.tmp])
//         self.__ifVerbose("   Running MarkDuplicates.")
//         self.__CallCommand('MarkDuplicates', ['java', '-Xmx8g', '-jar', self.__picard, 'MarkDuplicates',  
//                            'INPUT='+ GATKdir +'/GATK_s.bam', 'OUTPUT='+ GATKdir +'/GATK_sdr.bam',
//                            'METRICS_FILE='+ GATKdir +'/MarkDupes.metrics', 'ASSUME_SORTED=true', 
//                            'REMOVE_DUPLICATES=false', 'VALIDATION_STRINGENCY=LENIENT'])
//         self.__ifVerbose("   Running BuildBamIndex.")
//         self.__CallCommand('BuildBamIndex', ['java', '-Xmx8g', '-jar', self.__picard, 'BuildBamIndex',  
//                            'INPUT='+ GATKdir +'/GATK_sdr.bam', 'VALIDATION_STRINGENCY=LENIENT'])
//         self.__CallCommand(['samtools view', samDir + '/unmapped.txt'],['samtools', 'view', '-c', GATKdir +'/GATK_sdr.bam'])
      
//         """ Filter out unmapped reads """
//         self.__finalBam = self.fOut + '/'+ self.name + '_sdrcsm.bam'
//         self.__ifVerbose("   Running samtools view.")
//         self.__CallCommand('samtools view', ['samtools', 'view', '-bhF', '4', '-o', self.__finalBam, 
//                            GATKdir +'/GATK_sdr.bam'])
//         self.__ifVerbose("   Running BuildBamIndex.")
//         self.__CallCommand('BuildBamIndex', ['java', '-Xmx8g', '-jar', self.__picard, 'BuildBamIndex', 'INPUT='+ self.__finalBam, 
//                            'VALIDATION_STRINGENCY=LENIENT'])
//         self.__ifVerbose("")
//         self.__CallCommand('rm', ['rm', '-r', self.tmp])
//         self.__CallCommand(['samtools view', samDir + '/mapped.txt'],['samtools', 'view', '-c', self.__finalBam])

    
//     def runCoverage(self):
//         """ Run genome Coverage Statistics """

//         self.__ifVerbose("Running target Coverage Statistics")
//         samDir = self.outdir + "/SamTools"
//         i = datetime.now()

//         self.__CallCommand(['samtools depth', samDir + '/coverage.txt'],['samtools','depth', '-a', self.__finalBam])
//         #self.__CallCommand(['bedtools coverage', samDir + '/bed_amp_coverage.txt' ],
//                            #[self.__bedtools, 'coverage', '-abam', self.__finalBam, '-b', self.__bedlist_amp])
//         #self.__CallCommand(['sort', samDir + '/bed_amp_sorted_coverage.txt' ],['sort', '-nk', '6', samDir + '/bed_amp_coverage.txt'])
//         self.__CallCommand(['bedtools coverage', samDir + '/bed_1_coverage.txt' ],
//                            [self.__bedtools, 'coverage', '-abam', self.__finalBam, '-b', self.__bedlist_one])
//         self.__CallCommand(['bedtools coverage', samDir + '/bed_2_coverage.txt' ],
//                            [self.__bedtools, 'coverage', '-abam', self.__finalBam, '-b', self.__bedlist_two])
//         self.__CallCommand(['sort', samDir + '/bed_1_sorted_coverage.txt' ],['sort', '-nk', '6', samDir + '/bed_1_coverage.txt'])
//         self.__CallCommand(['sort', samDir + '/bed_2_sorted_coverage.txt' ],['sort', '-nk', '6', samDir + '/bed_2_coverage.txt'])
//         self.__CallCommand(['target region coverage estimator', samDir + '/target_region_coverage_amp.txt'],
//                             ['python', self.__target_estimator, self.__bedlist_amp, samDir + '/coverage.txt', self.name])
//         self.__CallCommand(['sort', self.fOut + "/" + self.name + '_target_region_coverage.txt' ],
//                            ['sort', '-nk', '3', samDir + '/target_region_coverage_amp.txt'])
//         self.__CallCommand(['genome stats estimator', samDir + "/" + self.name + '_genome_stats.txt'],
//                             ['python', self.__genome_stats_estimator, samDir + '/coverage.txt', self.name])
//         self.__CallCommand(['genome region coverage estimator', samDir + '/genome_region_coverage_1.txt'],
//                             ['python', self.__genome_coverage_estimator, samDir + '/bed_1_sorted_coverage.txt', samDir + '/coverage.txt', self.name])
//         self.__CallCommand(['genome region coverage estimator', samDir + '/genome_region_coverage_2.txt'],
//                             ['python', self.__genome_coverage_estimator, samDir + '/bed_2_sorted_coverage.txt', samDir + '/coverage.txt', self.name])
//         self.__CallCommand(['cat' , samDir + '/genome_region_coverage.txt'],
//                            ['cat', samDir + '/genome_region_coverage_1.txt', samDir + '/genome_region_coverage_2.txt'])
//         self.__CallCommand(['sort', self.fOut + "/" + self.name + '_genome_region_coverage.txt' ],
//                            ['sort', '-nk', '3', samDir + '/genome_region_coverage.txt'])
//         self.__CallCommand('sed',['sed', '-i', '1d', self.fOut + "/" + self.name + '_genome_region_coverage.txt'])
//         self.__CallCommand(['structural variant detector', self.fOut + "/" + self.name + '_structural_variants.txt'],
//                            ['python', self.__structparser, self.__bedstruct, self.fOut + "/" + self.name + '_genome_region_coverage.txt',
//                             samDir + '/coverage.txt', self.name])
//         self.__CallCommand(['stats estimator', self.fOut + "/" + self.name + '_stats.txt'],
//                            ['python', self.__stats_estimator, samDir + '/unmapped.txt', samDir + '/mapped.txt',
//                             self.fOut + "/" + self.name + '_target_region_coverage.txt', self.name, samDir + "/" + self.name + '_genome_stats.txt'])
//         statsOut = self.fOut + "/" + self.name + '_stats.txt'
//         fh20 = open(statsOut, 'r')
//         for lines in fh20:
//             if lines.startswith('Sample ID'):
//                continue
//             fields = lines.rstrip("\r\n").split("\t")
//             if float(fields[2]) < 90.0 or int(fields[3]) < 30 or float(fields[4]) < 90.0:
//                self.__logFH2.write(i.strftime('%Y/%m/%d %H:%M:%S') + "\t" + "Sample:" + "\t" + self.name + "\t" + "failed QC checks\n")
//                self.__CallCommand('rm', ['rm', '-r', self.outdir])
//                self.__CallCommand('rm', ['rm',  self.__finalBam])
//                self.__CallCommand('rm', ['rm', '-r', self.trimmomatic])
//                self.__CallCommand('rm', ['rm', '-r', self.clockwork])
//                self.__CallCommand('mv', ['mv', self.fOut, self.qlog])
//                sys.exit(1)


//     """ Callers """

//     def runGATK(self):
//         if os.path.isfile(self.__finalBam):
//             self.__ifVerbose("Calling SNPs/InDels with GATK.")
//             self.__logFH.write("########## Calling SNPs/InDels with Mutect2. ##########\n")
//             GATKdir = os.path.join(self.outdir, "GATK")
//             samDir = os.path.join(self.outdir, "SamTools")

//             """ Call SNPs/InDels with Mutect2 """
//             self.__ifVerbose("   Running Mutect2.")
//             self.__CallCommand('Mutect2', [self.__gatk, 'Mutect2',
//                                '-R', self.reference, '-I', self.__finalBam, '-O',  GATKdir +'/mutect.vcf',
//                                '--max-mnp-distance', '2','-L', self.__included])
//             self.__CallCommand('Mutect2', [self.__gatk, 'Mutect2',
//                                '-R', self.reference, '-I', self.__finalBam, '-O',  GATKdir +'/full_mutect.vcf',
//                                '--max-mnp-distance', '2'])
//             self.__CallCommand('LeftAlignAndTrimVariants', [self.__gatk, 'LeftAlignAndTrimVariants', 
//                                '-R', self.reference, '-V', GATKdir +'/mutect.vcf', '-O', GATKdir +'/gatk_mutect.vcf', '--split-multi-allelics'])
//             self.__CallCommand('mv', ['mv', GATKdir +'/gatk_mutect.vcf', GATKdir +'/mutect.vcf'])
//             self.__CallCommand('LeftAlignAndTrimVariants', [self.__gatk, 'LeftAlignAndTrimVariants', 
//                                '-R', self.reference, '-V', GATKdir +'/full_mutect.vcf', '-O', GATKdir +'/full_gatk_mutect.vcf', '--split-multi-allelics'])
//             self.__CallCommand('mv', ['mv', GATKdir +'/full_gatk_mutect.vcf', GATKdir +'/full_mutect.vcf'])
//             self.__CallCommand('FilterMutectCalls', [self.__gatk, 'FilterMutectCalls',
//                                '-R', self.reference, '-V', GATKdir +'/mutect.vcf', '--min-reads-per-strand', '1', '--min-median-read-position', '10', 
//                                '--min-allele-fraction', '0.01', '--microbial-mode', 'true', '-O', GATKdir + "/" + self.name + '_filter.vcf'])
//             self.__CallCommand('FilterMutectCalls', [self.__gatk, 'FilterMutectCalls',
//                                '-R', self.reference, '-V', GATKdir +'/full_mutect.vcf', '--min-reads-per-strand', '1', '--min-median-read-position', '10',
//                                '--min-allele-fraction', '0.01', '--microbial-mode', 'true', '-O', GATKdir + "/" + self.name + '_full_filter.vcf'])


//             """ Set final VCF file. """
            
//             if not self.__finalVCF: 
//                 self.__finalVCF = GATKdir + "/" + self.name + '_filter.vcf'
//             if not self.__fullVCF:
//                 self.__fullVCF = GATKdir + "/" + self.name + '_full_filter.vcf'
//         else:
//             # print error
//             pass  
       
//     def annotateVCF(self):
//         """ Annotate the final VCF file """
//         cwd = os.getcwd()
//         if self.__finalVCF:
//            self.__ifVerbose("Annotating final VCF.")
//            self.__CallCommand(['SnpEff', self.fOut + "/" + self.name +'_DR_loci_raw_annotation.txt'],
//                                 ['java', '-Xmx4g', '-jar', self.__annotator, 'NC_000962', self.__finalVCF])
//            self.__annotation = self.fOut + "/" + self.name +'_DR_loci_raw_annotation.txt'
//         if self.__fullVCF:
//            self.__ifVerbose("Annotating full VCF.")
//            self.__CallCommand(['SnpEff', self.fOut + "/" + self.name +'_full_raw_annotation.txt'],
//                                 ['java', '-Xmx4g', '-jar', self.__annotator, 'NC_000962', self.__fullVCF])
//            self.__full_annotation = self.fOut + "/" + self.name +'_full_raw_annotation.txt'
//            self.__ifVerbose("Parsing final Annotation.")
//            self.__CallCommand(['create annotation', self.fOut + "/" + self.name +'_DR_loci_annotation.txt'],
//                               ['python', self.__creater, self.__annotation, self.name])
//            self.__CallCommand(['create annotation', self.fOut + "/" + self.name +'_full_annotation.txt'],
//                               ['python', self.__creater, self.__full_annotation, self.name])
//            self.__CallCommand(['parse annotation', self.fOut + "/" + self.name +'_DR_loci_Final_annotation.txt'],
//                               ['python', self.__parser, self.__annotation, self.mutationloci, self.name])
//            self.__CallCommand(['parse annotation', self.fOut + "/" + self.name +'_full_Final_annotation.txt'],
//                               ['python', self.__parser, self.__full_annotation, self.mutationloci, self.name])
//         else:
//             self.__ifVerbose("Use SamTools, GATK, or Freebayes to annotate the final VCF.")
//         self.__CallCommand('rm', ['rm',  cwd + "/snpEff_genes.txt"])
//         self.__CallCommand('rm', ['rm',  cwd + "/snpEff_summary.html"])

//     def runLineage(self):
//         """ Run lineage Analysis """
//         self.__ifVerbose("Running Lineage Analysis")
//         self.__full_final_annotation = self.fOut + "/" + self.name +'_full_Final_annotation.txt'
//         self.__CallCommand(['lineage parsing', self.fOut + "/" + self.name +'_Lineage.txt'],
//                               ['python', self.__lineage_parser, self.__lineages, self.__full_final_annotation, self.__lineage, self.name])

    
 
//     def runPrint(self):
//         """ Print analysis report """
//         self.__ifVerbose("Printing report")           
//         self.__CallCommand(['create summary report', self.fOut + "/" + self.name + '_summary.txt'],
//                              ['python', self.__create_report, self.fOut + "/" + self.name + '_stats.txt',
//                               self.fOut + "/" + self.name + '_target_region_coverage.txt', self.fOut + "/" + self.name + '_DR_loci_Final_annotation.txt'])
//         self.__CallCommand(['run interpretation report', self.fOut + "/" + self.name + '_interpretation.txt'],
//                              ['python', self.__interpreter, self.__reported, self.fOut + "/" + self.name + '_summary.txt',
//                               self.fOut + "/" + self.name + '_structural_variants.txt', self.fOut + "/" + self.name + '_DR_loci_annotation.txt',
//                               self.fOut + "/" + self.name + '_target_region_coverage.txt', self.name]) 
//         self.__CallCommand('print pdf report',
//                              ['python', self.__print_report, self.fOut + "/" + self.name + '_summary.txt', self.fOut + "/" + self.name + '_report.pdf'])

//     def cleanUp(self):
//         """ Clean up the temporary files, and move them to a proper folder. """
//         i = datetime.now()
//         cwd = os.getcwd()
//         self.__CallCommand('rm', ['rm', '-r', self.outdir])
//         self.__CallCommand('rm', ['rm',  self.fOut + '/' + self.name + '_sdrcsm.bai'])
//         self.__CallCommand('rm', ['rm',  self.__finalBam])
//         self.__CallCommand('rm', ['rm', '-r', self.trimmomatic])
//         self.__CallCommand('rm', ['rm', '-r', self.clockwork])
//         self.__CallCommand('rm', ['rm', os.path.join(cwd, "config.yml")])
//         self.__CallCommand('rm', ['rm', '-r', os.path.join(cwd, "work")])
//         if os.path.isfile(self.fOut + "/" + self.name + '.log'):
//            self.__logFH.close()
//            self.__logged = False
//         fh4 = open(self.__log,'r')
//         for line in fh4:
//             lines = line.rstrip("\r\n")
//             if "Exception" in lines:
//                self.__exception = "positive"
//         fh4.close()        
//         if self.__exception == "positive":
//            self.__logFH2.write(i.strftime('%Y/%m/%d %H:%M:%S') + "\t" + "Input:" + "\t" + self.name + "\t" + "Exception in analysis\n")
//            self.__CallCommand('mv', ['mv', self.fOut, self.qlog])
        
//         if os.path.isfile(self.fOut + "/" + self.name + '_stats.txt'):
//            if os.path.getsize(self.fOut + "/" + self.name + '_stats.txt') < 1:
//               self.__logFH2.write(i.strftime('%Y/%m/%d %H:%M:%S') + "\t" + "Input:" + "\t" + self.name + "\t" + "Exception in analysis\n")
//               self.__CallCommand('mv', ['mv', self.fOut, self.qlog])        

//     def __ifVerbose(self, msg):
//         """ If verbose print a given message. """
//         if self.verbose: print msg



//   call fastqc.task_fastqc {
//     input:
//     forwardReads = task_repair.repaired_out1,
//     reverseReads = task_repair.repaired_out2
//   }

//   Boolean filter = task_fastqc.numberForwardReads > minNumberReads
//   if ( filter ) {
    
//     call trimmomatic.task_trimmomatic {
//       input:
//       read1 = task_repair.repaired_out1,
//       read2 = task_repair.repaired_out2,
//       samplename = samplename,
//       disk_size = disk_size_gb
//     }

//     if ( run_decontamination ) {
//       call bbduk.task_bbduk {
// 	input:
// 	read1_trimmed = task_trimmomatic.read1_trimmed,
// 	read2_trimmed = task_trimmomatic.read2_trimmed,
// 	samplename = samplename,
// 	disk_size = disk_size_gb
//       }
//     }

//     if ( run_clockwork_decontamination ) {
//       call cd.wf_clockwork_decontamination {
// 	input:
// 	reference = clockwork_contaminants,
// 	samplename = samplename,
//         metadata_file = clockwork_decontamination_metadata,
// 	input_reads_1 = select_first([task_bbduk.read1_clean, task_trimmomatic.read1_trimmed]),
// 	input_reads_2 = select_first([task_bbduk.read2_clean, task_trimmomatic.read2_trimmed]),
// 	disk_size = disk_size_gb
//       }
//     }

//     if ( run_ptrimmer ) {
//       call ptrimmer.task_ptrimmer {
// 	input:
// 	read1 = select_first([wf_clockwork_decontamination.clean_reads_1, task_bbduk.read1_clean, task_trimmomatic.read1_trimmed]),
// 	read2 = select_first([wf_clockwork_decontamination.clean_reads_2, task_bbduk.read2_clean, task_trimmomatic.read2_trimmed]),
// 	trim1 = "${samplename}_ptrimmer_1.fq.gz",
// 	trim2 = "${samplename}_ptrimmer_2.fq.gz",
// 	keep = ptrimmer_keep,
// 	seqtype = seqtype,
// 	ampfile = ampfile,
// 	summary = ptrimmer_summary,
// 	minqual = minqual,
// 	kmer = kmer,
// 	mismatch = mismatch
//       }
//     }
    
//     if ( run_fastqc_after_cleanup ) {
//       call fastqc.task_fastqc as task_fastqc_after_cleanup {
// 	input:
// 	forwardReads = select_first([task_ptrimmer.trimmedRead1, wf_clockwork_decontamination.clean_reads_1, task_bbduk.read1_clean, task_trimmomatic.read1_trimmed]),
// 	reverseReads = select_first([task_ptrimmer.trimmedRead2, wf_clockwork_decontamination.clean_reads_2, task_bbduk.read2_clean, task_trimmomatic.read2_trimmed]),
//       }
//     }
    
//     call varpipe.task_varpipe {
//       input:
//       read1 = select_first([wf_clockwork_decontamination.clean_reads_1, task_bbduk.read1_clean, task_trimmomatic.read1_trimmed]),
//       read2 = select_first([wf_clockwork_decontamination.clean_reads_2, task_bbduk.read2_clean, task_trimmomatic.read2_trimmed]),
//       reference = reference,
//       samplename = samplename,
//       config = config,
//       snpEff_config = snpEff_config,
//       snpEff_database = snpEff_data_dir,
//       outdir = outdir,
//       genome = genome,
//       keep = keep,
//       no_trim = no_trim,
//       whole_genome = whole_genome,
//       verbose = verbose,
//       disk_size = disk_size_gb
//     }

//     if ( run_delly ) {
//       call sv.wf_structural_variants {
// 	input:
// 	bam = task_varpipe.bam,
// 	bai = task_varpipe.bai,
// 	reference = reference,
// 	dataDir = snpEff_data_dir,
// 	config = snpEff_config,
// 	genome = genome
//       }

//       if (defined(wf_structural_variants.vcf_annotated)) {
// 	call concat.task_concat_2_vcfs {
// 	  input:
// 	  vcf1 = task_varpipe.DR_loci_raw_annotation,
// 	  vcf2 = select_first([wf_structural_variants.vcf_annotated]),
// 	  output_vcf_name = output_vcf_name
// 	}
//       }
//     }
    
//     if ( run_bamQC ) {
//       call bamQC.task_collect_multiple_metrics {
// 	input:
// 	bam = task_varpipe.bam,
// 	reference = reference
//       }
//       call wgsQC.task_collect_wgs_metrics {
// 	input:
// 	bam = task_varpipe.bam,
// 	reference = reference,
// 	bed = bed
//       }
//       call tpcrm.wf_collect_targeted_pcr_metrics {
// 	input:
// 	bam = task_varpipe.bam,
// 	reference = reference,
// 	amplicon_bed = bed,
// 	target_bed = bed
//       }
//       call doc.task_depth_of_coverage {
// 	input:
// 	bam = task_varpipe.bam,
// 	reference = reference,
// 	intervals = bed
//       }
//     }

//     call lineage.wf_lineage {
//       input:
//       vcf = task_varpipe.DR_loci_raw_annotation,
//       lineage_markers = lineage_markers,
//       samplename = samplename
//     }

//     if ( run_variant_interpretation ) {
//       call vi.wf_interpretation {
// 	input:
// 	vcf = select_first([task_concat_2_vcfs.concatenated_vcf, task_varpipe.DR_loci_raw_annotation]),
// 	bam = task_varpipe.bam,
// 	bai = task_varpipe.bai,
// 	bed = bed,
// 	json = json,
// 	samplename = samplename,
// 	lineage_information = wf_lineage.lineage_report
//       }
//     }
//   }