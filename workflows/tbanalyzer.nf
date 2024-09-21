/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                            } from '../modules/local/fastqc'
include { CLOCKWORK                         } from '../modules/local/clockwork'
include { TRIMM                             } from '../modules/local/trim'
include { FASTQ_SCREEN                      } from '../modules/local/fastq_screen'
include { SAMTOOLS_FLAGSTAT                 } from '../modules/local/samtools_flag'
include { BWA_INDEX                         } from '../modules/local/bwa_index'
include { BWA_MEM                           } from '../modules/local/bwa_mem'


// include { PICARD_CREATESEQUENCEDICTIONARY   } from '../modules/local/picard_seq_dict'
// include { BBDUK                  } from '../modules/local/bbduk'
// include { CENTRIFUGE_DB          } from '../modules/local/centrifuge_db'
// include { CENTRIFUGE             } from '../modules/local/centrifuge'
// include { TB_PROFILER            } from '../modules/local/tb_profiler'
// include { BCF2VCF                } from '../modules/local/bcf2vcf'
// include { MULTIQC                } from '../modules/nf-core/multiqc/main'
// include { paramsSummaryMap       } from 'plugin/nf-validation'


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

    // MODULE: CLOCKWORK
    CLOCKWORK(
        ch_reads,
        params.clockContam,
        params.clockMetadata
    )
    ch_bam=CLOCKWORK.out.deconBam

    // MODULE: Repair
    TRIMM (
        CLOCKWORK.out.reads
    )
    ch_trimmed=TRIMM.out.reads

    // MODULE: FASTQ Screen
    FASTQ_SCREEN (
        ch_trimmed,
        params.fastq_screen_configuration,
        params.contam_dir
    )

    // MODULE: SAMTOOLS SORT, FLAGSTAT
    SAMTOOLS_FLAGSTAT (
        ch_bam
    )

    // Index
    ch_fasta=params.fasta_ref
    BWA_INDEX(
        ch_fasta
    )
    ch_index=BWA_INDEX.out.index
    
    // Align
    BWA_MEM(
        ch_trimmed,
        ch_index,
        ch_fasta,
        params.sort_bam
    )


    // // PICARD_CREATESEQUENCEDICTIONARY
    // PICARD_CREATESEQUENCEDICTIONARY(
    //     ch_trimmed
    // )

    // GATK

    // SAMTOOLS

    // ANNOTATE

    // LINEAGE

    
    // //
    // // MODULE: Repair
    // BBDUK (
    //     ch_trimmed,
    //     params.contam_dir
    // )

    // // MODULE: CENTRIFUGE
    // CENTRIFUGE(
    //     ch_trimmed,
    //     params.hpv_tar
    // )

    // // MODULE: TB PROFILER
    // TB_PROFILER(
    //     BBDUK.out.repaired_reads,
    //     params.report,
    //     params.mapper,
    //     params.caller,
    //     params.min_depth,
    //     params.min_af,
    //     params.min_af_pred,
    //     params.cov_frac_threshold,
    //     params.threads
    // )

    // // MODULE: BCF2VCF
    // BCF2VCF (
    //     TB_PROFILER.out.bcf
    // )

    // MODULE: SNPEFF
    // SNPEFF (

    // )

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

// version 1.0

// DONE
// import "./task_concatenate_fastq.wdl" as concatenate_fastq DONE
// import "./task_fastqc.wdl" as fastqc DONE
// import "./task_fastq_screen.wdl" as fastq_screen DONE
// import "./task_trimmomatic.wdl" as trimmomatic DONE
// import "./task_bbduk.wdl" as bbduk
// import "./wf_clockwork_decontamination.wdl" as cd
// import "./task_tbprofiler.wdl" as tbprofiler
// import "./task_bcf2vcf.wdl" as bcf2vcf
// import "./task_snpEff.wdl" as snpEff


// TODO
// import "./task_concat_2_vcfs.wdl" as concat

// import "./task_collect_multiple_metrics.wdl" as bamQC
// import "./task_collect_wgs_metrics.wdl" as wgsQC
// import "./wf_collect_targeted_pcr_metrics.wdl" as tpcrm
// import "./task_depth_of_coverage.wdl" as doc

// import "./wf_lineage.wdl" as lineage
// import "./wf_interpretation.wdl" as vi

// import "./task_multiqc.wdl" as multiQC

    
//     Boolean run_bamQC = true

//     # snpEff
//     File snpEff_data_dir
//     File snpEff_config
//     String genome = "Mycobacterium_tuberculosis_h37rv"
//     String annotated_structural_variants_name = "annotated_structural_variants.vcf"
    
//     # concat vcfs
//     String output_vcf_name = "concatenated.vcf"

//     # variant interpretation
//     File bed
//     File json
//     File lineage_markers
//   }



    
//     call tbprofiler.task_tbprofiler {
//       input:
//       read1 = select_first([task_bbduk.read1_clean, task_trimmomatic.read1_trimmed]),
//       read2 = select_first([task_bbduk.read2_clean, task_trimmomatic.read1_trimmed]),
//       samplename = samplename
//     }

//     # structural variants
//     call bcf2vcf.task_bcf2vcf {
//       input:
//       bcf_file = task_tbprofiler.bcf
//     }
    
//     call snpEff.task_snpEff {
//       input:
//       vcf = task_bcf2vcf.vcf_file,
//       genome = genome,
//       config = snpEff_config,
//       dataDir = snpEff_data_dir,
//       outputPath = annotated_structural_variants_name
//     }

//     call concat.task_concat_2_vcfs {
//       input:
//       vcf1 = task_tbprofiler.vcf,
//       vcf2 = task_snpEff.outputVcf,
//       output_vcf_name = output_vcf_name
//     }

//     # bam QC
//     if ( run_bamQC ) {
//       call bamQC.task_collect_multiple_metrics {
// 	input:
// 	bam = task_tbprofiler.bam,
// 	reference = reference
//       }
//       call wgsQC.task_collect_wgs_metrics {
// 	input:
// 	bam = task_tbprofiler.bam,
// 	reference = reference,
// 	bed = bed
//       }
//       call tpcrm.wf_collect_targeted_pcr_metrics {
// 	input:
// 	bam = task_tbprofiler.bam,
// 	reference = reference,
// 	amplicon_bed = bed,
// 	target_bed = bed
//       }
//       call doc.task_depth_of_coverage {
// 	input:
// 	bam = task_tbprofiler.bam,
// 	reference = reference,
// 	intervals = bed
//       }
//     }

//     call lineage.wf_lineage {
//       input:
//       vcf = task_tbprofiler.vcf,
//       lineage_markers = lineage_markers,
//       samplename = samplename
//     }

//     call vi.wf_interpretation {
//       input:
//       vcf = select_first([task_concat_2_vcfs.concatenated_vcf, task_tbprofiler.vcf]),
//       bam = task_tbprofiler.bam,
//       bai = task_tbprofiler.bai,
//       bed = bed,
//       json = json,
//       samplename = samplename,
//       lineage_information = wf_lineage.lineage_report
//     }

//   }
//   # end filter
  
//   Array[File] allReports = flatten([
//   select_all([task_trimmomatic.trim_err,
//   task_fastq_screen.txt,
//   task_fastqc.forwardData,
//   task_fastqc.reverseData,
//   task_bbduk.adapter_stats,
//   task_bbduk.phiX_stats,
//   task_collect_wgs_metrics.collectMetricsOutput,
//   wf_collect_targeted_pcr_metrics.output_metrics ]),
//   flatten(select_all([task_collect_multiple_metrics.collectMetricsOutput,[]]))
//   ])
//   if ( length(allReports) > 0 ) {
//     call multiQC.task_multiqc {
//       input:
//       inputFiles = allReports,
//       outputPrefix = samplename
//     }
//   }
  
//   output {
//     File? csv = task_tbprofiler.csv
//     File? bam = task_tbprofiler.bam
//     File? bai = task_tbprofiler.bai
//     File? vcf = task_tbprofiler.vcf
//     File? svs = task_bcf2vcf.vcf_file
//     # output from trimmer trimmomatic
//     File? trim_stats = task_trimmomatic.trim_stats
//     # output from bbduk decontamination
//     File? phiX_stats = task_bbduk.phiX_stats
//     File? adapter_stats = task_bbduk.adapter_stats
//     File? polyA_stats   = task_bbduk.polyA_stats
//     File? Ecoli_stats   = task_bbduk.Ecoli_stats
//     File? Covid19_stats = task_bbduk.Covid19_stats
//     # output from clockwork decontamination
//     File? clockwork_decontamination_stats = wf_clockwork_decontamination.stats
//     # output from fastqc
//     File forwardHtml = task_fastqc.forwardHtml
//     File reverseHtml = task_fastqc.reverseHtml
//     File forwardZip = task_fastqc.forwardZip
//     File reverseZip = task_fastqc.reverseZip
//     File forwardSummary = task_fastqc.forwardSummary
//     File reverseSummary = task_fastqc.reverseSummary
//     File forwardData = task_fastqc.forwardData
//     File reverseData = task_fastqc.reverseData
//     # output from fastq_screen
//     File? fastq_screen_html = task_fastq_screen.html
//     File? fastq_screen_txt = task_fastq_screen.txt
//     File? fastq_screen_tagged = task_fastq_screen.tagged
//     File? fastq_screen_tagged_filter = task_fastq_screen.tagged_filter
//     # output from bam QC
//     Array[File]? multiple_metrics_outputs = task_collect_multiple_metrics.collectMetricsOutput
//     Array[File]? depth_of_coverage_outputs = task_depth_of_coverage.outputs
//     File? collect_wgs_output_metrics = task_collect_wgs_metrics.collectMetricsOutput
//     File? collect_targeted_pcr_metrics = wf_collect_targeted_pcr_metrics.output_metrics
//     # all annotated variants = variant valler + SV caller delly
//     File? concatenated_vcf = task_concat_2_vcfs.concatenated_vcf
//     # lineage
//     File? lineage_report = wf_lineage.lineage_report
//     # variant interpretation
//     File? lab_log = wf_interpretation.lab_log
//     File? lab_report = wf_interpretation.lab_report
//     File? lims_report = wf_interpretation.lims_report
//     # multiqc
//     File? multiqc_report = task_multiqc.report
//   }

//   meta {
//     author: "Dieter Best"
//     email: "Dieter.Best@cdph.ca.gov"
//     description: "## variant pipeline \n This is the London TB profiler: https://github.com/jodyphelan/TBProfiler.\n\n This also runs fastq QC, decontamination, and alignment QC."
//   }
  
//   parameter_meta {
//     read1: {
//       description: "List of fastq files with forward reads.",
//       category: "required"
//     }
//     read2: {
//       description: "List of fastq files with reverse reads.",
//       category: "required"
//     }
//     reference: {
//       description: "Reference sequence to align to.",
//       category: "required"
//     }
//     samplename: {
//       description: "Name of the sample.",
//       category: "required"
//     }
//     bed: {
//       description: "bed file with genomic intervals of interest. Note: reference name in case of London TB profiler is 'Chromosome', make sure to use correct bed file",
//       category: "required"
//     }
//     json: {
//       description: "json file with drug information for variants.",
//       category: "required"
//     }
//     report: {
//       description: "Name for output tsv file.",
//       category: "optional"
//     }
//     run_decontamination: {
//       description: "Flag, turn on if decontamination of fastq files should be run.",
//       category: "optional"
//     }
//     run_bamQC: {
//       description: "Flag for performing alignment bam QC.",
//       category: "optional"
//     }
//     # output
//     adapter_stats: {description: "Name file where decontamination procedure writes adapter contamination statistics to."}
//     phiX_stats: {description: "phiX contamination report from bbduk decontamination task."}
//     polyA_stats: {description: "polyA contamination report from bbduk decontamination task."}
//     Ecoli_stats: {description: "Ecoli contamination report from bbduk decontamination task."}
//     Covid19_stats: {description: "Covid19 contamination report from bbduk decontamination task."}

//     bam: {description: "Output alignement file of alignment procedure, aligner is bwa."}
//     bai: {description: "Index file for output alignement file of alignment procedure, aligner is bwa."}
//     csv: {description: "Ouput variant call file in csv format."}
//     vcf: {description: "Ouput variant call file in vcf format."}
//     svs: {description: "Ouput structural variants call file in vcf format."}

//     collectMetricsOutput: {description: "Array of output files from alignment bam QC."}

//     forwardData: {description: "Fastqc output data for forward reads."}
//     forwardHtml: {description: "Fastqc output html file for forward reads."}
//     forwardSummary: {description: "Fastqc output summary file for forward reads."}
//     forwardZip: {description: "Fastqc output zip file for forward reads."}
    
//     reverseData: {description: "Fastqc output data for reverse reads."}
//     reverseHtml: {description: "Fastqc output html file for reverse reads."}
//     reverseSummary: {description: "Fastqc output summary file for reverse reads."}
//     reverseZip: {description: "Fastqc output zip file for reverse reads."}

//     trim_stats: {description: "Output text file for read trimming statistics."}
//     interpretation_report: {description: "Output tsv file from variant interpretation."}
//     multiqc_report: {description: "Output html file with QC summary report."}
//   }

// }