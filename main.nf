#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/tbAnalyzer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/tbAnalyzer
    Website: https://nf-co.re/tbAnalyzer
    Slack  : https://nfcore.slack.com/channels/tbAnalyzer
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Main workflows
include { tbAnalyzer              } from './workflows/tbanalyzer'
include { testPrep                } from './workflows/testprep'

// Subworkflows
include { CREATE_INPUT_CHANNEL    } from './subworkflows/local/create_input_channel'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_tbAnalyzer_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_tbAnalyzer_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_tbAnalyzer_pipeline'

// Modules
include { PREPSRR                 } from './modules/local/srr'
include { BASESPACE                 } from './modules/local/basespace'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
// params.fasta = getGenomeAttribute('fasta')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Workflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// WORKFLOW: PREPARE TEST DATA
workflow OhioTestPrep {
    // set test samplesheet
    samplesheet = file(params.input)
    if(params.isTest==false) {exit 1, "YEP"}
    main:
        // read input
        CREATE_INPUT_CHANNEL(
            samplesheet
        )

        // Download test data
        BASESPACE(CREATE_INPUT_CHANNEL.out.reads)
        ch_reads = BASESPACE.out.reads
        
        // RUN ANALYZER
        tbAnalyzer (
            ch_reads
        )
}

// WORKFLOW: PREPARE TEST DATA
workflow OhioSRRPrep {
    // set test samplesheet
    samplesheet = file(params.input)
    if(params.isTest==false) {exit 1, "YEP"}
    main:
        // read input
        CREATE_INPUT_CHANNEL(
            samplesheet
        )

        // Download test data
        PREPSRR(CREATE_INPUT_CHANNEL.out.reads)
        ch_reads = PREPSRR.out.reads
        
        // RUN ANALYZER
        tbAnalyzer (
            ch_reads
        )
}

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow OhioTBAnalyzer {
    if (params.input) { ch_input = file(params.input) } else { exit 1, 'For -entry NFCORE_OhioTBGenomics: Input samplesheet not specified!' }
    samplesheet = file(params.input)

    main:

    // // Initialize
    //  PIPELINE_INITIALISATION (
    //     params.version,
    //     params.help,
    //     params.validate_params,
    //     params.monochrome_logs,
    //     args,
    //     params.outdir,
    //     params.input
    // )

    // read input
    CREATE_INPUT_CHANNEL(
        samplesheet
    )
    ch_reads=CREATE_INPUT_CHANNEL.out.reads

    // Download test data
    BASESPACE(CREATE_INPUT_CHANNEL.out.reads)
    ch_reads = BASESPACE.out.reads
        
    // RUN ANALYZER
    tbAnalyzer (
        ch_reads
    )

    // PIPELINE_COMPLETION (
    //     params.email,
    //     params.email_on_fail,
    //     params.plaintext_email,
    //     params.outdir,
    //     params.monochrome_logs,
    //     params.hook_url,
    //     tbAnalyzer.out.multiqc_report
    // )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

