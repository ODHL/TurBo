/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include { PREP_DB               } from '../../modules/local/clockwork_prepDB.nf'
// include { BUILD_SIF             } from '../../modules/local/clockwork_build.nf'
include { CLOCKWORK_REFS                  } from '../../modules/local/clockwork_refs.nf'
include { CLOCKWORK_MAP as MAP_RAW        } from '../../modules/local/clockwork_map.nf'
include { CLOCKWORK_MAP as MAP_CLEAN      } from '../../modules/local/clockwork_map.nf'
include { CLOCKWORK_REMOVE                } from '../../modules/local/clockwork_remove.nf'
include { CLOCKWORK_PREP                  } from '../../modules/local/clockwork_prep.nf'
include { CLOCKWORK_VARIANTS              } from '../../modules/local/clockwork_variants.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CLOCKWORK {

    take:
        ch_reads // channel: 

    main:

        // prep empty channels
        ch_versions = Channel.empty()

        // Set references
        // CLOCKWORK_REFS (
        //     params.clockwork_contaminants,
        //     params.clockwork_ntm
        // )
        // ch_metadata=CLOCKWORK_REFS.out.metadata

        // map reads
        MAP_RAW(
            ch_reads,
            params.clockContam
        )
        ch_unsortedSam=MAP_RAW.out.deconSam

        // decontaminate reads
        CLOCKWORK_REMOVE(
            ch_unsortedSam,
            params.clockMetadata
        )
        ch_deconFq=CLOCKWORK_REMOVE.out.deconFq

        // map cleaned reads
        MAP_CLEAN(
            ch_deconFq,
            params.clockContam
        )
        ch_unsortedCleanedSam=MAP_CLEAN.out.deconSam

        // // Index reference
        // CLOCKWORK_PREP(
        //     params.fasta_ref
        // )

        // // Call Variants
        // CLOCKWORK_VARIANTS(
        //     ch_deconFq,
        //     CLOCKWORK_PREP.out.indexedRefs
        // )

    emit:
        ch_unsortedSam        = ch_unsortedSam
        ch_deconFq            = ch_deconFq
        ch_unsortedCleanedSam = ch_unsortedCleanedSam
}
