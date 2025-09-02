//
// UShER-based phylogenetic placement for incremental updates
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Local modules
//
include { SNP_SITES     } from "../../modules/local/snp_sites/main"
include { USHER_BUILD   } from "../../modules/local/usher_build/main"
include { USHER_PLACE   } from "../../modules/local/usher_place/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN USHER PLACEMENT WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow USHER_PLACEMENT {

    take:
    ch_global_alignment // channel: [ val(meta), path(alignment) ]
    ch_reference        // channel: path(reference)
    ch_new_samples_vcf  // channel: path(new_samples.vcf) - optional for placement

    main:
    ch_versions = Channel.empty()

    // PROCESS: Convert alignment to VCF format
    SNP_SITES (
        ch_global_alignment
    )
    ch_versions = ch_versions.mix(SNP_SITES.out.versions)

    // PROCESS: Build initial MAT from VCF
    USHER_BUILD (
        SNP_SITES.out.vcf,
        ch_reference
    )
    ch_versions = ch_versions.mix(USHER_BUILD.out.versions)

    // PROCESS: Place new samples if provided
    if (ch_new_samples_vcf) {
        USHER_PLACE (
            USHER_BUILD.out.mat,
            ch_new_samples_vcf
        )
        ch_versions = ch_versions.mix(USHER_PLACE.out.versions)
        
        ch_final_mat = USHER_PLACE.out.updated_mat
        ch_placement_stats = USHER_PLACE.out.placement_stats
    } else {
        ch_final_mat = USHER_BUILD.out.mat
        ch_placement_stats = Channel.empty()
    }

    emit:
    versions         = ch_versions
    vcf              = SNP_SITES.out.vcf
    mat              = ch_final_mat
    placement_stats  = ch_placement_stats
}