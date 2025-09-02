nextflow.enable.dsl=2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSNPS.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet or directory not specified!' }
if (params.ref) { ch_ref_input = file(params.ref) } else { ch_ref_input = [] }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Local modules
//
include { INFILE_HANDLING_UNIX                             } from "../modules/local/infile_handling_unix/main"
include { INFILE_HANDLING_UNIX as REF_INFILE_HANDLING_UNIX } from "../modules/local/infile_handling_unix/main"

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                                      } from "../subworkflows/local/input_check"
include { INPUT_CHECK as REF_INPUT_CHECK                   } from "../subworkflows/local/input_check"
include { CLUSTERING                                       } from "../subworkflows/local/clustering"
include { FILTER_ATCG_COLUMNS                              } from "../modules/local/filter_atcg_columns/main.nf"
include { GUBBINS_CLUSTER_DIAGNOSTIC                       } from "../modules/local/gubbins_cluster_diagnostic/main.nf"
include { INTEGRATE_RESULTS                                } from "../subworkflows/local/integrate_results"
include { USHER_PLACEMENT                                  } from "../subworkflows/local/usher_placement"
include { SKA_ALIGN                                        } from "../modules/local/ska_align"
include { IQTREE_FAST                                      } from "../modules/local/iqtree_fast"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CREATE CHANNELS FOR INPUT PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Workflow mode selection
def workflow_mode = params.workflow_mode ?: 'cluster'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Convert input to lowercase
def toLower(it) {
    it.toString().toLowerCase()
}

// Check QC filechecks for a failure
def qcfilecheck(process, qcfile, inputfile) {
    qcfile.map{ meta, file -> [ meta, [file] ] }
            .join(inputfile)
            .map{ meta, qc, input ->
                data = []
                qc.flatten().each{ data += it.readLines() }

                if ( data.any{ it.contains('FAIL') } ) {
                    line = data.last().split('\t')
                    if (line.first() != "NaN") {
                        log.warn("${line[1]} QC check failed during process ${process} for sample ${line.first()}")
                    } else {
                        log.warn("${line[1]} QC check failed during process ${process}")
                    }
                } else {
                    [ meta, input ]
                }
            }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASSEMBLY_SNPS_SCALABLE {

    // SETUP: Define empty channels to concatenate certain outputs
    ch_versions             = Channel.empty()
    ch_qc_filecheck         = Channel.empty()
    ch_output_summary_files = Channel.empty()

    /*
    ================================================================================
                            Preprocess input data
    ================================================================================
    */

    // SUBWORKFLOW: Check input for samplesheet or pull inputs from directory
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Check input files meet size criteria
    INFILE_HANDLING_UNIX (
        INPUT_CHECK.out.input_files
    )
    ch_versions     = ch_versions.mix(INFILE_HANDLING_UNIX.out.versions)
    ch_qc_filecheck = ch_qc_filecheck.concat(INFILE_HANDLING_UNIX.out.qc_filecheck)

    // Create channel for assemblies
    ch_assemblies = INFILE_HANDLING_UNIX.out.input_files
        .flatten()
        .map { file -> 
            def sample_id = file.getBaseName().split('\\.')[0]
            tuple(sample_id, file)
        }

    // Handle reference genome
    if (params.ref) {
        REF_INPUT_CHECK (
            ch_ref_input
        )
        ch_versions = ch_versions.mix(REF_INPUT_CHECK.out.versions)

        REF_INFILE_HANDLING_UNIX (
            REF_INPUT_CHECK.out.input_files
        )
        ch_versions        = ch_versions.mix(REF_INFILE_HANDLING_UNIX.out.versions)
        ch_qc_filecheck    = ch_qc_filecheck.concat(REF_INFILE_HANDLING_UNIX.out.qc_filecheck)

        ch_reference = REF_INFILE_HANDLING_UNIX.out.input_files
            .flatten()
            .first()
    } else {
        // Use largest assembly as reference
        ch_reference = ch_assemblies
            .map { sample_id, assembly -> assembly }
            .toSortedList { a, b -> b.size() <=> a.size() }
            .map { it[0] }
    }

    /*
    ================================================================================
                            Workflow mode selection
    ================================================================================
    */

    if (workflow_mode == 'cluster') {
        /*
        ================================================================================
                                Clustering mode (default)
        ================================================================================
        */

        log.info "Running in clustering mode for scalable analysis"

        // SUBWORKFLOW: Pre-cluster genomes using Mash distances
        CLUSTERING (
            ch_assemblies
        )
        ch_versions = ch_versions.mix(CLUSTERING.out.versions)

        // Invoke SKA_ALIGN for each cluster
        SKA_ALIGN (
            CLUSTERING.out.clustered_assemblies.map { cluster_id, sample_ids, assemblies -> tuple(cluster_id, assemblies) }
        )
        ch_versions = ch_versions.mix(SKA_ALIGN.out.versions)

        // Invoke IQTREE_FAST for each alignment
        IQTREE_FAST (
            SKA_ALIGN.out.alignment
        )
        ch_versions = ch_versions.mix(IQTREE_FAST.out.versions)

        // Invoke FILTER_ATCG_COLUMNS for each alignment
        FILTER_ATCG_COLUMNS(
            SKA_ALIGN.out.alignment
        )
        ch_versions = ch_versions.mix(FILTER_ATCG_COLUMNS.out.versions)

        // Join filtered alignments with trees by cluster_id
        ch_for_gubbins = FILTER_ATCG_COLUMNS.out.filtered_alignment
            .join(IQTREE_FAST.out.tree, by: 0)   // -> (cluster_id, filtered.fa, tree)

        // Use multiMap to separate inputs for GUBBINS_CLUSTER_DIAGNOSTIC
        ch_for_gubbins
            .multiMap { cluster_id, alignment, tree ->
                cluster_ids: cluster_id
                alignments: alignment
                trees: tree
            }
            .set { gubbins_input }

        GUBBINS_CLUSTER_DIAGNOSTIC(
            gubbins_input.cluster_ids,
            gubbins_input.alignments,
            gubbins_input.trees
        )
        ch_versions = ch_versions.mix(GUBBINS_CLUSTER_DIAGNOSTIC.out.versions)

        // Collect cluster results for integration
        ch_cluster_alignments = GUBBINS_CLUSTER_DIAGNOSTIC.out.filtered_alignment
        ch_cluster_trees     = GUBBINS_CLUSTER_DIAGNOSTIC.out.final_tree

        // SUBWORKFLOW: Integrate results across clusters (if enabled)
        if (params.integrate_results) {
            log.info "Integrating core SNPs and grafting trees from all clusters"
            
            INTEGRATE_RESULTS (
                ch_cluster_alignments,
                ch_cluster_trees,
                CLUSTERING.out.clusters,
                CLUSTERING.out.cluster_summary
            )
            ch_versions = ch_versions.mix(INTEGRATE_RESULTS.out.versions)
        } else {
            log.info "Skipping result integration (disabled by --integrate_results false)"
        }

        // Log integration results (if integration was performed)
        if (params.integrate_results) {
            INTEGRATE_RESULTS.out.integrated_alignment
                .subscribe { alignment ->
                    log.info "Created integrated core SNPs alignment: ${alignment}"
                }

            INTEGRATE_RESULTS.out.supertree
                .subscribe { supertree ->
                    log.info "Created grafted supertree: ${supertree}"
                }

            INTEGRATE_RESULTS.out.integrated_tree
                .subscribe { tree ->
                    log.info "Built phylogenetic tree from integrated core SNPs: ${tree}"
                }
        }

        // Optional: Build UShER MAT for incremental updates
        if (params.build_usher_mat && params.integrate_results) {
            ch_global_alignment = INTEGRATE_RESULTS.out.integrated_alignment
                .map { alignment ->
                    def meta = [id: 'integrated_global']
                    tuple(meta, alignment)
                }

            USHER_PLACEMENT (
                ch_global_alignment,
                ch_reference,
                Channel.empty() // No new samples for initial build
            )
            ch_versions = ch_versions.mix(USHER_PLACEMENT.out.versions)
        }

    } else if (workflow_mode == 'place') {
        /*
        ================================================================================
                                Placement mode (incremental updates)
        ================================================================================
        */

        log.info "Running in placement mode for incremental updates"

        if (!params.existing_mat) {
            error "Placement mode requires --existing_mat parameter"
        }

        // Load existing MAT and place new samples
        ch_existing_mat = Channel.fromPath(params.existing_mat)
            .map { mat -> 
                def meta = [id: 'existing']
                tuple(meta, mat)
            }

        // Convert new assemblies to VCF format (simplified - would need proper implementation)
        // This would require aligning new samples to reference and calling SNPs
        log.warn "Placement mode implementation requires additional SNP calling steps"

    } else {
        /*
        ================================================================================
                                Global mode (original approach)
        ================================================================================
        */

        log.info "Running in global mode (original Parsnp-based approach)"
        log.warn "Global mode may not scale well beyond 200 samples"

        // Fall back to original workflow logic
        // This would include the original Parsnp-based approach
        // (Implementation would mirror the original assembly_snps.nf)
    }

    /*
    ================================================================================
                        Collect QC information
    ================================================================================
    */

    // Collect QC file check information
    ch_qc_filecheck = ch_qc_filecheck
                        .map{ meta, file -> file }
                        .collectFile(
                            name:       "Summary.QC_File_Checks.tsv",
                            keepHeader: true,
                            storeDir:   "${params.outdir}/Summaries",
                            sort:       'index'
                        )

    ch_output_summary_files = ch_output_summary_files.mix(ch_qc_filecheck.collect())

    /*
    ================================================================================
                        Collect version information
    ================================================================================
    */

    // Collect version information
    ch_versions
        .unique()
        .collectFile(
            name:     "software_versions.yml",
            storeDir: params.tracedir
        )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/