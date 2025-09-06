#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RECOMBINATION-AWARE ASSEMBLY SNPs WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Goal: Recombination-aware per-cluster trees for B. pseudomallei, 
          then a single overall tree via grafting.
    
    Production sequence:
    1. Cluster with Mash
    2. Per-cluster whole/core alignment (keep invariant A/T/C/G)
    3. Gubbins on the WGA
    4. Per-cluster final ML tree
    5. Select representative per cluster
    6. Backbone tree on representatives
    7. Graft subtrees onto backbone
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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
// MODULES: Local modules for input handling
//
include { INFILE_HANDLING_UNIX                             } from "../modules/local/infile_handling_unix/main"
include { INFILE_HANDLING_UNIX as REF_INFILE_HANDLING_UNIX } from "../modules/local/infile_handling_unix/main"

//
// MODULES: Step 1 - Clustering with Mash
//
include { MASH_SKETCH                                      } from "../modules/local/mash_sketch/main"
include { MASH_DIST                                        } from "../modules/local/mash_dist/main"
include { CLUSTER_GENOMES                                  } from "../modules/local/cluster_genomes/main"
include { MASH_TAB_TO_MATRIX                               } from "../modules/local/mash_tab_to_matrix/main"

//
// MODULES: Step 2 - Per-cluster whole/core alignment
//
include { SELECT_CLUSTER_REPRESENTATIVE                    } from "../modules/local/select_cluster_representative/main"
include { SNIPPY_ALIGN                                     } from "../modules/local/snippy_align/main"
include { CORE_GENOME_ALIGNMENT_PARSNP                     } from "../modules/local/core_genome_alignment_parsnp/main"
include { KEEP_INVARIANT_ATCG                              } from "../modules/local/keep_invariant_atcg/main"

//
// MODULES: Step 3 - Gubbins on WGA
//
//include { GUBBINS_CLUSTER                                  } from "../modules/local/gubbins_cluster/main"
include { GUBBINS_CLUSTER_DIAGNOSTIC } from '../modules/local/gubbins_cluster_diagnostic/main'
//
// MODULES: Step 4 - Per-cluster final ML tree
//
include { IQTREE_FAST                                      } from "../modules/local/iqtree_fast/main"
include { IQTREE_ASC                                       } from "../modules/local/iqtree_asc/main"

//
// MODULES: Step 5-7 - Representatives, backbone, and grafting
//
include { COLLECT_REPRESENTATIVES                          } from "../modules/local/collect_representatives/main"
include { BUILD_BACKBONE_TREE                              } from "../modules/local/build_backbone_tree/main"
include { GRAFT_SUBTREES                                   } from "../modules/local/graft_subtrees/main"

//
// SUBWORKFLOWS
//
include { INPUT_CHECK                                      } from "../subworkflows/local/input_check"
include { INPUT_CHECK as REF_INPUT_CHECK                   } from "../subworkflows/local/input_check"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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

workflow RECOMBINATION_AWARE_SNPS {

    // SETUP: Define empty channels to concatenate certain outputs
    ch_versions             = Channel.empty()
    ch_qc_filecheck         = Channel.empty()
    ch_output_summary_files = Channel.empty()

    /*
    ================================================================================
                            Preprocess input data
    ================================================================================
    */

    log.info "Starting recombination-aware SNP analysis workflow"
    log.info "Goal: Per-cluster recombination detection + tree grafting"

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

    // Create channel for assemblies: [ val(sample_id), path(assembly) ]
    ch_assemblies = INFILE_HANDLING_UNIX.out.input_files
        .flatten()
        .map { file -> 
            def sample_id = file.getBaseName().split('\\.')[0]
            tuple(sample_id, file)
        }

    // Handle reference genome (optional)
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
        ch_reference = Channel.empty()
    }

    /*
    ================================================================================
                        STEP 1: Cluster with Mash
    ================================================================================
    */

    log.info "STEP 1: Clustering genomes with Mash distances"

    // Create Mash sketches for each assembly
    MASH_SKETCH (
        ch_assemblies
    )
    ch_versions = ch_versions.mix(MASH_SKETCH.out.versions)

    // Calculate pairwise distances
    MASH_DIST (
        MASH_SKETCH.out.sketch.map{ sample_id, sketch -> sketch }.collect()
    )
    ch_versions = ch_versions.mix(MASH_DIST.out.versions)

    // Convert Mash tabular output to square matrix
    MASH_TAB_TO_MATRIX (
        MASH_DIST.out.distances
    )
    ch_versions = ch_versions.mix(MASH_TAB_TO_MATRIX.out.versions)

    // Cluster genomes based on matrix
    CLUSTER_GENOMES (
        MASH_TAB_TO_MATRIX.out.matrix
    )
    ch_versions = ch_versions.mix(CLUSTER_GENOMES.out.versions)

    // Create clustered assemblies channel: [ val(cluster_id), val(sample_ids), path(assemblies) ]
    ch_cluster_assignments = CLUSTER_GENOMES.out.clusters
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            def sample_id = row.sample_id.split('\\.')[0]
            tuple(row.cluster_id, sample_id)
        }

    ch_clustered_assemblies = ch_cluster_assignments
        .map { cluster_id, sample_id -> tuple(sample_id, cluster_id) }
        .join(ch_assemblies, by: 0)
        .map { sample_id, cluster_id, assembly -> tuple(cluster_id, sample_id, assembly) }
        .groupTuple(by: 0)
        .filter { cluster_id, sample_ids, assemblies -> sample_ids.size() >= 3 }

    /*
    ================================================================================
                    STEP 2: Per-cluster whole/core alignment
    ================================================================================
    */

    log.info "STEP 2: Creating per-cluster whole genome alignments"

    // Select cluster representative (medoid/high-quality)
    SELECT_CLUSTER_REPRESENTATIVE (
        ch_clustered_assemblies,
        MASH_TAB_TO_MATRIX.out.matrix
    )
    ch_versions = ch_versions.mix(SELECT_CLUSTER_REPRESENTATIVE.out.versions)

    // Create channel for alignment: [ cluster_id, sample_ids, assemblies, representative_id, reference ]
    ch_for_alignment = ch_clustered_assemblies
        .join(SELECT_CLUSTER_REPRESENTATIVE.out.representative.map { cluster_id, rep_id, rep_file -> 
            tuple(cluster_id, rep_id, rep_file) 
        }, by: 0)
        .map { cluster_id, sample_ids, assemblies, rep_id, rep_file ->
            tuple(cluster_id, sample_ids, assemblies, rep_id, rep_file)
        }

    // Choose alignment method: Snippy (preferred) or Parsnp (fallback)
    if (params.alignment_method == 'snippy' || !params.alignment_method) {
        log.info "Using Snippy for per-cluster whole genome alignment"
        
        SNIPPY_ALIGN (
            ch_for_alignment
        )
        ch_versions = ch_versions.mix(SNIPPY_ALIGN.out.versions)
        ch_core_alignments = SNIPPY_ALIGN.out.core_alignment
        
    } else if (params.alignment_method == 'parsnp') {
        log.info "Using Parsnp for per-cluster core genome alignment"
        
        // Prepare for Parsnp: [ meta, assemblies ], [ meta, reference ]
        ch_parsnp_input = ch_for_alignment
            .map { cluster_id, sample_ids, assemblies, rep_id, rep_file ->
                def meta = [snp_package: cluster_id]
                tuple(meta, assemblies)
            }
        
        ch_parsnp_ref = ch_for_alignment
            .map { cluster_id, sample_ids, assemblies, rep_id, rep_file ->
                def meta = [snp_package: cluster_id]
                tuple(meta, rep_file)
            }
        
        CORE_GENOME_ALIGNMENT_PARSNP (
            ch_parsnp_input,
            ch_parsnp_ref
        )
        ch_versions = ch_versions.mix(CORE_GENOME_ALIGNMENT_PARSNP.out.versions)
        
        // Convert Parsnp output to expected format
        ch_core_alignments = CORE_GENOME_ALIGNMENT_PARSNP.out.output
            .map { meta, files ->
                def cluster_id = meta.snp_package
                def alignment_file = files.find { it.name.endsWith('.fa.gz') || it.name.endsWith('.fa') }
                tuple(cluster_id, alignment_file)
            }
    } else {
        error "Unknown alignment method: ${params.alignment_method}. Use 'snippy' or 'parsnp'"
    }

    // Keep invariant A/T/C/G sites (guardrail: do not feed SNP-only to Gubbins)
    KEEP_INVARIANT_ATCG (
        ch_core_alignments
    )
    ch_versions = ch_versions.mix(KEEP_INVARIANT_ATCG.out.versions)

    /*
    ================================================================================
                        STEP 3: Gubbins on the WGA
    ================================================================================
    */

    log.info "STEP 3: Running Gubbins for recombination detection"

    // Create starting trees for Gubbins using IQ-TREE fast mode
    ch_starting_trees = KEEP_INVARIANT_ATCG.out.core_alignment
        .map { cluster_id, alignment ->
            tuple(cluster_id, alignment)
        }

    // Use IQTREE_FAST to create starting trees
    IQTREE_FAST (
        ch_starting_trees
    )
    ch_versions = ch_versions.mix(IQTREE_FAST.out.versions)

    // Prepare input for Gubbins: join alignments with starting trees
    ch_for_gubbins = KEEP_INVARIANT_ATCG.out.core_alignment
        .join(IQTREE_FAST.out.tree, by: 0)

    GUBBINS_CLUSTER_DIAGNOSTIC(
        ch_for_gubbins.map { cluster_id, alignment, tree -> cluster_id },
        ch_for_gubbins.map { cluster_id, alignment, tree -> alignment },
        ch_for_gubbins.map { cluster_id, alignment, tree -> tree }
    )
    ch_versions = ch_versions.mix(GUBBINS_CLUSTER_DIAGNOSTIC.out.versions)

    /*
    ================================================================================
                    STEP 4: Per-cluster final ML tree
    ================================================================================
    */

    log.info "STEP 4: Building per-cluster final ML trees with ASC correction"

    // Combine Gubbins filtered SNPs with representative info for IQ-TREE
    ch_for_final_tree = GUBBINS_CLUSTER_DIAGNOSTIC.out.filtered_alignment
        .join(SELECT_CLUSTER_REPRESENTATIVE.out.representative.map { cluster_id, rep_id, rep_file ->
            tuple(cluster_id, rep_id)
        }, by: 0)

    IQTREE_ASC (
        ch_for_final_tree
    )
    ch_versions = ch_versions.mix(IQTREE_ASC.out.versions)

    /*
    ================================================================================
                    STEP 5: Select representatives per cluster
    ================================================================================
    */

    log.info "STEP 5: Collecting cluster representatives"

    // Collect all representative files
    COLLECT_REPRESENTATIVES (
        SELECT_CLUSTER_REPRESENTATIVE.out.representative.map { cluster_id, rep_id, rep_file -> rep_file }.collect(),
        CLUSTER_GENOMES.out.clusters
    )
    ch_versions = ch_versions.mix(COLLECT_REPRESENTATIVES.out.versions)

    /*
    ================================================================================
                    STEP 6: Backbone tree on representatives
    ================================================================================
    */

    log.info "STEP 6: Building backbone tree from representatives"

    BUILD_BACKBONE_TREE (
        COLLECT_REPRESENTATIVES.out.representatives_fasta
    )
    ch_versions = ch_versions.mix(BUILD_BACKBONE_TREE.out.versions)

    /*
    ================================================================================
                    STEP 7: Graft subtrees onto backbone
    ================================================================================
    */

    log.info "STEP 7: Grafting cluster subtrees onto backbone tree"

    GRAFT_SUBTREES (
        BUILD_BACKBONE_TREE.out.backbone_tree,
        IQTREE_ASC.out.final_tree.map { cluster_id, tree, rep_id -> tree }.collect(),
        COLLECT_REPRESENTATIVES.out.representatives_mapping
    )
    ch_versions = ch_versions.mix(GRAFT_SUBTREES.out.versions)

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

    /*
    ================================================================================
                        Log completion
    ================================================================================
    */

    GRAFT_SUBTREES.out.grafted_tree
        .subscribe { tree ->
            log.info "✅ WORKFLOW COMPLETE!"
            log.info "Final grafted tree: ${tree}"
            log.info "Recombination-aware per-cluster analysis with tree grafting finished successfully"
        }
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