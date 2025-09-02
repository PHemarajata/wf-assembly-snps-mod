//
// Integrate results from cluster analysis
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Local modules
//
include { EXTRACT_CORE_SNPS      } from "../../modules/local/extract_core_snps/main"
include { INTEGRATE_CORE_SNPS    } from "../../modules/local/integrate_core_snps/main"
include { GRAFT_TREES            } from "../../modules/local/graft_trees/main"
include { BUILD_INTEGRATED_TREE  } from "../../modules/local/build_integrated_tree/main"
include { CREATE_FINAL_SUMMARY   } from "../../modules/local/create_final_summary/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN INTEGRATION WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow INTEGRATE_RESULTS {

    take:
    ch_cluster_alignments // channel: [ val(cluster_id), path(alignment) ]
    ch_cluster_trees      // channel: [ val(cluster_id), path(tree) ]
    ch_clusters_file      // channel: path(clusters.tsv)
    ch_cluster_summary    // channel: path(cluster_summary.txt) - optional

    main:
    ch_versions = Channel.empty()

    // PROCESS: Extract core SNPs from each cluster alignment
    EXTRACT_CORE_SNPS (
        ch_cluster_alignments,
        ch_clusters_file
    )
    ch_versions = ch_versions.mix(EXTRACT_CORE_SNPS.out.versions)

    // PROCESS: Integrate core SNPs across all clusters
    INTEGRATE_CORE_SNPS (
        EXTRACT_CORE_SNPS.out.core_snps.map { cluster_id, snps -> snps }.collect(),
        EXTRACT_CORE_SNPS.out.snp_positions.map { cluster_id, positions -> positions }.collect(),
        ch_clusters_file
    )
    ch_versions = ch_versions.mix(INTEGRATE_CORE_SNPS.out.versions)

    // PROCESS: Graft cluster trees into supertree
    GRAFT_TREES (
        ch_cluster_trees.map { cluster_id, tree -> tree }.collect(),
        ch_clusters_file,
        INTEGRATE_CORE_SNPS.out.integrated_alignment
    )
    ch_versions = ch_versions.mix(GRAFT_TREES.out.versions)

    // PROCESS: Build phylogenetic tree from integrated core SNPs
    BUILD_INTEGRATED_TREE (
        INTEGRATE_CORE_SNPS.out.integrated_alignment,
        INTEGRATE_CORE_SNPS.out.sample_mapping
    )
    ch_versions = ch_versions.mix(BUILD_INTEGRATED_TREE.out.versions)

    // PROCESS: Create comprehensive final summary
    CREATE_FINAL_SUMMARY (
        ch_clusters_file,
        ch_cluster_summary.ifEmpty(Channel.empty()),
        INTEGRATE_CORE_SNPS.out.summary,
        GRAFT_TREES.out.report,
        BUILD_INTEGRATED_TREE.out.report,
        INTEGRATE_CORE_SNPS.out.sample_mapping,
        GRAFT_TREES.out.tree_summary
    )
    ch_versions = ch_versions.mix(CREATE_FINAL_SUMMARY.out.versions)

    emit:
    versions                = ch_versions
    integrated_alignment    = INTEGRATE_CORE_SNPS.out.integrated_alignment
    integrated_positions    = INTEGRATE_CORE_SNPS.out.integrated_positions
    core_snp_summary        = INTEGRATE_CORE_SNPS.out.summary
    sample_mapping          = INTEGRATE_CORE_SNPS.out.sample_mapping
    supertree              = GRAFT_TREES.out.supertree
    tree_grafting_report   = GRAFT_TREES.out.report
    tree_summary           = GRAFT_TREES.out.tree_summary
    integrated_tree        = BUILD_INTEGRATED_TREE.out.tree
    integrated_tree_log    = BUILD_INTEGRATED_TREE.out.log
    phylogeny_report       = BUILD_INTEGRATED_TREE.out.report
    final_html_report      = CREATE_FINAL_SUMMARY.out.html_report
    final_text_report      = CREATE_FINAL_SUMMARY.out.text_report
    analysis_statistics    = CREATE_FINAL_SUMMARY.out.statistics
}