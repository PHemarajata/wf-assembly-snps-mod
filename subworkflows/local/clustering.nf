//
// Pre-clustering using Mash distances
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Local modules
//
include { MASH_SKETCH      } from "../../modules/local/mash_sketch/main"
include { MASH_DIST        } from "../../modules/local/mash_dist/main"
include { CLUSTER_GENOMES  } from "../../modules/local/cluster_genomes/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN CLUSTERING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CLUSTERING {

    take:
    ch_assemblies // channel: [ val(sample_id), path(assembly) ]

    main:
    ch_versions = Channel.empty()

    // PROCESS: Create Mash sketches for each assembly
    MASH_SKETCH (
        ch_assemblies
    )
    ch_versions = ch_versions.mix(MASH_SKETCH.out.versions)

    // PROCESS: Calculate pairwise distances
    MASH_DIST (
        MASH_SKETCH.out.sketch.map{ sample_id, sketch -> sketch }.collect()
    )
    ch_versions = ch_versions.mix(MASH_DIST.out.versions)

    // PROCESS: Cluster genomes based on distances
    CLUSTER_GENOMES (
        MASH_DIST.out.distances
    )
    ch_versions = ch_versions.mix(CLUSTER_GENOMES.out.versions)

    // Create channel with cluster assignments
    ch_cluster_assignments = CLUSTER_GENOMES.out.clusters
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            // Extract sample ID without extension to match ch_assemblies format
            // Use the same logic as in the scalable workflow
            def sample_id = row.sample_id.split('\\.')[0]
            tuple(row.cluster_id, sample_id)
        }

    // Join with original assemblies to create clustered channel
    ch_clustered_raw = ch_cluster_assignments
        .map { cluster_id, sample_id -> tuple(sample_id, cluster_id) }
        .join(ch_assemblies, by: 0)
        .map { sample_id, cluster_id, assembly -> tuple(cluster_id, sample_id, assembly) }
        .groupTuple(by: 0)
        .branch { cluster_id, sample_ids, assemblies ->
            phylo_viable: sample_ids.size() >= 3  // At least 3 samples for phylogenetic analysis
            small_cluster: sample_ids.size() == 2  // Exactly 2 samples - too small for trees
            singleton: sample_ids.size() == 1      // Single sample clusters
        }

    // Log small clusters that will be skipped
    ch_clustered_raw.small_cluster
        .subscribe { cluster_id, sample_ids, assemblies ->
            log.info "Skipping cluster ${cluster_id} with ${sample_ids.size()} samples (${sample_ids.join(', ')}) - phylogenetic analysis requires at least 3 samples"
        }

    ch_clustered_raw.singleton
        .subscribe { cluster_id, sample_ids, assemblies ->
            log.info "Skipping singleton cluster ${cluster_id} (sample: ${sample_ids[0]}) - phylogenetic analysis requires at least 3 samples"
        }

    // Handle small clusters and singletons based on parameter
    if (params.merge_singletons) {
        // Merge all small clusters (singletons + 2-sample clusters) into one large cluster
        ch_small_clusters_combined = ch_clustered_raw.singleton
            .mix(ch_clustered_raw.small_cluster)
            .map { cluster_id, sample_ids, assemblies -> 
                // Flatten the lists in case of 2-sample clusters
                [sample_ids, assemblies].transpose()
            }
            .collect()
            .map { all_pairs ->
                if (all_pairs.size() >= 3) {
                    def merged_sample_ids = all_pairs.collect { it[0] }
                    def merged_assemblies = all_pairs.collect { it[1] }
                    tuple("merged_small_clusters", merged_sample_ids, merged_assemblies)
                } else {
                    // If fewer than 3 total samples, skip
                    null
                }
            }
            .filter { it != null }

        // Combine phylogenetically viable clusters with merged small clusters
        ch_clustered_final = ch_clustered_raw.phylo_viable
            .mix(ch_small_clusters_combined)
    } else {
        // Process only phylogenetically viable clusters (3+ samples)
        ch_clustered_final = ch_clustered_raw.phylo_viable
    }

    // Transform to final format - keep sample_ids and assemblies separate for SKA_BUILD
    ch_clustered_assemblies = ch_clustered_final
        .map { cluster_id, sample_ids, assemblies ->
            tuple(cluster_id, sample_ids, assemblies)
        }

    // Count clusters for summary
    ch_clustered_assemblies
        .count()
        .subscribe { count ->
            if (count == 0) {
                log.warn "No multi-sample clusters found for phylogenetic analysis."
                log.warn "This could be due to:"
                log.warn "  - All samples being singletons (each in its own cluster)"
                log.warn "  - Sample ID mismatch between clustering and assembly files"
                log.warn "Consider:"
                log.warn "  - Adjusting --mash_threshold (current: ${params.mash_threshold ?: 0.03}) to allow more similar samples to cluster together"
                log.warn "  - Using --merge_singletons to combine all singletons into one large cluster"
                log.warn "  - Checking that sample names in clusters.tsv match assembly file names"
                log.warn "Phylogenetic analysis requires at least 3 samples per cluster."
            } else {
                log.info "Found ${count} clusters for phylogenetic analysis"
            }
        }

    emit:
    versions             = ch_versions
    clusters             = CLUSTER_GENOMES.out.clusters
    cluster_summary      = CLUSTER_GENOMES.out.summary
    clustered_assemblies = ch_clustered_assemblies // channel: [ val(cluster_id), val(sample_ids), path(assemblies) ]
}