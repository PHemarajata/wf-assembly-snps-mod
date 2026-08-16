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
include { MASH_SKETCH_BATCH                                } from "../modules/local/mash_sketch_batch/main"
include { MASH_PASTE                                       } from "../modules/local/mash_paste/main"
include { MASH_TRIANGLE                                    } from "../modules/local/mash_triangle/main"
include { CLUSTER_GENOMES                                  } from "../modules/local/cluster_genomes/main"
// Alternative clustering front end: PopPUNK bgmm+refine, for collections Mash
// single-linkage cannot partition (--clustering_method poppunk).
include { POPPUNK_CLUSTER                                  } from "../modules/local/poppunk_cluster/main"

//
// MODULES: Step 2 - Per-cluster whole/core alignment
//
include { SELECT_CLUSTER_REPRESENTATIVE                    } from "../modules/local/select_cluster_representative/main"
// Picks a COMPLETE per-cluster mapping reference (medoid != reference; see the
// module header). Enabled with --pick_complete_references.
include { PICK_CLUSTER_REFERENCES                          } from "../modules/local/pick_cluster_references/main"
include { SPLIT_REFERENCE_REPLICONS                        } from "../modules/local/split_reference_replicons/main"
// Snippy is now scatter (one task per cluster-sample) + snippy-core gather. The old
// monolithic SNIPPY_ALIGN ran every sample serially inside a single task.
include { SNIPPY_SCATTER                                   } from "../modules/local/snippy_align/main"
include { SNIPPY_CORE_GATHER                               } from "../modules/local/snippy_align/main"
// SKA2 reference-anchored path (alignment_method='ska'): fast low-spec alternative.
// NOTE: `ska map`, not `ska align`. ska's own help calls `align` an "unordered
// alignment" and `map` an "ordered alignment using a reference sequence"; Gubbins
// scans for recombination spatially along the genome, so unordered columns are
// invalid input. Measured here with ska 0.5.1 on 30 draft assemblies: map 0.49 s ->
// 376,564 columns, align 0.29 s -> 284,078 columns (both ~99.9% constant, so column
// COUNT is not the discriminator -- ordering is).
include { SKA_BUILD_SAMPLE                                 } from "../modules/local/ska_map_align/main"
include { SKA_MAP_ALIGN                                    } from "../modules/local/ska_map_align/main"
include { CORE_GENOME_ALIGNMENT_PARSNP                     } from "../modules/local/core_genome_alignment_parsnp/main"
include { KEEP_INVARIANT_ATCG                              } from "../modules/local/keep_invariant_atcg/main"

//
// MODULES: Step 3 - Gubbins on WGA
//
include { GUBBINS_CLUSTER                                  } from "../modules/local/gubbins_cluster/main"
//
// MODULES: Step 4 - Per-cluster final ML tree
//
include { IQTREE_FAST                                      } from "../modules/local/iqtree_fast/main"
include { ASC_PREFLIGHT                                    } from "../modules/local/asc_preflight/main"
include { IQTREE_ASC                                       } from "../modules/local/iqtree_asc/main"

//
// MODULES: Step 5-7 - Representatives, backbone, and grafting
//
include { COLLECT_REPRESENTATIVES                          } from "../modules/local/collect_representatives/main"
include { BUILD_BACKBONE_TREE                              } from "../modules/local/build_backbone_tree/main"
include { SELECT_UNIT_MEDOID                               } from "../modules/local/select_unit_medoid/main"
include { GLOBAL_CORE_ALIGNMENT                            } from "../modules/local/global_core_alignment/main"
include { GLOBAL_ML_TREE                                   } from "../modules/local/global_ml_tree/main"
include { SUMMARIZE_CLUSTER_PHYLOGENY                      } from "../modules/local/summarize_cluster_phylogeny/main"
include { GRAFT_TREES                                      } from "../modules/local/graft_trees/main"

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

    // Handle reference genome (optional). When provided, ch_reference carries the
    // path; when absent, ch_reference emits a NO_FILE sentinel so downstream
    // processes can branch on file existence without requiring optional inputs.
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
        ch_reference = Channel.value(file("${projectDir}/assets/NO_FILE"))
    }

    /*
    ================================================================================
        STEP 1 (CURATED MODE): supplied partition + per-cluster references
    ================================================================================
    When --cluster_assignments AND --cluster_references are BOTH given, Mash
    clustering and medoid selection are skipped entirely; the analysis is driven
    by the externally-computed partition and the reference decisions that go with
    it, exactly as reference_sensitivity_bp.py does. This is the faithful-
    reproduction path for a collection that Mash single-linkage cannot partition:
    measured, all 2802 B. pseudomallei collapse into ONE component at every
    threshold >= 0.007, and shatter into hundreds of singletons below that, so a
    PopPUNK/fastbaps partition supplied here is the only valid clustering.

    clusters.tsv     : header `cluster_id<TAB>sample_id`, one row per genome.
    references.tsv   : header `cluster_id<TAB>reference_path`, one row per cluster;
                       the reference is external (complete or borrowed), NOT a
                       cluster member, matching the manual --reference contract.
    */
    def curated_mode = (params.cluster_assignments && params.cluster_references)

    if (curated_mode) {

        log.info "CURATED MODE: partition=${params.cluster_assignments}, references=${params.cluster_references}. Mash clustering and medoid selection are SKIPPED."

        ch_curated_assignments = Channel.fromPath(params.cluster_assignments, checkIfExists: true)
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.sample_id.toString().split('\\.')[0], row.cluster_id.toString()) }

        ch_curated_refs = Channel.fromPath(params.cluster_references, checkIfExists: true)
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.cluster_id.toString(), file(row.reference_path, checkIfExists: true)) }

        ch_clustered_assemblies = ch_curated_assignments
            .join(ch_assemblies, by: 0)
            .map { sample_id, cluster_id, assembly -> tuple(cluster_id, sample_id, assembly) }
            .groupTuple(by: 0)

        // Same <3-taxa accounting as the default path: a cluster of <3 cannot
        // yield a tree. Report rather than silently drop.
        ch_clustered_assemblies
            .filter { cluster_id, sample_ids, assemblies -> sample_ids.size() < 3 }
            .map { cluster_id, sample_ids, assemblies -> sample_ids.size() }
            .sum()
            .subscribe { dropped ->
                if (dropped > 0) {
                    log.warn "EXCLUDED_FROM_ANALYSIS: ${dropped} genome(s) are in supplied clusters of <3 taxa and cannot produce a tree."
                }
            }

        ch_clustered_assemblies = ch_clustered_assemblies
            .filter { cluster_id, sample_ids, assemblies -> sample_ids.size() >= 3 }

        // ch_for_alignment shape: [cluster_id, sample_ids, assemblies, rep_id, reference].
        // rep_id := cluster_id: the reference is external, so there is no medoid
        // cluster member to name here.
        ch_for_alignment = ch_clustered_assemblies
            .join(ch_curated_refs, by: 0)
            .map { cluster_id, sample_ids, assemblies, reference ->
                tuple(cluster_id, sample_ids, assemblies, cluster_id, reference)
            }

        // (cluster_id, rep_id) label channel, consumed once in STEP 4. Re-read
        // the file rather than ch_curated_refs (a queue channel already consumed
        // above). rep_id := cluster_id since the reference is external.
        ch_cluster_repid = Channel.fromPath(params.cluster_references, checkIfExists: true)
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.cluster_id.toString(), row.cluster_id.toString()) }

        // The supplied partition file IS the clusters.tsv the per-cluster summary
        // needs (same header cluster_id<TAB>sample_id).
        ch_clusters_file = Channel.fromPath(params.cluster_assignments, checkIfExists: true)

    } else {

    /*
    ================================================================================
                        STEP 1: Cluster with Mash
    ================================================================================
    */

    log.info "STEP 1: Clustering genomes with Mash distances"

    // Sketch in BATCHES rather than one task per genome. The original scattered
    // MASH_SKETCH over every assembly, so 2000 genomes meant 2000 tasks (and 2000
    // container starts) to do work that `mash sketch -p` parallelizes internally.
    // params.mash_batch_size genomes per task; ~10 tasks at n=2000, batch=200.
    ch_sketch_batches = ch_assemblies
        .collate( params.mash_batch_size ?: 200 )
        .toList()
        .flatMap { batches ->
            batches.withIndex().collect { batch, i ->
                tuple( i, batch.collect { it[0] }, batch.collect { it[1] } )
            }
        }

    MASH_SKETCH_BATCH (
        ch_sketch_batches
    )
    ch_versions = ch_versions.mix(MASH_SKETCH_BATCH.out.versions.first())

    // Merge the per-batch sketches into one .msh
    MASH_PASTE (
        MASH_SKETCH_BATCH.out.sketch.collect()
    )
    ch_versions = ch_versions.mix(MASH_PASTE.out.versions)

    // `mash triangle -p` computes each pair once and writes lower-triangular
    // Phylip directly, replacing MASH_DIST + MASH_TAB_TO_MATRIX (the old path
    // never passed -p despite cpus=4, and built an n^2-row TSV intermediate).
    MASH_TRIANGLE (
        MASH_PASTE.out.sketch
    )
    ch_versions = ch_versions.mix(MASH_TRIANGLE.out.versions)

    // Cluster genomes. Two front ends, same `cluster_id<TAB>sample_id` contract.
    //
    // Mash single-linkage CANNOT partition every collection. Measured on the
    // 2,802-genome B. pseudomallei set: one connected component at threshold
    // >= 0.007, 786/1,094 singletons below it, nothing workable in between -- so
    // what comes out is a size-capped chop of one component, an imposed
    // partition rather than a found one. `--clustering_method poppunk` runs
    // PopPUNK's bgmm+refine fit instead, which gives 271 strains on the same
    // input. Mash still runs either way: MASH_TRIANGLE's matrix is what
    // SELECT_CLUSTER_REPRESENTATIVE uses to pick a medoid, and it is cheap next
    // to PopPUNK sketching.
    def clustering_method = (params.clustering_method ?: 'mash').toString().toLowerCase()
    if (!(clustering_method in ['mash', 'poppunk'])) {
        error "clustering_method must be 'mash' or 'poppunk', got '${clustering_method}'"
    }

    if (clustering_method == 'poppunk') {
        log.info "STEP 1: Clustering with PopPUNK (bgmm + refine)"
        // ch_assemblies carries [ sample_id, assembly_path ]; PopPUNK builds its
        // own rfile from the staged files, so pass the paths only.
        POPPUNK_CLUSTER (
            ch_assemblies.map { it[1] }.collect()
        )
        ch_versions = ch_versions.mix(POPPUNK_CLUSTER.out.versions)
        ch_clusters_tsv = POPPUNK_CLUSTER.out.clusters
        // PopPUNK emits no per-cluster submatrices; the ch_rep_input mapping
        // below already falls back to the full MASH_TRIANGLE matrix when a
        // cluster has none, so nothing else needs to change.
        ch_submatrices = Channel.empty()
    } else {
        CLUSTER_GENOMES (
            MASH_TRIANGLE.out.matrix
        )
        ch_versions = ch_versions.mix(CLUSTER_GENOMES.out.versions)
        ch_clusters_tsv = CLUSTER_GENOMES.out.clusters
        ch_submatrices  = CLUSTER_GENOMES.out.submatrices
    }

    // Create clustered assemblies channel: [ val(cluster_id), val(sample_ids), path(assemblies) ]
    ch_cluster_assignments = ch_clusters_tsv
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            def sample_id = row.sample_id.split('\\.')[0]
            tuple(row.cluster_id, sample_id)
        }

    ch_grouped_clusters = ch_cluster_assignments
        .map { cluster_id, sample_id -> tuple(sample_id, cluster_id) }
        .join(ch_assemblies, by: 0)
        .map { sample_id, cluster_id, assembly -> tuple(cluster_id, sample_id, assembly) }
        .groupTuple(by: 0)

    // Clusters with <3 taxa cannot yield a tree. Previously they were dropped by a
    // bare .filter with no accounting, so genomes could vanish from the analysis
    // silently -- at a tight --mash_threshold that is a large fraction of the input
    // (measured: 37 of 112 real genomes at 0.002 vs 16 at 0.003). Report it.
    ch_grouped_clusters
        .filter { cluster_id, sample_ids, assemblies -> sample_ids.size() < 3 }
        .map { cluster_id, sample_ids, assemblies -> sample_ids.size() }
        .sum()
        .subscribe { dropped ->
            if (dropped > 0) {
                log.warn "EXCLUDED_FROM_ANALYSIS: ${dropped} genome(s) are in clusters of <3 taxa " +
                         "and cannot produce a tree (mash_threshold=${params.mash_threshold}, " +
                         "merge_singletons=${params.merge_singletons})."
                if (!params.merge_singletons) {
                    log.warn "  These genomes are DROPPED entirely. Raise --mash_threshold so they " +
                             "join real clusters, or set --merge_singletons true to pool them."
                } else {
                    log.warn "  --merge_singletons=true pools them into one 'merged_small_clusters' " +
                             "bin regardless of similarity; on real data that bin spanned nearly the " +
                             "full diversity of the collection, so treat its recombination calls with " +
                             "caution (see AUDIT_REPORT.md section G.3)."
                }
            }
        }

    ch_clustered_assemblies = ch_grouped_clusters
        .filter { cluster_id, sample_ids, assemblies -> sample_ids.size() >= 3 }

    /*
    ================================================================================
                    STEP 2: Per-cluster whole/core alignment
    ================================================================================
    */

    log.info "STEP 2: Creating per-cluster whole genome alignments"

    // Select cluster representative (medoid/high-quality)
    // The module now takes ONE channel whose 4th element is the cluster's own
    // submatrix, written by CLUSTER_GENOMES --emit-submatrices. The original
    // passed the full n x n matrix as a second channel and re-parsed it in every
    // per-cluster task. Clusters with no submatrix (e.g. merged_small_clusters)
    // fall back to the full matrix, which the script reads with pandas usecols.
    // cluster_mash.py writes cluster_matrices/cluster_<id>.matrix.tsv, and
    // Nextflow's simpleName strips at the first dot, giving exactly the
    // cluster_id used in clusters.tsv (verified against real output).
    ch_submatrix_by_cluster = ch_submatrices
        .flatten()
        .map { f -> tuple(f.simpleName, f) }

    ch_rep_input = ch_clustered_assemblies
        .join( ch_submatrix_by_cluster, by: 0, remainder: true )
        .combine( MASH_TRIANGLE.out.matrix )
        .map { row ->
            // join(remainder:true) yields [cluster_id, sample_ids, assemblies, submatrix]
            // with submatrix == null when the cluster has none (e.g. merged_small_clusters);
            // combine() appends the full matrix as the last element.
            def cluster_id  = row[0]
            def sample_ids  = row[1]
            def assemblies  = row[2]
            def submatrix   = row[3]
            def full_matrix = row[-1]
            tuple(cluster_id, sample_ids, assemblies, submatrix ?: full_matrix)
        }
        .filter { cluster_id, sample_ids, assemblies, matrix ->
            // remainder:true also emits submatrices with no matching cluster; drop those
            sample_ids != null && assemblies != null
        }

    SELECT_CLUSTER_REPRESENTATIVE (
        ch_rep_input
    )
    ch_versions = ch_versions.mix(SELECT_CLUSTER_REPRESENTATIVE.out.versions)

    // Create channel for alignment: [ cluster_id, sample_ids, assemblies, representative_id, reference ]
    // If params.use_global_reference is enabled and params.ref was supplied, swap
    // the per-cluster medoid for the global reference so every cluster aligns
    // against the same anchor (e.g., K96243). Default keeps medoid behavior.
    // A complete per-cluster reference, chosen by completeness-gate +
    // centrality, replacing the medoid for MAPPING only. The medoid still
    // serves as the cluster's backbone representative -- they are different
    // objects (see PICK_CLUSTER_REFERENCES). Required for --split_replicons on
    // a draft-heavy collection, where the medoid is usually multi-contig.
    if (params.pick_complete_references) {
        // Staged as a file, not interpolated as a path, so its contents join the
        // task hash: adding a newly-identified bad reference must invalidate the
        // cached selection rather than silently reuse a pick made without it.
        ch_ref_blocklist = params.reference_blocklist
            ? Channel.fromPath(params.reference_blocklist, checkIfExists: true)
            : Channel.fromPath("${projectDir}/assets/NO_BLOCKLIST")

        PICK_CLUSTER_REFERENCES (
            ch_clusters_tsv,
            MASH_TRIANGLE.out.matrix,
            ch_assemblies.map { it[1] }.collect(),
            ch_ref_blocklist
        )
        ch_versions = ch_versions.mix(PICK_CLUSTER_REFERENCES.out.versions)

        ch_picked_ref = PICK_CLUSTER_REFERENCES.out.references
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.cluster_id, file(row.reference_path, checkIfExists: true)) }

        ch_for_alignment = ch_clustered_assemblies
            .join(SELECT_CLUSTER_REPRESENTATIVE.out.representative.map { cluster_id, rep_id, rep_file ->
                tuple(cluster_id, rep_id)
            }, by: 0)
            .join(ch_picked_ref, by: 0)
            .map { cluster_id, sample_ids, assemblies, rep_id, picked_ref ->
                tuple(cluster_id, sample_ids, assemblies, rep_id, picked_ref)
            }
    } else if (params.use_global_reference && params.ref) {
        log.info "use_global_reference=true: overriding per-cluster medoid with params.ref for alignment"
        ch_for_alignment = ch_clustered_assemblies
            .join(SELECT_CLUSTER_REPRESENTATIVE.out.representative.map { cluster_id, rep_id, rep_file ->
                tuple(cluster_id, rep_id)
            }, by: 0)
            .combine(ch_reference)
            .map { cluster_id, sample_ids, assemblies, rep_id, global_ref ->
                tuple(cluster_id, sample_ids, assemblies, rep_id, global_ref)
            }
    } else {
        ch_for_alignment = ch_clustered_assemblies
            .join(SELECT_CLUSTER_REPRESENTATIVE.out.representative.map { cluster_id, rep_id, rep_file ->
                tuple(cluster_id, rep_id, rep_file)
            }, by: 0)
            .map { cluster_id, sample_ids, assemblies, rep_id, rep_file ->
                tuple(cluster_id, sample_ids, assemblies, rep_id, rep_file)
            }
    }

    // (cluster_id, rep_id) label channel, mirrored by the curated branch so
    // STEP 4 works in both modes. Process outputs are re-readable, so this does
    // not disturb the reads at ch_for_alignment above.
    ch_cluster_repid = SELECT_CLUSTER_REPRESENTATIVE.out.representative
        .map { cluster_id, rep_id, rep_file -> tuple(cluster_id, rep_id) }

    ch_clusters_file = ch_clusters_tsv

    }  // end of default (Mash-clustering) vs curated-mode branch; ch_for_alignment is set either way

    /*
    ================================================================================
        STEP 1b (optional): split each cluster's reference into replicons
    ================================================================================
    When --split_replicons is set, the per-cluster reference is split into one
    FASTA per replicon (contig) and each cluster fans out into one analysis unit
    per replicon, keyed cluster_id__<replicon>. Everything downstream (alignment,
    Gubbins, tree) then runs per replicon, matching the manual pipeline
    (close__ska_map__chr1 / chr2). This is why Gubbins never sees a contig
    junction. Requires a COMPLETE reference (<= max_replicons contigs), so pair it
    with curated mode / supplied complete references; the split module fails loudly
    on a draft.
    */
    if (params.split_replicons) {
        log.info "STEP 1b: splitting per-cluster references into replicons (max_replicons=${params.max_replicons ?: 4}); analysis fans out per replicon."

        // Fork so the reference can be split while the per-cluster metadata is
        // preserved for the rejoin (ch_for_alignment is a queue channel and
        // cannot be read twice).
        ch_for_alignment.multiMap { cid, sample_ids, assemblies, rep_id, ref ->
            meta: tuple(cid, sample_ids, assemblies, rep_id)
            ref:  tuple(cid, ref)
        }.set { ch_fa_forked }

        SPLIT_REFERENCE_REPLICONS ( ch_fa_forked.ref )
        ch_versions = ch_versions.mix(SPLIT_REFERENCE_REPLICONS.out.versions.first())

        ch_replicons = SPLIT_REFERENCE_REPLICONS.out.replicons
            .transpose()   // one emission per replicon fasta

        // Cross-join metadata (1 per cluster) with replicons (N per cluster) by
        // cluster_id; re-key each unit as cluster_id__<replicon>.
        ch_split_units = ch_fa_forked.meta
            .combine(ch_replicons, by: 0)
            .map { cid, sample_ids, assemblies, rep_id, replicon_fa ->
                def compound = "${cid}__${replicon_fa.simpleName}"
                tuple(compound, sample_ids, assemblies, rep_id, replicon_fa)
            }

        // Tee the fanned-out units into the alignment channel and the STEP 4
        // label channel. multiMap forks one source into two so neither read
        // splits the other's emissions.
        ch_split_units.multiMap { cid, sample_ids, assemblies, rep_id, ref ->
            aln:   tuple(cid, sample_ids, assemblies, rep_id, ref)
            label: tuple(cid, rep_id)
        }.set { ch_split_forked }

        ch_for_alignment = ch_split_forked.aln
        ch_cluster_repid = ch_split_forked.label
    }

    // Fork off the per-unit assemblies for STEP 4b's medoid selection BEFORE any
    // alignment consumer reads ch_for_alignment.
    //
    // ch_for_alignment is a queue channel: reading it twice SPLITS its emissions
    // between the readers rather than duplicating them, which is why every other
    // multi-consumer point in this workflow forks with multiMap. Reading it again
    // further down (where STEP 4b lives) silently stole emissions from the
    // alignment path -- measured as KEEP_INVARIANT_ATCG re-executing all 164
    // units on a resume that should have been fully cached.
    if (params.global_ml_tree == null ? true : params.global_ml_tree) {
        ch_for_alignment.multiMap { cid, sample_ids, assemblies, rep_id, ref ->
            aln:    tuple(cid, sample_ids, assemblies, rep_id, ref)
            medoid: tuple(cid, sample_ids, assemblies)
        }.set { ch_fa_gml }
        ch_for_alignment   = ch_fa_gml.aln
        ch_unit_assemblies = ch_fa_gml.medoid
    } else {
        ch_unit_assemblies = Channel.empty()
    }

    // Choose alignment method: Snippy (scatter/gather), SKA2 (fast low-spec), or Parsnp
    if (params.alignment_method == 'snippy' || !params.alignment_method) {
        log.info "Using Snippy (scattered per sample) for per-cluster whole genome alignment"

        // SCATTER: explode [cluster, [sample_ids], [assemblies], rep_id, ref] into one
        // emission per (cluster, sample). sample_id and assembly are transposed
        // TOGETHER, so the sample->file binding is carried by the channel. This
        // replaces the old in-task `[[ "$file" == *"$sample"* ]]` substring search,
        // which mis-bound any sample id that was a substring of another.
        ch_snippy_scatter = ch_for_alignment
            .flatMap { cluster_id, sample_ids, assemblies, rep_id, ref ->
                def ids   = sample_ids  instanceof List ? sample_ids  : [sample_ids]
                def files = assemblies  instanceof List ? assemblies  : [assemblies]
                assert ids.size() == files.size() :
                    "cluster ${cluster_id}: ${ids.size()} sample ids but ${files.size()} assemblies"
                // Sorting by sample id makes snippy-core's column order deterministic
                // across resumes, which matters because Gubbins output is order-sensitive.
                def pairs = [ids, files].transpose().sort { it[0] }
                pairs
                    .findAll { sid, asm -> sid != rep_id }   // reference is not aligned to itself
                    .collect { sid, asm -> tuple(cluster_id, sid, asm, rep_id, ref) }
            }

        SNIPPY_SCATTER (
            ch_snippy_scatter
        )
        ch_versions = ch_versions.mix(SNIPPY_SCATTER.out.versions.first())

        // GATHER: regroup per-sample snippy dirs by cluster, rejoin rep_id/reference.
        //
        // The per-cluster size is carried by groupKey, NOT just asserted afterwards.
        // A bare groupTuple(by: 0) cannot know a group is finished until the whole
        // upstream channel closes, so EVERY cluster's gather waits for the LAST
        // snippy task in the entire run. Measured on the 82-unit L1 run: 132 of 164
        // replicon-units were fully mapped with zero gathers submitted, which
        // serialises Gubbins behind all 4,140 mappings instead of letting it start
        // on the units that are already done.
        //
        // groupKey(key, size) tells Nextflow the expected group size up front, so a
        // group emits the moment its last member arrives. The assert is kept: it now
        // guards the groupKey size rather than substituting for it.
        ch_cluster_sizes = ch_snippy_scatter
            .map { cluster_id, sid, asm, rep_id, ref -> tuple(cluster_id, 1) }
            .groupTuple()
            .map { cluster_id, ones -> tuple(cluster_id, ones.size()) }

        ch_snippy_gather = SNIPPY_SCATTER.out.sample_dir
            .combine(ch_cluster_sizes, by: 0)
            .map { cluster_id, sid, dir, expected ->
                tuple(groupKey(cluster_id, expected), sid, dir)
            }
            .groupTuple(by: 0)
            .map { key, sample_ids, dirs ->
                def cluster_id = key.toString()
                assert sample_ids.size() == key.size :
                    "cluster ${cluster_id}: expected ${key.size} snippy samples, got ${sample_ids.size()}"
                tuple(cluster_id, sample_ids, dirs)
            }
            .join(
                ch_for_alignment.map { cluster_id, sids, asms, rep_id, ref ->
                    tuple(cluster_id, rep_id, ref)
                }, by: 0
            )
            .map { cluster_id, sample_ids, dirs, rep_id, ref ->
                tuple(cluster_id, sample_ids, dirs, rep_id, ref)
            }

        SNIPPY_CORE_GATHER (
            ch_snippy_gather
        )
        ch_versions = ch_versions.mix(SNIPPY_CORE_GATHER.out.versions.first())
        ch_core_alignments = SNIPPY_CORE_GATHER.out.core_alignment

    } else if (params.alignment_method == 'ska') {
        log.info "Using SKA2 (ska build per sample + ska map) for per-cluster whole genome alignment"

        // Measured, not theoretical -- see modules/local/ska_map_align/main.nf.
        // Split k-mers cannot call a SNP whose flank contains another SNP, so SKA
        // recovers only ~11% of SNPs within 10 bp of a neighbour while matching
        // Snippy beyond ~100 bp. Gubbins detects recombination as elevated SNP
        // density, so this removes the signal: on the 112-genome set SKA produced
        // 54% fewer recombination blocks (2,388 vs 5,148) and disagreed with Snippy
        // on 3 of 4 non-trivial topologies. The result looks plausible, which is
        // what makes it dangerous.
        if (params.run_gubbins) {
            log.warn "ALIGNMENT METHOD 'ska' UNDERCALLS RECOMBINATION. Split k-mers miss " +
                     "clustered SNPs (~11% recovery within 10 bp), which is precisely the " +
                     "signal Gubbins uses. Measured on 112 B. pseudomallei genomes: 54% " +
                     "fewer recombination blocks than snippy and different per-cluster " +
                     "topologies. Use --alignment_method snippy for recombination-aware " +
                     "analysis; see modules/local/ska_map_align/main.nf for the measurements."
        }

        // Same scatter shape as Snippy: per-sample k-mer counting, then a per-cluster
        // reference-anchored map. `ska map` yields a full-length, genome-ORDERED
        // alignment with invariant sites retained (measured: 30 taxa x 376,564
        // columns); `ska align` is unordered and cannot be used as Gubbins input.
        ch_ska_scatter = ch_for_alignment
            .flatMap { cluster_id, sample_ids, assemblies, rep_id, ref ->
                def ids   = sample_ids instanceof List ? sample_ids : [sample_ids]
                def files = assemblies instanceof List ? assemblies : [assemblies]
                assert ids.size() == files.size() :
                    "cluster ${cluster_id}: ${ids.size()} sample ids but ${files.size()} assemblies"
                [ids, files].transpose().sort { it[0] }
                    .collect { sid, asm -> tuple(cluster_id, sid, asm) }
            }

        SKA_BUILD_SAMPLE (
            ch_ska_scatter
        )
        ch_versions = ch_versions.mix(SKA_BUILD_SAMPLE.out.versions.first())

        ch_ska_gather = SKA_BUILD_SAMPLE.out.skf
            .groupTuple(by: 0)
            .join(
                ch_for_alignment.map { cluster_id, sids, asms, rep_id, ref ->
                    tuple(cluster_id, rep_id, ref)
                }, by: 0
            )
            .map { cluster_id, sample_ids, skfs, rep_id, ref ->
                tuple(cluster_id, sample_ids, skfs, rep_id, ref)
            }

        SKA_MAP_ALIGN (
            ch_ska_gather
        )
        ch_versions = ch_versions.mix(SKA_MAP_ALIGN.out.versions.first())
        ch_core_alignments = SKA_MAP_ALIGN.out.core_alignment

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
        error "Unknown alignment method: ${params.alignment_method}. Use 'snippy', 'ska', or 'parsnp'"
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

    // Prepare input for Gubbins: join alignments with starting trees.
    // With gubbins_skip_starting_tree, substitute the 0-byte assets/NO_FILE so
    // GUBBINS_CLUSTER's `[ -s "$starting_tree" ]` guard falls through to the
    // no-starting-tree branch and Gubbins builds its own first tree, which is
    // what the production analysis does (it passes no --starting-tree).
    ch_for_gubbins = params.gubbins_skip_starting_tree
        ? KEEP_INVARIANT_ATCG.out.core_alignment
              .map { cluster_id, alignment ->
                  tuple(cluster_id, alignment, file("${projectDir}/assets/NO_FILE"))
              }
        : KEEP_INVARIANT_ATCG.out.core_alignment
              .join(IQTREE_FAST.out.tree, by: 0)

    GUBBINS_CLUSTER(
        ch_for_gubbins
    )
    ch_versions = ch_versions.mix(GUBBINS_CLUSTER.out.versions)

    /*
    ================================================================================
                    STEP 4: Per-cluster final ML tree
    ================================================================================
    */

    log.info "STEP 4: Building per-cluster final ML trees with ASC correction"

    // Combine Gubbins filtered SNPs with representative info for IQ-TREE
    ch_for_final_tree = GUBBINS_CLUSTER.out.filtered_alignment
        .join(ch_cluster_repid, by: 0)

    // The ASC decision is computed in its own process because the IQ-TREE 2.2.6
    // container ships no python; see modules/local/asc_preflight/main.nf.
    ASC_PREFLIGHT (
        ch_for_final_tree
    )
    ch_versions = ch_versions.mix(ASC_PREFLIGHT.out.versions)

    IQTREE_ASC (
        ch_for_final_tree.join(ASC_PREFLIGHT.out.decision, by: 0)
    )
    ch_versions = ch_versions.mix(IQTREE_ASC.out.versions)

    /*
    ================================================================================
                    STEP 4b: Global ML tree across units, with branch support
    ================================================================================
    */

    // One medoid per unit -> parsnp core alignment -> IQ-TREE with UFBoot and
    // SH-aLRT. This runs in BOTH modes, including curated.
    //
    // The older comment below said curated mode "has no medoids". Only half of
    // that was true: its REFERENCE is external and cannot be a backbone tip, but
    // every unit still has a most-typical member. What curated mode lacks is the
    // Mash matrix the default path used to find one -- and it does not need it,
    // because by this point each unit has a Gubbins filtered polymorphic-sites
    // alignment, and centrality measured on the clonal frame is both cheaper and
    // more apt than whole-genome k-mer distance.
    //
    // Separate from BUILD_BACKBONE_TREE, which runs parsnp with --use-fasttree
    // and therefore yields a backbone with NO support values while every
    // per-cluster tree has them.
    if (params.global_ml_tree == null ? true : params.global_ml_tree) {
        // Reuse the same tuple IQ-TREE consumes: (cluster_id, filtered_aln) joined
        // back to the unit's sample ids and assemblies, so the medoid can be
        // resolved to a genome without re-reading the samplesheet.
        // ch_unit_assemblies was forked off ch_for_alignment before any
        // alignment consumer read it (see STEP 1b). Do NOT read ch_for_alignment
        // here: it is a queue channel and a second reader steals its emissions.
        ch_medoid_in = GUBBINS_CLUSTER.out.filtered_alignment
            .join(ch_unit_assemblies, by: 0)

        SELECT_UNIT_MEDOID ( ch_medoid_in )
        ch_versions = ch_versions.mix(SELECT_UNIT_MEDOID.out.versions.first())

        // Replicon-split runs emit one unit per replicon, and both replicons of a
        // unit would contribute the SAME genome under a different tip name. Take
        // the first replicon's medoid per base cluster so each unit is one tip.
        ch_global_medoids = SELECT_UNIT_MEDOID.out.medoid
            .map { cid, fa -> tuple(cid.replaceFirst(/__.*$/, ''), fa) }
            .groupTuple(by: 0)
            .map { base, fas -> fas.flatten().sort { it.name }.first() }
            .collect()

        GLOBAL_CORE_ALIGNMENT ( ch_global_medoids )
        ch_versions = ch_versions.mix(GLOBAL_CORE_ALIGNMENT.out.versions)

        GLOBAL_ML_TREE ( GLOBAL_CORE_ALIGNMENT.out.alignment )
        ch_versions = ch_versions.mix(GLOBAL_ML_TREE.out.versions)
    }

    /*
    ================================================================================
                    STEP 5: Select representatives per cluster
    ================================================================================
    */

    // STEP 5 + 6 build the parsnp/FastTree backbone used for GRAFTING. It runs
    // only in the default Mash path: grafting needs a backbone whose tips are
    // cluster representatives from the same distance space as the subtrees.
    // The supported global ML tree above is a separate, better-supported product
    // and is not a graft target.
    if (!curated_mode) {

    log.info "STEP 5: Collecting cluster representatives"

    // Collect all representative files. The rep-id text files travel alongside the
    // FASTAs because the module needs the real representative label for each
    // sequence header -- it used to derive that from the filename, which produced
    // one tip named `representative` per cluster.
    COLLECT_REPRESENTATIVES (
        SELECT_CLUSTER_REPRESENTATIVE.out.representative.map { cluster_id, rep_id_file, rep_fasta -> rep_fasta }.collect(),
        SELECT_CLUSTER_REPRESENTATIVE.out.representative.map { cluster_id, rep_id_file, rep_fasta -> rep_id_file }.collect(),
        ch_clusters_file
    )
    ch_versions = ch_versions.mix(COLLECT_REPRESENTATIVES.out.versions)

    /*
    ================================================================================
                    STEP 6: Backbone tree on representatives
    ================================================================================
    */

    log.info "STEP 6: Building backbone tree from representatives"

    BUILD_BACKBONE_TREE (
        COLLECT_REPRESENTATIVES.out.representatives_fasta,
        ch_reference
    )
    ch_versions = ch_versions.mix(BUILD_BACKBONE_TREE.out.versions)

    }  // end STEP 5+6 (default mode only)

    /*
    ================================================================================
                    STEP 6b: Roll up per-cluster phylogeny status
    ================================================================================
    */

    log.info "STEP 6b: Aggregating per-cluster diagnostics into cluster_phylogeny_summary.csv"

    SUMMARIZE_CLUSTER_PHYLOGENY (
        ch_clusters_file,
        GUBBINS_CLUSTER.out.diagnostics.collect(),
        IQTREE_ASC.out.final_tree.map { cluster_id, treefile, rep_id -> treefile }.collect(),
        IQTREE_ASC.out.log.map { cluster_id, iqtree_log -> iqtree_log }.collect(),
        GUBBINS_CLUSTER.out.filtered_alignment.map { cluster_id, aln -> aln }.collect(),
        GUBBINS_CLUSTER.out.recombination_gff.map { cluster_id, gff -> gff }.collect()
    )
    ch_versions = ch_versions.mix(SUMMARIZE_CLUSTER_PHYLOGENY.out.versions)

    /*
    ================================================================================
                    STEP 7: Graft cluster trees onto backbone
    ================================================================================
    */

    if (curated_mode) {
        // No backbone in curated mode, so nothing to graft onto. The per-cluster
        // recombination-corrected ML trees are the final product; expose them as
        // ch_final_tree for the completion log.
        log.info "STEP 7: Grafting SKIPPED in curated mode; per-cluster trees are the output."
        ch_final_tree = IQTREE_ASC.out.final_tree.map { cluster_id, treefile, rep_id -> treefile }
    } else if (params.enable_grafting) {
        log.info "STEP 7: Grafting per-cluster trees onto backbone via graft_trees.py"

        // Build cluster_representatives.tsv from SELECT_CLUSTER_REPRESENTATIVE
        // outputs so graft_trees.py has a deterministic cluster_id -> rep_label
        // mapping (instead of inferring from cluster tip labels).
        ch_reps_tsv = SELECT_CLUSTER_REPRESENTATIVE.out.representative
            .map { cluster_id, rep_id_file, rep_fasta ->
                def rep_label = rep_id_file.text.trim()
                "${cluster_id}\t${rep_label}\n"
            }
            .collectFile(name: 'cluster_representatives.tsv', newLine: false)

        GRAFT_TREES (
            BUILD_BACKBONE_TREE.out.backbone_tree,
            IQTREE_ASC.out.final_tree.map { cluster_id, treefile, rep_id -> treefile }.collect(),
            ch_reps_tsv
        )
        ch_versions   = ch_versions.mix(GRAFT_TREES.out.versions)
        ch_final_tree = GRAFT_TREES.out.grafted_tree
    } else {
        log.info "STEP 7: Tree grafting disabled (enable_grafting=false) — using backbone tree as final tree"
        ch_final_tree = BUILD_BACKBONE_TREE.out.backbone_tree
    }

    /*
    ================================================================================
                        Collect QC information
    ================================================================================
    */

    // Collect QC file check information
    ch_qc_filecheck = ch_qc_filecheck
                        .map{ meta_ignore, file -> file }
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

    // Both strings used to be hardcoded, so a completed run with grafting ENABLED
    // still reported "(grafting disabled)" and called its grafted output a
    // backbone tree -- the log said the opposite of what the run had done.
    ch_final_tree
        .subscribe { tree ->
            log.info "✅ WORKFLOW COMPLETE!"
            if (params.enable_grafting) {
                log.info "Final tree (grafted): ${tree}"
                log.info "Recombination-aware per-cluster analysis finished (per-cluster trees grafted onto backbone)"
            } else {
                log.info "Final tree (backbone): ${tree}"
                log.info "Recombination-aware per-cluster analysis finished (grafting disabled)"
            }
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