#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * GLOBAL_ML_TREE
 *
 * Maximum-likelihood tree over the analysis units, WITH branch support:
 * IQ-TREE on GLOBAL_CORE_ALIGNMENT's core SNP alignment, UFBoot + SH-aLRT.
 *
 * WHY NOT BUILD_BACKBONE_TREE. That process runs parsnp with --use-fasttree, so
 * the backbone was the only tree in the output carrying no support values while
 * every per-cluster tree had them. This is the ML equivalent.
 *
 * **THE RESULT IS NOT RECOMBINATION-CORRECTED, AND MUST NOT BE.** Gubbins finds
 * recombination as regions of unusually dense SNPs against a clonal background.
 * Across dozens of divergent lineages there is no shared clonal background: the
 * between-lineage differences ARE the dense regions, and Gubbins would call most
 * of the alignment recombinant. That is precisely the failure the per-unit
 * partition exists to prevent. This tree shows how units relate, its branch
 * lengths include recombination, and NO r/m may be derived from it -- r/m is a
 * within-unit measurement.
 *
 * THE +ASC RETRY is not defensive padding. IQ-TREE refuses an ascertainment-bias
 * model when the alignment contains invariant columns, and parsnp's SNP
 * alignment does contain a few once the anchor's own label is dropped: measured
 * 87 invariant of 82,601 columns across 82 medoids. IQ-TREE writes a
 * variable-sites-only file when it refuses, so the retry uses that and records
 * that it happened rather than silently switching model.
 */

process GLOBAL_ML_TREE {

    tag "global_ml_tree"
    label 'process_high'
    container "quay.io/biocontainers/iqtree:2.2.6--h21ec9f0_0"

    publishDir "${params.outdir}", mode: params.publish_dir_mode, pattern: "global_ml_tree.*"

    input:
    path alignment

    output:
    path "global_ml_tree.treefile", emit: tree
    path "global_ml_tree.log"     , emit: log
    path "global_ml_tree.iqtree"  , emit: report, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ufboot = (params.iqtree_ufboot ?: 1000) as int
    // Determinism is governed by either parameter. `deterministic` is the current
    // name; `gubbins_deterministic` predates it and is honoured so the recipe in
    // PROVENANCE.md keeps working. Measured 2026-09-04 on three real unit
    // alignments: -seed alone still gives a different tree, -T 1 alone still
    // gives a different tree, and both together give an identical one.
    def deterministic = [params.deterministic, params.gubbins_deterministic].any {
        it != null && it.toString().toLowerCase() in ['true', '1', 'yes'] }
    def iq_threads  = deterministic ? 1 : task.cpus
    def iq_seed_arg = (params.iqtree_seed == null) ? '' : "-seed ${params.iqtree_seed}"
    def alrt   = (params.iqtree_alrt   ?: 1000) as int
    def model  = (params.global_ml_model ?: 'GTR+ASC').toString()
    """
    set -euo pipefail
    LOG=global_ml_tree.log
    : > "\$LOG"

    ntax=\$(grep -c '^>' ${alignment})
    echo "taxa: \$ntax   model: ${model}   support: SH-aLRT ${alrt} / UFBoot ${ufboot}" | tee -a "\$LOG"
    if [ "\$ntax" -lt 4 ]; then
        # UFBoot needs at least 4 taxa for a resolvable internal branch.
        echo "ERROR: \$ntax taxa is too few for a supported ML tree." | tee -a "\$LOG" >&2
        exit 1
    fi

    set +e
    iqtree2 -s ${alignment} -st DNA -m ${model} -T ${iq_threads} ${iq_seed_arg} \\
            --prefix gml -bb ${ufboot} -alrt ${alrt} >> "\$LOG" 2>&1
    rc=\$?
    set -e

    if [ "\$rc" -ne 0 ] && [ -s gml.varsites.phy ]; then
        echo "NOTE: ${model} rejected because the alignment holds invariant columns;" | tee -a "\$LOG"
        echo "      retrying on gml.varsites.phy (variable sites only)." | tee -a "\$LOG"
        set +e
        iqtree2 -s gml.varsites.phy -st DNA -m ${model} -T ${iq_threads} ${iq_seed_arg} \\
                --prefix gml -bb ${ufboot} -alrt ${alrt} >> "\$LOG" 2>&1
        rc=\$?
        set -e
    fi

    if [ "\$rc" -ne 0 ] || [ ! -s gml.treefile ]; then
        echo "ERROR: IQ-TREE exited \$rc and produced no tree." | tee -a "\$LOG" >&2
        tail -30 "\$LOG" >&2
        exit 1
    fi

    cp gml.treefile global_ml_tree.treefile
    [ -s gml.iqtree ] && cp gml.iqtree global_ml_tree.iqtree || true
    echo "NOT recombination-corrected -- do not derive r/m from this tree" | tee -a "\$LOG"

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    iqtree: \$(iqtree2 --version 2>&1 | head -1)
END_VERSIONS
    """
}
