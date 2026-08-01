#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IQTREE_FAST -- THROWAWAY starting tree for Gubbins
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    This tree is consumed by `run_gubbins.py --starting-tree` and then discarded:
    Gubbins re-estimates the tree at every iteration. It does not need to be a good
    tree, only a reasonable topology.

    WHAT CHANGED
    1. The old script computed `def model = params.iqtree_model ?: 'GTR+ASC'` and then
       never used it -- it hardcoded `-m MFP`. On a 40 taxa x 3000 site alignment,
       IQ-TREE reported "ModelFinder will test up to 968 DNA models" and then selected
       plain JC: 968 model fits to seed a tree that gets thrown away.
       Now the model comes from params.iqtree_starting_model (default 'GTR+G'),
       or the whole tree can be delegated to a distance/parsimony builder via
       params.starting_tree_builder ('iqtree' | 'veryfasttree' | 'rapidnj' | 'fasttree').
    2. `-nt AUTO` -> `-T ${task.cpus}`. AUTO spends ~57 s per task probing thread
       counts, and it logged "auto-detect threads (22 CPU cores detected)" -- it sizes
       to the HOST, not to task.cpus, so under maxForks>1 every concurrent task tries
       to use the whole machine. Measured on 40 taxa x 3000 sites (22-core host):
           -m MFP   -nt AUTO  (old)   60.60 s
           -m MFP   -T 4               3.96 s
           -m GTR+G -nt AUTO          60.73 s
           -m GTR+G -T 4  (new)        0.68 s
    3. FIVE identical publishDir directives collapsed to one.
    4. The `|| { touch treefile; }` fallback is gone for the failure case: an empty
       starting tree makes GUBBINS_CLUSTER silently take its no-starting-tree branch,
       so a failure here changes what Gubbins does without any error surfacing. The
       >=3-sequence guard is KEPT and still emits an empty treefile, because that is
       the intended, documented signal for a cluster too small to seed a tree.
*/

process IQTREE_FAST {
    tag "cluster_${cluster_id}"
    label 'process_medium'
    container "quay.io/biocontainers/iqtree:2.2.6--h21ec9f0_0"

    publishDir "${params.outdir}/Clusters/cluster_${cluster_id}", mode: params.publish_dir_mode, pattern: "*.{treefile,iqtree}"

    input:
    tuple val(cluster_id), path(alignment)

    output:
    tuple val(cluster_id), path("${cluster_id}.treefile"), emit: tree
    tuple val(cluster_id), path("${cluster_id}.iqtree"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    // Cheap FIXED model. This tree is a throwaway Gubbins seed; ModelFinder (-m MFP)
    // tested up to 968 models here and picked JC anyway.
    def model   = (params.iqtree_starting_model ?: 'GTR+G').toString().trim()
    def builder = (params.starting_tree_builder ?: 'iqtree').toString().trim()
    """
    set -euo pipefail

    # Check if alignment has at least 3 sequences (minimum for tree building).
    # An empty treefile is the documented signal for "cluster too small to seed";
    # GUBBINS_CLUSTER branches on `[ -s starting_tree ]`.
    seq_count=\$(grep -c "^>" ${alignment} || echo 0)
    if [ "\$seq_count" -lt 3 ]; then
        echo "WARNING: alignment has only \$seq_count sequences; IQ-TREE needs >=3."
        echo "Emitting an empty starting tree for cluster ${cluster_id}; Gubbins will build its own."
        : > ${cluster_id}.treefile
        : > ${cluster_id}.iqtree
        cat > versions.yml <<END_VERSIONS
"${task.process}":
    starting_tree_builder: ${builder}
END_VERSIONS
        exit 0
    fi

    BUILDER="${builder}"
    echo "Starting-tree builder: \$BUILDER (model ${model} when builder=iqtree)"

    case "\$BUILDER" in
      iqtree)
        # -T \${task.cpus}: NOT -nt AUTO. AUTO costs ~57 s/task probing thread counts
        # and sizes itself to the HOST core count, ignoring task.cpus, which
        # oversubscribes the machine under maxForks > 1.
        iqtree2 \\
            -s ${alignment} \\
            -st DNA \\
            -m ${model} \\
            --fast \\
            -T ${task.cpus} \\
            --prefix ${cluster_id} \\
            ${args}
        ;;
      veryfasttree|fasttree)
        # Distance/ML-lite alternative. Measured inside Gubbins as --first-tree-builder,
        # veryfasttree was ~14% faster than rapidnj on n=50 (371.1 s vs 431.0 s).
        BIN=\$(command -v VeryFastTree || command -v veryfasttree || command -v FastTree || command -v fasttree)
        "\$BIN" -nt -gtr ${args} ${alignment} > ${cluster_id}.treefile
        echo "starting tree built by \$BIN" > ${cluster_id}.iqtree
        ;;
      rapidnj)
        rapidnj ${alignment} -i fa -t d -x ${cluster_id}.treefile ${args}
        echo "starting tree built by rapidnj" > ${cluster_id}.iqtree
        ;;
      *)
        echo "ERROR: unknown params.starting_tree_builder='\$BUILDER' (expected iqtree|veryfasttree|fasttree|rapidnj)" >&2
        exit 1
        ;;
    esac

    # A non-empty starting tree is required, because an EMPTY one silently flips
    # GUBBINS_CLUSTER onto its no-starting-tree branch instead of reporting a failure.
    if [ ! -s "${cluster_id}.treefile" ]; then
        echo "ERROR: starting tree for cluster ${cluster_id} is empty after a successful build" >&2
        exit 1
    fi
    [ -f "${cluster_id}.iqtree" ] || : > ${cluster_id}.iqtree

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    starting_tree_builder: ${builder}
    iqtree: \$(iqtree2 --version 2>&1 | head -n1 | sed 's/^IQ-TREE multicore version //;s/ for .*//')
END_VERSIONS
    """
}
