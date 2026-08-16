#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * GLOBAL_CORE_ALIGNMENT
 *
 * Core-genome SNP alignment across the analysis units: parsnp over one medoid
 * per unit. Feeds GLOBAL_ML_TREE.
 *
 * Split from the tree step because the containers are disjoint -- the parsnp
 * biocontainer carries parsnp, harvesttools and fasttree but no IQ-TREE -- and
 * because an alignment is worth publishing on its own.
 *
 * `--skip-phylogeny` is deliberate. Only the alignment is wanted, and parsnp's
 * own tree step was measured hanging for over twenty minutes after the alignment
 * was already complete on 82 medoids. The tree comes from IQ-TREE downstream,
 * with support values parsnp's FastTree does not provide.
 */

process GLOBAL_CORE_ALIGNMENT {

    tag "global_core_alignment"
    label 'process_high'
    container "quay.io/biocontainers/parsnp:1.7.4--hdcf5f25_2"

    publishDir "${params.outdir}", mode: params.publish_dir_mode, pattern: "global_ml_alignment.fasta"
    publishDir "${params.outdir}", mode: params.publish_dir_mode, pattern: "global_core_alignment.log"

    input:
    path medoids, stageAs: 'medoids/*'

    output:
    path "global_ml_alignment.fasta"  , emit: alignment
    path "global_core_alignment.log"  , emit: log
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail
    LOG=global_core_alignment.log
    : > "\$LOG"

    n=\$(ls medoids/ | wc -l)
    echo "unit medoids supplied: \$n" | tee -a "\$LOG"
    if [ "\$n" -lt 3 ]; then
        echo "ERROR: only \$n unit medoids; a global tree needs at least 3." | tee -a "\$LOG" >&2
        exit 1
    fi

    # Anchor on the medoid with the fewest contigs so the core is as large as the
    # data allow: parsnp's core shrinks toward the most fragmented member.
    REF=\$(for f in medoids/*; do printf '%s\\t%s\\n' "\$(grep -c '^>' "\$f")" "\$f"; done | sort -n | head -1 | cut -f2)
    echo "parsnp reference (fewest contigs): \$REF" | tee -a "\$LOG"

    set +e
    parsnp -r "\$REF" -d medoids -o parsnp_out -p ${task.cpus} -c --skip-phylogeny >> "\$LOG" 2>&1
    rc=\$?
    set -e

    if [ ! -s parsnp_out/parsnp.snps.mblocks ]; then
        echo "ERROR: parsnp exited \$rc and produced no SNP alignment." | tee -a "\$LOG" >&2
        tail -30 parsnp_out/parsnpAligner.log >&2 2>/dev/null || true
        exit 1
    fi

    # parsnp labels each tip with the staged FILENAME, and the anchor with a
    # trailing '.ref'. Strip the extension (.fa or .fasta) and any .ref suffix so
    # the alignment -- and every tree built from it -- is labelled by unit alone.
    # Getting this wrong is not cosmetic: the label lands in the published tree.
    sed -e 's/\\.ref\$//' -e 's/^>\\(.*\\)\\.\\(fa\\|fasta\\|fna\\)\$/>\\1/' \\
        parsnp_out/parsnp.snps.mblocks > global_ml_alignment.fasta

    ntax=\$(grep -c '^>' global_ml_alignment.fasta)
    nsit=\$(awk '/^>/{if(n)exit;next}{n+=length(\$0)}END{print n+0}' global_ml_alignment.fasta)
    echo "core SNP alignment: \$ntax taxa, \$nsit sites" | tee -a "\$LOG"
    if [ "\$ntax" -lt 3 ]; then
        echo "ERROR: alignment has \$ntax taxa; parsnp dropped too many medoids." | tee -a "\$LOG" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parsnp: \$(parsnp --version 2>&1 | head -1)
    END_VERSIONS
    """
}
