#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process COLLECT_REPRESENTATIVES {
    tag "collect_reps"
    label 'process_low'
    container "ubuntu:jammy"
    // Publish the representative FASTA and mapping at the run outdir root so
    // downstream scripts can find `cluster_representatives.tsv` and
    // `representatives.fa` directly under the pipeline output directory.
    publishDir "${params.outdir}", mode: params.publish_dir_mode, pattern: "representatives.fa"
    publishDir "${params.outdir}", mode: params.publish_dir_mode, pattern: "cluster_representatives.tsv"

    input:
    // Both are collected across every cluster, so they MUST arrive under
    // cluster-scoped names (`<cluster_id>.rep.fa` / `<cluster_id>.rep_id.txt`).
    // SELECT_CLUSTER_REPRESENTATIVE emitted bare `representative.fa` /
    // `representative_id.txt` until this was fixed, which made Nextflow abort with
    // "input file name collision" as soon as there were two clusters.
    path rep_fastas
    path rep_id_files
    path cluster_info

    output:
    path "representatives.fa", emit: representatives_fasta
    path "cluster_representatives.tsv", emit: representatives_mapping
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    echo "Collecting cluster representatives into single FASTA file"

    > representatives.fa
    printf 'cluster_id\\trepresentative_id\\tsource_file\\n' > cluster_representatives.tsv

    rep_count=\$(ls *.rep.fa 2>/dev/null | wc -l)
    echo "Found \$rep_count representative FASTA file(s)"

    if [ "\$rep_count" -eq 0 ]; then
        echo "ERROR: no *.rep.fa staged -- SELECT_CLUSTER_REPRESENTATIVE produced nothing." >&2
        echo "staged files: \$(ls -A | tr '\\n' ' ')" >&2
        exit 1
    fi

    for rep_file in *.rep.fa; do
        # cluster_id is carried by the filename, which SELECT_CLUSTER_REPRESENTATIVE
        # controls; it is not guessed from the sample name the way it used to be.
        cluster_id="\${rep_file%.rep.fa}"
        id_file="\${cluster_id}.rep_id.txt"

        if [ ! -s "\$id_file" ]; then
            echo "ERROR: missing or empty \$id_file for \$cluster_id" >&2
            exit 1
        fi

        # The REAL representative label, as written by select_representative.py.
        # This must match the rep_label that the workflow feeds GRAFT_TREES in
        # cluster_representatives.tsv, otherwise no backbone tip can be matched to
        # its cluster subtree. Deriving it from the filename (the old behaviour)
        # relabelled every representative `representative`.
        rep_id=\$(head -n1 "\$id_file" | tr -d '[:space:]')
        if [ -z "\$rep_id" ]; then
            echo "ERROR: \$id_file contains no representative id for \$cluster_id" >&2
            exit 1
        fi

        if [ ! -s "\$rep_file" ]; then
            echo "ERROR: empty representative FASTA for \$cluster_id (\$rep_file)" >&2
            exit 1
        fi

        echo "Cluster \$cluster_id -> representative \$rep_id"

        # Every contig header becomes the representative label. Draft assemblies
        # here carry a median of ~91 contigs; BUILD_BACKBONE_TREE splits this file
        # back apart on the header, so all contigs of one representative must share
        # that one label for them to regroup into a single per-rep FASTA.
        awk -v rep_id="\$rep_id" '
            /^>/ { print ">" rep_id; next }
            { print }
        ' "\$rep_file" >> representatives.fa

        printf '%s\\t%s\\t%s\\n' "\$cluster_id" "\$rep_id" "\$rep_file" >> cluster_representatives.tsv
    done

    # Distinct labels, not sequence count: a multi-contig representative writes one
    # '>' per contig, so grep -c '^>' would far exceed the number of clusters.
    final_count=\$(grep '^>' representatives.fa | sort -u | wc -l)
    echo "Collected \$final_count distinct representative(s) from \$rep_count cluster(s)"

    if [ "\$final_count" -ne "\$rep_count" ]; then
        # Previously this fabricated a '>dummy_representative / ATCG' entry and
        # exited 0, which poisoned the backbone tree the same way the star-tree
        # fallback poisoned cluster trees. Failures fail now.
        echo "ERROR: \$rep_count clusters produced \$final_count distinct labels -- duplicate or missing representative ids." >&2
        grep '^>' representatives.fa | sort | uniq -d >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ubuntu: \$(awk -F ' ' '{print \$2,\$3}' /etc/issue | tr -d '\n')
END_VERSIONS
    """
}
