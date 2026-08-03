#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * SELECT_CLUSTER_REPRESENTATIVE
 *
 * Changes relative to the original module:
 *
 * 1. No more runtime `pip install pandas numpy`.  Container pinned to
 *    quay.io/biocontainers/pyseer:1.4.2--pyhdfd78af_0, whose contents were
 *    verified by pulling the image manifest from the quay.io registry API and
 *    listing the conda-meta/ records in its layers:
 *        python 3.8.20, pandas 2.0.3, numpy 1.22.4, scipy 1.7.0,
 *        scikit-learn 1.3.2, mash 2.3
 *    (biocontainers/python:3.9--1, used before, ships none of them.)
 *
 * 2. The 100-line inline Python heredoc moved to bin/select_representative.py.
 *    Heredocs interpolate Nextflow's ${...} into Python source, so every `$` in
 *    the script had to be escaped and the script could not be run or tested
 *    outside a pipeline execution.
 *
 * 3. The matrix input is now the CLUSTER'S OWN submatrix, written once by
 *    CLUSTER_GENOMES (--emit-submatrices).  The original staged and re-parsed
 *    the full n x n mash_matrix.tsv in every per-cluster task: with n = 2000 and
 *    ~40-60 clusters that is 40-60 full parses of a 2000-column TSV to answer a
 *    question about <= max_cluster_size rows.  If a full matrix is passed
 *    instead (e.g. the merged_small_clusters pseudo-cluster, which has no
 *    submatrix), the script falls back to reading it with pandas usecols so only
 *    the cluster's columns are materialized.
 *
 * The medoid computation itself is unchanged: argmin of the NaN-skipping row
 * sums, O(k^2) per cluster, which is negligible at k <= max_cluster_size.
 */

process SELECT_CLUSTER_REPRESENTATIVE {
    tag "cluster_${cluster_id}"
    label 'process_single'
    conda "conda-forge::python=3.10 conda-forge::numpy conda-forge::pandas conda-forge::scipy"
    container "quay.io/biocontainers/pyseer:1.4.2--pyhdfd78af_0"

    input:
    tuple val(cluster_id), val(sample_ids), path(assemblies), path(cluster_matrix)

    output:
    // EXACT filenames, not path("*.fa"): input assemblies are staged into this same
    // work dir, and any input already named *.fa would also match the glob, making
    // `reference` a LIST of two files. Downstream snippy then interpolates both
    // into --ref and fails with "reference is missing or empty: a.fa representative.fa".
    //
    // The names are also CLUSTER-SCOPED. They used to be the bare `representative.fa`
    // / `representative_id.txt`, which is unique within one task dir but collides the
    // moment COLLECT_REPRESENTATIVES stages every cluster's copy into a single dir
    // ("input file name collision"). It also silently corrupted the backbone tree:
    // that module derived each sequence header from `basename <file> .fa`, so every
    // representative was relabelled `representative` and no tip could ever match the
    // cluster_id -> rep_label map GRAFT_TREES joins on.
    tuple val(cluster_id), path("${cluster_id}.rep_id.txt"), path("${cluster_id}.rep.fa"), emit: representative
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${projectDir}/bin/select_representative.py \\
        --cluster-id '${cluster_id}' \\
        --sample-ids '${sample_ids}' \\
        --matrix '${cluster_matrix}' \\
        --out-id ${cluster_id}.rep_id.txt

    # Emit the representative's ASSEMBLY under the cluster-scoped name declared in
    # the output block. Downstream SNIPPY_SCATTER / SKA_MAP_ALIGN consume it as the
    # reference FASTA -- if it does not exist Nextflow binds the rep-id text file in
    # its place and snippy is handed a one-line text file as its reference.
    REP_ID=\$(head -n1 ${cluster_id}.rep_id.txt | tr -d '[:space:]')
    REP_SRC=""
    for f in *.fasta *.fa *.fna *.fas; do
        [ -e "\$f" ] || continue
        case "\$f" in ${cluster_id}.rep.fa) continue ;; esac
        # strip one or two extensions to match the normalized sample id
        b=\$(basename "\$f"); b1="\${b%.*}"; b2="\${b1%.*}"
        if [ "\$b" = "\$REP_ID" ] || [ "\$b1" = "\$REP_ID" ] || [ "\$b2" = "\$REP_ID" ]; then
            REP_SRC="\$f"; break
        fi
    done
    if [ -z "\$REP_SRC" ]; then
        echo "ERROR: no staged assembly matches representative id '\$REP_ID' for ${cluster_id}" >&2
        echo "staged files: \$(ls *.fasta *.fa *.fna *.fas 2>/dev/null | tr '\\n' ' ')" >&2
        exit 1
    fi
    cp "\$REP_SRC" ${cluster_id}.rep.fa
    echo "Representative for ${cluster_id}: \$REP_ID (\$REP_SRC)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
