process IQTREE_ASC {
  tag "cluster_${cluster_id}"
  label 'process_high'
  container "quay.io/biocontainers/iqtree:2.2.6--h21ec9f0_0"

  input:
    tuple val(cluster_id), path(filtered_snps), val(representative_id)

  output:
    tuple val(cluster_id), path("${cluster_id}.final.treefile"), val(representative_id), emit: final_tree
    tuple val(cluster_id), path("${cluster_id}.final.iqtree"), emit: log
    path "versions.yml", emit: versions

  when:
    task.ext.when == null || task.ext.when

  script:
    // Optional extra CLI flags (string)
    def args  = (task.ext.args ?: '').toString().trim()
    // Keep ASC via model string for SNP-only alignments
    def model = (params.iqtree_asc_model ?: 'GTR+ASC').toString().trim()
    // If you want to supply constant-site counts for full alignments, set params.iqtree_fconst="a,c,g,t" or counts "100,120,98,115"
    def fconst = (params.iqtree_fconst ?: '').toString().trim()

    """
    set -euo pipefail

    echo "Building final ML tree for cluster ${cluster_id}"
    echo "Representative: ${representative_id}"

    # Choose IQ-TREE binary (handles images exposing 'iqtree' instead of 'iqtree2')
    IQTREE=\$(command -v iqtree2 || command -v iqtree || true)
    if [ -z "\$IQTREE" ]; then
      echo "ERROR: iqtree/iqtree2 not found in PATH"
      : > ${cluster_id}.final.treefile
      : > ${cluster_id}.final.iqtree
      printf "\"%s\":\\n  iqtree:     N/A\\n" "${task.process}" > versions.yml
      exit 0
    fi

    # Basic sanity checks
    if [ ! -s "${filtered_snps}" ]; then
      echo "WARNING: Empty alignment for ${cluster_id} — emitting minimal outputs"
      : > ${cluster_id}.final.treefile
      : > ${cluster_id}.final.iqtree
      IQVER=\$( "\$IQTREE" -version 2>&1 | head -1 )
      printf "\"%s\":\\n  iqtree:     %s\\n" "${task.process}" "    \$IQVER" > versions.yml
      exit 0
    fi

    seq_count=\$(grep -c "^>" ${filtered_snps} || echo 0)
    if [ "\$seq_count" -lt 3 ]; then
      echo "WARNING: only \$seq_count sequences — IQ-TREE requires >=3"
      : > ${cluster_id}.final.iqtree
      names=\$(grep "^>" ${filtered_snps} | sed 's/^>//' | tr '\\n' ',' | sed 's/,\$//')
      echo "(\$names);" > ${cluster_id}.final.treefile
      IQVER=\$( "\$IQTREE" -version 2>&1 | head -1 )
      printf "\"%s\":\\n  iqtree:     %s\\n" "${task.process}" "    \$IQVER" > versions.yml
      exit 0
    fi

    # Build extra flags as a plain string (POSIX sh-friendly)
    EXTRA="${args}"
    if [ -n "${fconst}" ]; then
      EXTRA="\${EXTRA} -fconst ${fconst}"
    fi

    echo "Running IQ-TREE with model ${model}"
    if ! "\$IQTREE" \
        -s ${filtered_snps} \
        -st DNA \
        -m ${model} \
        -bb 1000 \
        -alrt 1000 \
        -nt AUTO \
        --prefix ${cluster_id}.final \
        \${EXTRA}; then
      echo "WARNING: IQ-TREE failed for ${cluster_id}; emitting minimal tree"
      : > ${cluster_id}.final.iqtree
      names=\$(grep "^>" ${filtered_snps} | sed 's/^>//' | tr '\\n' ',' | sed 's/,\$//')
      echo "(\$names);" > ${cluster_id}.final.treefile
    fi

    # Ensure expected files exist
    [ -f "${cluster_id}.final.treefile" ] || { names=\$(grep "^>" ${filtered_snps} | sed 's/^>//' | tr '\\n' ',' | sed 's/,\$//'); echo "(\$names);" > ${cluster_id}.final.treefile; }
    [ -f "${cluster_id}.final.iqtree" ]   || : > ${cluster_id}.final.iqtree

    # Versions (avoid command substitution inside heredoc)
    IQVER=\$( "\$IQTREE" -version 2>&1 | head -1 )
    printf "\"%s\":\\n  iqtree:     %s\\n" "${task.process}" "    \$IQVER" > versions.yml

    echo "Final ML tree construction completed for cluster ${cluster_id}"
    """
}
