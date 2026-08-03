#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * ASC_PREFLIGHT
 *
 * Decides how IQ-TREE should handle ascertainment-bias correction for one
 * cluster, and does it in a container that HAS PYTHON.
 *
 * Why this is a separate process at all: the work used to run inline inside
 * IQTREE_ASC, which declares
 *
 *     conda     "bioconda::iqtree=2.2.6 conda-forge::python=3.10 ..."
 *     container "quay.io/biocontainers/iqtree:2.2.6--h21ec9f0_0"
 *
 * The conda spec pulls python in; the container has NO python and no python3 at
 * all. So `asc_preflight.py` ran fine under -profile conda and died with
 * `env: can't execute 'python3': No such file or directory` under
 * -profile docker -- the module could never have worked in a container. Rather
 * than move IQ-TREE onto some image that happens to bundle python (the obvious
 * candidate, the Gubbins container, ships IQ-TREE 2.3.0 and would silently
 * change the version that builds the final publishable trees), the python step
 * gets its own container and IQ-TREE stays pinned at 2.2.6.
 *
 * The container here is the same pyseer image SELECT_CLUSTER_REPRESENTATIVE
 * already uses, whose contents were verified against the quay.io manifest.
 */

process ASC_PREFLIGHT {
  tag "cluster_${cluster_id}"
  label 'process_single'
  conda "conda-forge::python=3.10 conda-forge::numpy"
  container "quay.io/biocontainers/pyseer:1.4.2--pyhdfd78af_0"

  input:
    tuple val(cluster_id), path(filtered_snps), val(representative_id)

  output:
    // Deliberately NOT named `<cluster_id>.asc_preflight.txt`: IQTREE_ASC both
    // stages this file as an input and publishes a file of that name, and
    // Nextflow refuses to have an input and an output share a name.
    tuple val(cluster_id), path("${cluster_id}.asc_decision.txt"),
                           path("${cluster_id}.asc_stripped.fasta"), emit: decision
    path "versions.yml", emit: versions

  when:
    task.ext.when == null || task.ext.when

  script:
    def model  = (params.iqtree_asc_model ?: 'GTR+ASC').toString().trim()
    def asc_fb = (params.iqtree_asc_fallback ?: 'varsites').toString().trim()
    """
    set -euo pipefail

    if [ ! -s "${filtered_snps}" ]; then
      echo "ERROR: Gubbins filtered alignment for cluster ${cluster_id} is missing or empty: ${filtered_snps}" >&2
      exit 1
    fi

    seq_count=\$(grep -c "^>" "${filtered_snps}" || echo 0)

    if [ "\$seq_count" -lt 3 ]; then
      # <3 tips means no topology exists. Record the decision and let IQTREE_ASC
      # emit the degenerate Newick -- this is the one remaining place a degenerate
      # tree is legitimate, and DEGENERATE_TREE=1 is what lets the summary tell it
      # apart from a real inference.
      echo "WARNING: cluster ${cluster_id} has \$seq_count sequences; <3 tips means no topology exists."
      printf 'N_TAXA=%s\\nASC_ACTION=skipped_too_few_taxa\\nDEGENERATE_TREE=1\\nIQ_MODEL=\\nIQ_ALIGNMENT=%s\\nIQ_FCONST=\\n' \\
        "\$seq_count" "${filtered_snps}" > ${cluster_id}.asc_decision.txt
      : > ${cluster_id}.asc_stripped.fasta
    else
      asc_preflight.py \\
        --input "${filtered_snps}" \\
        --out "${cluster_id}.asc_decision.txt" \\
        --model "${model}" \\
        --strategy "${asc_fb}" \\
        --stripped-output "${cluster_id}.asc_stripped.fasta"

      # The script only writes the stripped alignment when it actually strips
      # columns, but the output block declares it unconditionally.
      [ -e "${cluster_id}.asc_stripped.fasta" ] || : > ${cluster_id}.asc_stripped.fasta
    fi

    cat > versions.yml <<END_VERSIONS
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
END_VERSIONS
    """
}
