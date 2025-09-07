#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GUBBINS_CLUSTER {
  tag "cluster_${cluster_id}"
  label 'process_high'
  container "quay.io/biocontainers/gubbins:3.3.5--py39pl5321he4a0461_0"

  publishDir "${params.outdir}/Clusters/cluster_${cluster_id}/Gubbins",
             mode: params.publish_dir_mode,
             pattern: "*.{fasta,gff,tre,log}"

  input:
    tuple val(cluster_id), path(alignment), path(starting_tree)

  output:
    tuple val(cluster_id), path("${cluster_id}.filtered_polymorphic_sites.fasta"), emit: filtered_alignment
    tuple val(cluster_id), path("${cluster_id}.recombination_predictions.gff"),    emit: recombination_gff
    tuple val(cluster_id), path("${cluster_id}.node_labelled.final_tree.tre"),     emit: final_tree
    path "${cluster_id}.diagnostics.log",                                          emit: diagnostics
    path "versions.yml",                                                           emit: versions

  when:
    task.ext.when == null || task.ext.when

  script:
  // -------- params / flags --------
  def args            = (task.ext.args ?: '').toString().trim()
  def iterations      = (params.gubbins_iterations ?: 3) as int
  def tree_builder    = (params.gubbins_tree_builder ?: 'iqtree').toString()
  def first_builder   = (params.gubbins_first_tree_builder ?: 'rapidnj').toString()
  def min_snps        = (params.gubbins_min_snps ?: 2) as int
  def use_hybrid_flag = ((params.gubbins_use_hybrid == null)
                          ? true
                          : params.gubbins_use_hybrid.toString().toLowerCase() in ['true','1','yes'])

  """
  set -euo pipefail

  DIAG="${cluster_id}.diagnostics.log"
  : > "\${DIAG}"  # start fresh

  # --- Diagnostics header ---
  {
    echo "=== Gubbins Diagnostics for cluster ${cluster_id} ==="
    echo "Container: gubbins 3.3.5 | CPUs: ${task.cpus}"
    echo "Params: iterations=${iterations}, tree_builder=${tree_builder}, first_builder=${first_builder}, min_snps=${min_snps}, use_hybrid=${use_hybrid_flag}"
  } >> "\${DIAG}"

  # --- Input checks ---
  if [ ! -s "${alignment}" ]; then
    echo "ERROR: Alignment file missing or empty: ${alignment}" >> "\${DIAG}"
    : > "${cluster_id}.filtered_polymorphic_sites.fasta"
    : > "${cluster_id}.recombination_predictions.gff"
    : > "${cluster_id}.node_labelled.final_tree.tre"
    cat > versions.yml <<END_VERSIONS
"${task.process}":
    gubbins: \$(run_gubbins.py --version | sed 's/^/    /')
END_VERSIONS
    exit 0
  fi

  seq_count=\$(grep -c "^>" "${alignment}" || echo 0)
  aln_len=\$(awk '/^[^>]/ {gsub(/[ \\t\\r\\n]/,"",\$0); sum += length(\$0)} END {print (sum ? sum : 0)}' "${alignment}")
  {
    echo "Sequence count: \$seq_count"
    echo "Alignment length (non-header chars): \$aln_len"
  } >> "\${DIAG}"

  if [ -e "${starting_tree}" ]; then
    st_bytes=\$(wc -c < "${starting_tree}" || echo 0)
    echo "Starting tree present: ${starting_tree} (\${st_bytes} bytes)" >> "\${DIAG}"
  else
    echo "Starting tree: (none provided)" >> "\${DIAG}"
  fi

  # --- Too-few sequences guard ---
  if [ "\$seq_count" -lt 3 ]; then
    echo "WARNING: Only \$seq_count sequences; skipping Gubbins." >> "\${DIAG}"
    : > "${cluster_id}.filtered_polymorphic_sites.fasta"
    : > "${cluster_id}.recombination_predictions.gff"
    : > "${cluster_id}.node_labelled.final_tree.tre"
    cat > versions.yml <<END_VERSIONS
"${task.process}":
    gubbins: \$(run_gubbins.py --version | sed 's/^/    /')
END_VERSIONS
    exit 0
  fi

  # --- Decide CLI based on starting tree + hybrid flag ---
  has_starting_tree="false"
  if [ -s "${starting_tree}" ]; then
    has_starting_tree="true"
  fi
  echo "has_starting_tree=\$has_starting_tree" >> "\${DIAG}"

  set +e
  if [ "\$has_starting_tree" = "true" ]; then
    echo "RUN: starting-tree + single tree_builder=${tree_builder}" >> "\${DIAG}"
    run_gubbins.py \\
      --starting-tree "${starting_tree}" \\
      --prefix "${cluster_id}" \\
      --tree-builder "${tree_builder}" \\
      --iterations ${iterations} \\
      --min-snps ${min_snps} \\
      --threads ${task.cpus} \\
      ${args} \\
      "${alignment}" >> "\${DIAG}" 2>&1
    gubbins_exit_code=\$?
  else
    if ${use_hybrid_flag}; then
      echo "RUN: no starting-tree + HYBRID first_builder=${first_builder}, tree_builder=${tree_builder}" >> "\${DIAG}"
      run_gubbins.py \\
        --prefix "${cluster_id}" \\
        --first-tree-builder "${first_builder}" \\
        --tree-builder "${tree_builder}" \\
        --iterations ${iterations} \\
        --min-snps ${min_snps} \\
        --threads ${task.cpus} \\
        ${args} \\
        "${alignment}" >> "\${DIAG}" 2>&1
      gubbins_exit_code=\$?
    else
      echo "RUN: no starting-tree + NON-HYBRID tree_builder=${first_builder}" >> "\${DIAG}"
      run_gubbins.py \\
        --prefix "${cluster_id}" \\
        --tree-builder "${first_builder}" \\
        --iterations ${iterations} \\
        --min-snps ${min_snps} \\
        --threads ${task.cpus} \\
        ${args} \\
        "${alignment}" >> "\${DIAG}" 2>&1
      gubbins_exit_code=\$?
    fi
  fi
  set -e

  echo "Gubbins exit code: \${gubbins_exit_code:-NA}" >> "\${DIAG}"
  if [ "\${gubbins_exit_code:-1}" -ne 0 ]; then
    echo "ERROR: Gubbins failed; creating placeholder outputs." >> "\${DIAG}"
    : > "${cluster_id}.filtered_polymorphic_sites.fasta"
    : > "${cluster_id}.recombination_predictions.gff"
    : > "${cluster_id}.node_labelled.final_tree.tre"
  fi

  # --- Output integrity checks ---
  for f in "${cluster_id}.filtered_polymorphic_sites.fasta" \\
           "${cluster_id}.recombination_predictions.gff" \\
           "${cluster_id}.node_labelled.final_tree.tre"; do
    if [ ! -s "\$f" ]; then
      echo "WARNING: Output file \$f is empty; leaving placeholder." >> "\${DIAG}"
      : > "\$f"
    else
      sz=\$(wc -c < "\$f" || echo 0)
      echo "Output file \$f size: \$sz bytes" >> "\${DIAG}"
    fi
  done

  # --- Versions file (Nextflow expands \${task.process}, Bash runs \$(...)) ---
  cat > versions.yml <<END_VERSIONS
"${task.process}":
    gubbins: \$(run_gubbins.py --version | sed 's/^/    /')
END_VERSIONS
  """
}
