/* -*- coding: utf-8 -*- */
nextflow.enable.dsl=2

process BUILD_BACKBONE_TREE {
  tag "backbone_tree"
  label 'process_high'

  // Override via params.backbone_container in your config if you’d like.
  // This biocontainer typically includes parsnp, harvesttools, fasttree.
  container "${params.backbone_container ?: 'quay.io/biocontainers/parsnp:1.7.4--pyhdfd78af_0'}"

  input:
  path representatives_fasta
  path mash_distances, optional: true

  output:
  path "backbone.treefile",       emit: backbone_tree
  path "backbone_alignment.fa",   emit: backbone_alignment
  path "backbone_report.txt",     emit: report
  path "versions.yml",            emit: versions

  script:
  def method   = params.backbone_method ?: 'parsnp'
  def ft_opts  = params.backbone_fasttree_opts ?: '-nt -gtr'
  def threads  = task.cpus
  """
  set -euo pipefail

  echo "== BUILD_BACKBONE_TREE =="
  echo "Method: ${method}"
  echo "CPUs: ${threads}"

  rep_count=\$(grep -c '^>' "${representatives_fasta}" || echo 0)
  echo "Representatives: \$rep_count"

  status="SUCCESS"

  # Handle tiny sets gracefully
  if [ "\$rep_count" -lt 3 ]; then
    echo "WARNING: <3 representatives found; emitting trivial tree."
    awk '/^>/{gsub(/^>/,""); n=n (n? ",":"") \$0} END{print "(" n ");"}' "${representatives_fasta}" > backbone.treefile || echo "(A:0.1,B:0.1);" > backbone.treefile
    cp "${representatives_fasta}" backbone_alignment.fa
    status="TRIVIAL"
  else
    if [ "${method}" = "parsnp" ]; then
      echo "Running Parsnp backbone..."
      mkdir -p reps
      awk '/^>/{fn="reps/" substr(\$0,2) ".fa"}{print > fn}' "${representatives_fasta}"
      ref=\$(ls reps/*.fa | head -n1 || true)
      if [ -z "\$ref" ]; then
        echo "No reference resolved; falling back to FastTree."
        cp "${representatives_fasta}" backbone_alignment.fa
        if command -v fasttree >/dev/null 2>&1; then
          fasttree ${ft_opts} backbone_alignment.fa > backbone.treefile || status="FALLBACK"
        else
          status="FALLBACK"
        fi
      else
        set +e
        parsnp --sequences reps --reference "\$ref" --output-dir parsnp_backbone --use-fasttree --threads ${threads} --verbose
        rc=\$?
        set -e
        if [ "\$rc" -eq 0 ] && [ -s parsnp_backbone/parsnp.tree ]; then
          cp parsnp_backbone/parsnp.tree backbone.treefile
          if [ -s parsnp_backbone/parsnp.xmfa ]; then
            cp parsnp_backbone/parsnp.xmfa backbone_alignment.fa
          else
            cp "${representatives_fasta}" backbone_alignment.fa
          fi
        else
          echo "Parsnp failed or produced no tree; falling back to FastTree."
          cp "${representatives_fasta}" backbone_alignment.fa
          if command -v fasttree >/dev/null 2>&1; then
            fasttree ${ft_opts} backbone_alignment.fa > backbone.treefile || status="FALLBACK"
          else
            status="FALLBACK"
          fi
        fi
      fi
    elif [ "${method}" = "fasttree" ]; then
      echo "Running FastTree backbone..."
      cp "${representatives_fasta}" backbone_alignment.fa
      if command -v fasttree >/dev/null 2>&1; then
        fasttree ${ft_opts} backbone_alignment.fa > backbone.treefile || status="FALLBACK"
      else
        status="FALLBACK"
      fi
    else
      echo "Unknown method '${method}'; using FastTree."
      cp "${representatives_fasta}" backbone_alignment.fa
      if command -v fasttree >/dev/null 2>&1; then
        fasttree ${ft_opts} backbone_alignment.fa > backbone.treefile || status="FALLBACK"
      else
        status="FALLBACK"
      fi
    fi
  fi

  # Ensure mandatory outputs exist
  [ -s backbone.treefile ] || awk '/^>/{gsub(/^>/,""); n=n (n? ",":"") \$0} END{print "(" n ");"}' "${representatives_fasta}" > backbone.treefile
  [ -s backbone_alignment.fa ] || cp "${representatives_fasta}" backbone_alignment.fa

  # Report
  {
    echo "BACKBONE TREE CONSTRUCTION REPORT"
    echo "================================="
    echo "Date (UTC): \$(date -u +%FT%TZ)"
    echo "Method: ${method}"
    echo "Representatives: \$rep_count"
    echo "Tree file: backbone.treefile"
    echo "Alignment file: backbone_alignment.fa"
    echo "Status: \$status"
    if [ -f "${mash_distances}" ]; then
      echo "Mash distances (optional input) present: yes"
    else
      echo "Mash distances (optional input) present: no"
    fi
  } > backbone_report.txt

  # Versions
  {
    echo "\"\${task.process}\":"
    echo "  parsnp: \$( (parsnp --version 2>/dev/null || echo 'N/A') | sed 's/^/    /')"
    echo "  fasttree: \$( (fasttree -expert 2>&1 | head -1 || echo 'N/A') | sed 's/^/    /')"
  } > versions.yml
  """
}
