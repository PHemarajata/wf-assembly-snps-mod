/* -*- coding: utf-8 -*- */
nextflow.enable.dsl=2

process BUILD_BACKBONE_TREE {
  tag "backbone_tree"
  label 'process_high'

  // Override via params.backbone_container in your config if you’d like.
  // This biocontainer typically includes parsnp, harvesttools, fasttree.
  container "${params.backbone_container ?: 'quay.io/biocontainers/parsnp:1.7.4--hdcf5f25_2'}"

    // Ensure the backbone tree is published at the run outdir root so downstream
    // tools can find it as `backbone.treefile` (the python graft script expects
    // this path relative to the output directory).
    publishDir "${params.outdir}", mode: params.publish_dir_mode, pattern: "backbone.treefile"
    // The alignment the tree was built from was emitted but never published, so it
    // survived only in work/ and vanished with routine cleanup -- leaving no way to
    // check what the backbone actually came from. The report was likewise unpublished.
    publishDir "${params.outdir}", mode: params.publish_dir_mode, pattern: "backbone_alignment.*"
    publishDir "${params.outdir}", mode: params.publish_dir_mode, pattern: "backbone_report.txt"

  input:
  path representatives_fasta
  path reference_genome    // Pass a real FASTA to anchor Parsnp, or "NO_FILE" sentinel to fall back to head -n1

  output:
  path "backbone.treefile",       emit: backbone_tree
  // Glob, not a fixed .fa: parsnp emits XMFA (multi-block, 4,168 LCBs x N records
  // on the 2,802-genome run), which is NOT a FASTA alignment. Naming it .fa made
  // every consumer -- including a careful reader -- misparse it. The extension now
  // states what the file is.
  path "backbone_alignment.*",    emit: backbone_alignment
  path "backbone_report.txt",     emit: report
  path "versions.yml",            emit: versions

  script:
  """
  set -euo pipefail

  method="${params.backbone_method ?: 'parsnp'}"
  ft_opts="${params.backbone_fasttree_opts ?: '-nt -gtr'}"
  threads=${task.cpus}

  echo "== BUILD_BACKBONE_TREE =="
  echo "Method: \${method}"
  echo "CPUs: \${threads}"

  # Count REPRESENTATIVES, not FASTA records. Every contig of a representative
  # carries that representative's label (COLLECT_REPRESENTATIVES rewrites headers),
  # so `grep -c '^>'` counted contigs: 7,608 for 76 representatives on the
  # 2,802-genome run. That was not merely a cosmetic figure in the report -- it
  # feeds the "<3 representatives" guard below, which with draft assemblies
  # (~100 contigs each) could never fire, so the tiny-set path was unreachable.
  contig_count=\$(grep -c '^>' "${representatives_fasta}" || echo 0)
  rep_count=\$(grep '^>' "${representatives_fasta}" | sort -u | wc -l)
  echo "Representatives: \$rep_count (across \$contig_count contigs)"

  status="SUCCESS"

  # Handle tiny sets gracefully
    if [ "\$rep_count" -lt 3 ]; then
        echo "WARNING: <3 representatives found; emitting trivial tree."
        # Distinct labels only. Every contig carries its representative's label, so
        # using every header here would emit one tip per CONTIG with duplicates.
        # This path is only reachable now that rep_count counts representatives.
        printf '(%s);\\n' "\$(grep '^>' "${representatives_fasta}" | sed 's/^>//' | sort -u | paste -sd, -)" > backbone.treefile
        cp "${representatives_fasta}" backbone_alignment.fa
        status="TRIVIAL"
    else
        if [ "\$method" = "parsnp" ]; then
            echo "Running Parsnp backbone..."
            mkdir -p reps
            awk '/^>/{fn="reps/" substr(\$0,2) ".fa"}{print > fn}' "${representatives_fasta}"

            # Prefer a user-supplied global reference (e.g., K96243); fall back to the
            # first representative if no real reference file was passed in.
            ref=""
            if [ -s "${reference_genome}" ] && [ "\$(basename "${reference_genome}")" != "NO_FILE" ]; then
                cp "${reference_genome}" reps/__global_reference.fa
                ref="reps/__global_reference.fa"
                echo "Using user-supplied reference for backbone: \$(basename "${reference_genome}")"
            else
                ref=\$(ls reps/*.fa | head -n1 || true)
                echo "No global reference provided; using first representative as backbone ref: \$ref"
            fi

            if [ -z "\$ref" ]; then
                echo "ERROR: no reference could be resolved for parsnp and there is no \
valid fallback -- FastTree on representatives.fa would be FastTree on unaligned \
concatenated assemblies. Supply --ref, or check that representatives were produced." >&2
                exit 1
            else
                set +e
                parsnp --sequences reps --reference "\$ref" --output-dir parsnp_backbone --use-fasttree --threads \$threads --verbose
                rc=\$?
                set -e
                if [ "\$rc" -eq 0 ] && [ -s parsnp_backbone/parsnp.tree ]; then
                    cp parsnp_backbone/parsnp.tree backbone.treefile
                    if [ -s parsnp_backbone/parsnp.xmfa ]; then
                        # XMFA, not FASTA. Extension says so; bin/xmfa_to_fasta.py
                        # converts it for anything wanting a single-block alignment.
                        cp parsnp_backbone/parsnp.xmfa backbone_alignment.xmfa
                    else
                        echo "WARNING: parsnp produced a tree but no XMFA; the published \
backbone_alignment.fa is RAW CONCATENATED GENOMES, not the alignment the tree came from." >&2
                        cp "${representatives_fasta}" backbone_alignment.fa
                    fi
                else
                    echo "ERROR: parsnp exited \$rc or produced no tree. The previous \
fallback ran FastTree on representatives.fa -- raw concatenated assemblies, not an \
alignment -- and, failing that, emitted a STAR TREE as the backbone. Neither is a \
result. See parsnp_backbone/parsnpAligner.log." >&2
                    exit 1
                fi
            fi
        else
            # FastTree requires an ALIGNMENT. representatives_fasta is raw concatenated
            # assemblies of differing lengths, so feeding it here produces either an
            # error or a meaningless tree. Fail loudly rather than either.
            echo "ERROR: backbone_method='fasttree' would run FastTree on \
representatives.fa, which is raw concatenated assemblies and NOT an alignment. \
Use backbone_method='parsnp' (the default), which aligns first." >&2
            exit 1
        fi
    fi

  # A star tree used to be written here whenever backbone.treefile was missing --
  # the same silent-success pattern that produced star trees per cluster (see
  # AUDIT_REPORT.md A.2). A backbone with no internal structure is not a backbone.
  if [ ! -s backbone.treefile ]; then
      echo "ERROR: no backbone tree was produced. Refusing to emit a star tree as \
the backbone." >&2
      exit 1
  fi
  if ! ls backbone_alignment.* >/dev/null 2>&1; then
      echo "ERROR: no backbone alignment was produced." >&2
      exit 1
  fi

  # Report
  {
    echo "BACKBONE TREE CONSTRUCTION REPORT"
    echo "================================="
    echo "Date (UTC): \$(date -u +%FT%TZ)"
    echo "Method: \${method}"
    echo "Representatives: \$rep_count"
    echo "Tree file: backbone.treefile"
    echo "Alignment file: \$(ls backbone_alignment.* 2>/dev/null | head -1)"
    echo "Alignment format: \$(ls backbone_alignment.xmfa >/dev/null 2>&1 && echo 'XMFA (parsnp multi-block; NOT a single-block FASTA alignment)' || echo 'FASTA (raw concatenated representatives, NOT an alignment)')"
    echo "Status: \$status"
    echo "Mash distances input: not required for backbone tree construction"
  } > backbone_report.txt

  # Versions
  {
    echo "\"${task.process}\":"
    echo "  parsnp: \$( (parsnp --version 2>/dev/null || echo 'N/A') | sed 's/^/    /')"
    echo "  fasttree: \$( (fasttree -expert 2>&1 | head -1 || echo 'N/A') | sed 's/^/    /')"
  } > versions.yml
  """
}