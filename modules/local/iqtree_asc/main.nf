#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IQTREE_ASC -- final per-cluster ML tree on the Gubbins-filtered SNP alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    WHAT CHANGED
    1. `-nt AUTO` -> `-T ${task.cpus}`. Measured on 40 taxa x 3000 sites (22-core host),
       AUTO cost ~57 s per task purely to probe thread counts, and it logged
       "auto-detect threads (22 CPU cores detected)" -- it sizes to the HOST core count
       regardless of task.cpus, so with maxForks=2 both tasks tried to take the machine.
    2. `-bb 1000 -alrt 1000` are now gated behind params.iqtree_support (default false).
       They were unconditional. Measured on 50 taxa x 2000 SNP sites, -m GTR+G -T 4:
           -bb 1000 -alrt 1000  (old, always on)   151.0 s
           -alrt 1000 only                         119.2 s
           no support values                        91.1 s
       UFBoot+SH-aLRT is a 1.66x multiplier on the final per-cluster tree; across
       40-52 clusters that is ~1.0 h vs ~1.7 h of IQ-TREE time. SH-aLRT alone saves
       little (119.2 s of the 151.0 s), so this is effectively both-or-neither -- hence
       one boolean rather than two. Branch support is a real scientific output, so OFF
       by default is a DEFAULT CHANGE; set --iqtree_support true for publication runs.
    3. CORRECTNESS DEFECT (confirmed, not hypothetical): the `($names);` star-tree
       fallback is gone. GTR+ASC HARD-FAILS on any invariant column --
           "ERROR: Invalid use of +ASC because of 255 invariant sites in the alignment"
       -- and no treefile is written. The old module caught that and wrote `($names);`,
       a star tree with ZERO internal nodes, to ${cluster_id}.final.treefile, where it
       is indistinguishable from a successful result and flows on into
       SUMMARIZE_CLUSTER_PHYLOGENY and GRAFT_TREES.
       This fires on nominally SNP-only input: a 30-taxon 1500-site "SNP" alignment
       still contained 242 invariant columns, because in clonal genomes many sites are
       constant across the SAMPLED taxa while polymorphic in the wider population.
       Gubbins' filtered_polymorphic_sites.fasta is therefore very likely to trigger it
       for small clonal clusters.
       Measured remedies (30 taxa, 1500 sites, 242 invariant):
           GTR+ASC as-is (old)                 FAILS, no treefile
           GTR+G (drop +ASC)                   40.4 s, valid tree, but ASC-uncorrected
                                               branch lengths
           GTR+ASC on variable sites only      14.9 s, valid tree -- 2.7x FASTER than
                                               dropping ASC, and statistically correct
       bin/asc_preflight.py therefore DETECTS constant columns before IQ-TREE runs.
       Validated end to end against IQ-TREE 2.4.0 on the 30x1500/611-constant case,
       all three strategies producing a fully resolved 27-internal-node tree where the
       old code produced a 0-node star:
           varsites (DEFAULT)  strip constant cols, KEEP +ASC, no -fconst   4.89 s
           fconst              strip constant cols, DROP +ASC, -fconst      6.47 s
           drop_asc            alignment unchanged, DROP +ASC              6.63 s
       IMPORTANT: +ASC and -fconst are MUTUALLY EXCLUSIVE. -fconst reconstitutes the
       constant sites, so IQ-TREE re-raises "Invalid use of +ASC" and writes no tree
       (verified: fails in 0.01 s). The preflight therefore never emits both.
       'varsites' is the default because it is the only strategy that keeps the
       ascertainment-bias correction, and it is also the fastest. The other two change
       branch lengths (no ASC correction) and are gated behind
       params.iqtree_asc_fallback. Any other IQ-TREE failure now fails the task loudly.
       The <3-sequence case still writes a degenerate tree, because with 1-2 tips there
       IS no topology to infer; it is flagged DEGENERATE_TREE=1 in the preflight file.
*/

process IQTREE_ASC {
  tag "cluster_${cluster_id}"
  label 'process_high'
  container "quay.io/biocontainers/iqtree:2.2.6--h21ec9f0_0"

  publishDir "${params.outdir}/Clusters/cluster_${cluster_id}",
             mode: params.publish_dir_mode,
             pattern: "*.{treefile,iqtree,asc_preflight.txt}"

  input:
    tuple val(cluster_id), path(filtered_snps), val(representative_id)

  output:
    tuple val(cluster_id), path("${cluster_id}.final.treefile"), val(representative_id), emit: final_tree
    tuple val(cluster_id), path("${cluster_id}.final.iqtree"), emit: log
    tuple val(cluster_id), path("${cluster_id}.asc_preflight.txt"), emit: asc_preflight
    path "versions.yml", emit: versions

  when:
    task.ext.when == null || task.ext.when

  script:
    def args     = (task.ext.args ?: '').toString().trim()
    def model    = (params.iqtree_asc_model ?: 'GTR+ASC').toString().trim()
    // User-supplied constant-site counts still win over the preflight's own counts.
    def fconst   = (params.iqtree_fconst ?: '').toString().trim()
    // Branch support OFF by default: measured cost on this workload, see notes.
    def support  = (params.iqtree_support == null ? false : params.iqtree_support.toString().toLowerCase() in ['true','1','yes'])
    def ufboot   = (params.iqtree_ufboot ?: 1000) as int
    def alrt     = (params.iqtree_alrt ?: 1000) as int
    // varsites = strip constant columns and KEEP +ASC (default; ASC-correct, fastest)
    // fconst   = strip constant columns, drop +ASC, pass constant counts via -fconst
    // drop_asc = leave alignment alone, drop +ASC  (branch lengths not ASC-corrected)
    def asc_fb   = (params.iqtree_asc_fallback ?: 'varsites').toString().trim()
    """
    set -euo pipefail

    echo "Building final ML tree for cluster ${cluster_id} (representative ${representative_id})"

    IQTREE=\$(command -v iqtree2 || command -v iqtree || true)
    if [ -z "\$IQTREE" ]; then
      echo "ERROR: neither iqtree2 nor iqtree is on PATH" >&2
      exit 1
    fi

    if [ ! -s "${filtered_snps}" ]; then
      echo "ERROR: Gubbins filtered alignment for cluster ${cluster_id} is missing or empty: ${filtered_snps}" >&2
      exit 1
    fi

    seq_count=\$(grep -c "^>" ${filtered_snps} || echo 0)
    if [ "\$seq_count" -lt 3 ]; then
      # With 1-2 tips there is genuinely no topology to infer. This is the ONLY
      # remaining case that writes a degenerate Newick string, and it is recorded
      # in the preflight file so the summary can distinguish it from a real tree.
      echo "WARNING: cluster ${cluster_id} has \$seq_count sequences; <3 tips means no topology exists."
      names=\$(grep "^>" ${filtered_snps} | sed 's/^>//' | tr '\\n' ',' | sed 's/,\$//')
      echo "(\$names);" > ${cluster_id}.final.treefile
      : > ${cluster_id}.final.iqtree
      printf 'N_TAXA=%s\\nASC_ACTION=skipped_too_few_taxa\\nDEGENERATE_TREE=1\\n' "\$seq_count" \\
        > ${cluster_id}.asc_preflight.txt
      printf '"%s":\\n    iqtree: %s\\n' "${task.process}" "\$("\$IQTREE" --version 2>&1 | head -1)" > versions.yml
      exit 0
    fi

    # ---- ASC preflight -------------------------------------------------------
    # +ASC aborts on constant columns and writes NO tree. Gubbins'
    # filtered_polymorphic_sites.fasta contains them routinely for clonal clusters
    # (a 30-taxon 1500-site SNP alignment had 611), so decide BEFORE running IQ-TREE
    # rather than catching the failure and inventing a topology.
    asc_preflight.py \\
      --input "${filtered_snps}" \\
      --out "${cluster_id}.asc_preflight.txt" \\
      --model "${model}" \\
      --strategy "${asc_fb}" \\
      --stripped-output "${cluster_id}.asc_stripped.fasta"

    # shellcheck disable=SC1090
    . ./${cluster_id}.asc_preflight.txt
    echo "ASC decision: action=\$ASC_ACTION model=\$IQ_MODEL alignment=\$IQ_ALIGNMENT fconst='\$IQ_FCONST'"

    EXTRA="${args}"
    # An explicit params.iqtree_fconst overrides the preflight's computed counts.
    if [ -n "${fconst}" ]; then
      EXTRA="\${EXTRA} -fconst ${fconst}"
    elif [ -n "\$IQ_FCONST" ]; then
      EXTRA="\${EXTRA} -fconst \$IQ_FCONST"
    fi

    # Branch support is opt-in. UFBoot + SH-aLRT dominate per-cluster cost on these
    # alignments; see ALIGNMENT_TREES_NOTES.md for the measurement.
    if ${support}; then
      EXTRA="\${EXTRA} -bb ${ufboot} -alrt ${alrt}"
      echo "Branch support ENABLED: -bb ${ufboot} -alrt ${alrt}"
    else
      echo "Branch support DISABLED (params.iqtree_support=false); no -bb / -alrt"
    fi

    # -T \${task.cpus}, never -nt AUTO.
    "\$IQTREE" \\
      -s "\$IQ_ALIGNMENT" \\
      -st DNA \\
      -m "\$IQ_MODEL" \\
      -T ${task.cpus} \\
      --prefix ${cluster_id}.final \\
      \${EXTRA}

    # No star-tree fallback: if IQ-TREE fails, `set -e` above fails the task.
    if [ ! -s "${cluster_id}.final.treefile" ]; then
      echo "ERROR: IQ-TREE exited 0 but produced no ${cluster_id}.final.treefile" >&2
      exit 1
    fi
    [ -f "${cluster_id}.final.iqtree" ] || : > ${cluster_id}.final.iqtree

    printf '"%s":\\n    iqtree: %s\\n' "${task.process}" "\$("\$IQTREE" --version 2>&1 | head -1)" > versions.yml
    echo "Final ML tree completed for cluster ${cluster_id}"
    """
}
