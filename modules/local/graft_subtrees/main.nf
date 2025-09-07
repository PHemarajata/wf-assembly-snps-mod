process GRAFT_SUBTREES {
  tag "tree_grafting"
  label 'process_medium'
  container "quay.io/biocontainers/gotree:0.4.4--h9ee0642_0"

  publishDir "${params.outdir}/Final_Results",
             mode: params.publish_dir_mode,
             pattern: "*.{tre,txt,pdf,log,yml}"

  input:
    path backbone_tree
    path cluster_trees
    path cluster_representatives

  output:
    path "global_grafted.treefile", emit: grafted_tree
    path "grafting_report.txt",    emit: report
    path "grafting_log.txt",       emit: log
    path "versions.yml",           emit: versions

  when:
    task.ext.when == null || task.ext.when

  script:
  """
  set -Eeuo pipefail

  # -- Let Nextflow expand these paths once, then use shell vars everywhere --
  BACKBONE="${backbone_tree}"
  CLUSTER_TREES="${cluster_trees}"
  REPS="${cluster_representatives}"
  PROC_NAME="${task.process}"

  : > grafting_log.txt
  printf 'Starting tree grafting process\n' >> grafting_log.txt
  printf 'Backbone tree: %s\n' "\$BACKBONE" >> grafting_log.txt

  # Collect cluster trees, EXCLUDING the backbone filename
  shopt -s nullglob
  backbone_basename=\$(basename "\$BACKBONE")
  declare -a cluster_list=()

  # Accept either staged individual *.treefile or a staged directory of trees
  for f in ./*.treefile "\$CLUSTER_TREES"/*.treefile; do
    [[ -e "\$f" ]] || continue
    [[ "\$(basename "\$f")" == "\$backbone_basename" ]] && continue
    # avoid duplicates
    already=false
    for g in "\${cluster_list[@]:-}"; do
      [[ "\$g" == "\$f" ]] && { already=true; break; }
    done
    \$already || cluster_list+=( "\$f" )
  done
  shopt -u nullglob

  cluster_count=\${#cluster_list[@]}
  printf 'Found %s cluster trees to graft (excluding backbone)\n' "\$cluster_count" >> grafting_log.txt

  if [[ "\$cluster_count" -eq 0 ]]; then
    printf 'WARNING: No cluster trees found for grafting\n' >> grafting_log.txt
    printf 'Copying backbone tree as final result\n' >> grafting_log.txt
    cp "\$BACKBONE" global_grafted.treefile

    {
      printf 'TREE GRAFTING REPORT\n'
      printf '===================\n'
      printf 'Status: No cluster trees to graft\n'
      printf 'Result: Backbone tree used as final tree\n'
    } > grafting_report.txt

    gotree_version=\$(gotree version 2>/dev/null | sed 's/^/    /' || echo '    unknown')
    printf '\"%s\":\n    %s\n' "\$PROC_NAME" "gotree: \${gotree_version}" > versions.yml
    exit 0
  fi

  # Work on a copy of the backbone
  cp "\$BACKBONE" current_tree.tre

  # --- Helpers: normalization + unquoting ---
  normalize() {
    # strip surrounding quotes, common extensions, spaces->underscores
    local s="\$1"
    s="\${s%\"}"; s="\${s#\"}"; s="\${s%'}"; s="\${s#'}"
    s="\${s%.fa}"; s="\${s%.fasta}"; s="\${s%.fa.gz}"; s="\${s%.fna}"
    s="\${s%.treefile}"
    s="\${s// /_}"
    printf '%s' "\$s"
  }
  strip_quotes() {
    local s="\$1"
    s="\${s%\"}"; s="\${s#\"}"; s="\${s%'}"; s="\${s#'}"
    printf '%s' "\$s"
  }

  # Build backbone label map: norm -> actual label (as printed by gotree)
  declare -A BBMAP
  mapfile -t bb_labels < <(gotree labels -i current_tree.tre || true)
  {
    printf 'Backbone labels (first 20):\n'
    idx=0
    for L in "\${bb_labels[@]}"; do
      printf '  %s\n' "\$L"
      idx=\$((idx+1)); [[ "\$idx" -ge 20 ]] && break
    done
  } >> grafting_log.txt
  for L in "\${bb_labels[@]}"; do
    n=\$(normalize "\$L")
    [[ -n "\${BBMAP[\$n]:-}" ]] || BBMAP["\$n"]="\$L"
  done

  # Representatives mapping: cluster_id <TAB> representative
  declare -A REPMAP
  if [[ -s "\$REPS" ]]; then
    while IFS=\$'\\t' read -r cid rep || [[ -n "\$cid" ]]; do
      [[ -z "\$cid" || -z "\$rep" ]] && continue
      REPMAP["\$cid"]="\$rep"
      REPMAP["\${cid}.final"]="\$rep"
    done < "\$REPS"
    printf 'Loaded representatives mapping from %s\n' "\$REPS" >> grafting_log.txt
  else
    printf 'No cluster representatives mapping provided; will infer from cluster tree labels.\n' >> grafting_log.txt
  fi

  grafted_count=0
  failed_count=0

  for cluster_tree in "\${cluster_list[@]}"; do
    [[ ! -s "\$cluster_tree" ]] && { printf 'Skipping empty file %s\n' "\$cluster_tree" >> grafting_log.txt; continue; }

    cluster_id=\$(basename "\$cluster_tree" .treefile)
    printf 'Processing cluster: %s\n' "\$cluster_id" >> grafting_log.txt

    # Prefer mapping file
    representative_id="\${REPMAP[\$cluster_id]:-}"

    # Index cluster labels by normalized form
    declare -A CLMAP
    mapfile -t cl_labels < <(gotree labels -i "\$cluster_tree" || true)
    for L in "\${cl_labels[@]}"; do
      n=\$(normalize "\$L")
      [[ -n "\${CLMAP[\$n]:-}" ]] || CLMAP["\$n"]="\$L"
    done

    if [[ -z "\$representative_id" ]]; then
      chosen=""
      for L in "\${cl_labels[@]}"; do
        n=\$(normalize "\$L")
        if [[ -n "\${BBMAP[\$n]:-}" ]]; then chosen="\$L"; break; fi
      done
      [[ -z "\$chosen" ]] && chosen="\${cl_labels[0]:-}"
      representative_id="\$chosen"
      printf 'Representative inferred from cluster labels: %s\n' "\$representative_id" >> grafting_log.txt
    else
      printf 'Representative from mapping: %s\n' "\$representative_id" >> grafting_log.txt
    fi

    # Find the backbone label to attach to
    rep_norm=\$(normalize "\$representative_id")
    bb_label="\${BBMAP[\$rep_norm]:-}"

    if [[ -z "\$bb_label" ]]; then
      for k in "\${!BBMAP[@]}"; do
        if [[ "\$k" == *"\$rep_norm"* || "\$rep_norm" == *"\$k"* ]]; then
          bb_label="\${BBMAP[\$k]}"
          printf 'Using fuzzy backbone match: [%s] for rep [%s]\n' "\$bb_label" "\$representative_id" >> grafting_log.txt
          break
        fi
      done
    fi

    if [[ -z "\$bb_label" ]]; then
      printf 'WARNING: No matching representative in backbone for cluster [%s] (rep=[%s]).\n' "\$cluster_id" "\$representative_id" >> grafting_log.txt
      printf 'Available backbone labels (first 10):\n' >> grafting_log.txt
      idx=0; for L in "\${bb_labels[@]}"; do printf '  %s\n' "\$L" >> grafting_log.txt; idx=\$((idx+1)); [[ "\$idx" -ge 10 ]] && break; done
      failed_count=\$((failed_count + 1))
      continue
    fi

    # Ensure representative exists in the cluster tree (by normalized form)
    cl_rep_label="\${CLMAP[\$rep_norm]:-}"
    if [[ -z "\$cl_rep_label" ]]; then
      for k in "\${!CLMAP[@]}"; do
        if [[ "\$k" == *"\$rep_norm"* || "\$rep_norm" == *"\$k"* ]]; then
          cl_rep_label="\${CLMAP[\$k]}"
          printf 'Using fuzzy cluster match: [%s] for rep [%s]\n' "\$cl_rep_label" "\$representative_id" >> grafting_log.txt
          break
        fi
      done
    fi

    if [[ -z "\$cl_rep_label" ]]; then
      printf 'WARNING: Representative [%s] not found in cluster tree [%s].\n' "\$representative_id" "\$cluster_tree" >> grafting_log.txt
      failed_count=\$((failed_count + 1))
      continue
    fi

    printf 'Representative resolved -> backbone:[%s]  cluster:[%s]\n' "\$bb_label" "\$cl_rep_label" >> grafting_log.txt

    # Unquote raw labels for gotree flags
    tip_label=\$(strip_quotes "\$bb_label")
    cl_rep_raw=\$(strip_quotes "\$cl_rep_label")

    # If strings differ, rename the cluster’s rep label to the (unquoted) backbone label
    cluster_for_graft="\$cluster_tree"
    if [[ "\$tip_label" != "\$cl_rep_raw" ]]; then
      printf '%s\\t%s\n' "\$cl_rep_raw" "\$tip_label" > "\${cluster_id}.rename.tsv"
      if gotree rename -i "\$cluster_tree" -m "\${cluster_id}.rename.tsv" -o "\${cluster_id}.renamed.tre" 2>> grafting_log.txt; then
        cluster_for_graft="\${cluster_id}.renamed.tre"
        printf 'Renamed cluster rep label to [%s] for grafting.\n' "\$tip_label" >> grafting_log.txt
      else
        printf 'WARNING: Could not rename cluster rep label; attempting graft with original names (may fail).\n' >> grafting_log.txt
      fi
    fi

    # Graft cluster tree onto backbone at the representative label
    if ! gotree graft -i current_tree.tre -c "\$cluster_for_graft" -l "\$tip_label" -o temp_grafted.tre 2>> grafting_log.txt; then
      printf 'WARNING: Grafting failed for cluster [%s] at label [%s].\n' "\$cluster_id" "\$tip_label" >> grafting_log.txt
      failed_count=\$((failed_count + 1))
      continue
    fi

    mv temp_grafted.tre current_tree.tre
    grafted_count=\$((grafted_count + 1))
    printf 'Successfully grafted cluster [%s] at [%s].\n' "\$cluster_id" "\$tip_label" >> grafting_log.txt
  done

  # Final result
  cp current_tree.tre global_grafted.treefile
  printf 'Grafting completed: %s successful, %s failed\n' "\$grafted_count" "\$failed_count" >> grafting_log.txt

  # Report
  {
    printf 'TREE GRAFTING REPORT\n'
    printf '===================\n'
    printf 'Backbone tree: %s\n' "\$BACKBONE"
    printf 'Total cluster trees: %s\n' "\$cluster_count"
    printf 'Successfully grafted: %s\n' "\$grafted_count"
    printf 'Failed grafting: %s\n' "\$failed_count"
    printf '\\n'
    if [[ "\$grafted_count" -gt 0 ]]; then
      printf 'Status: SUCCESS (partial or complete)\n'
      printf 'Final tree: global_grafted.treefile\n'
      final_leaves=\$(gotree labels -i global_grafted.treefile | wc -l 2>/dev/null || echo "unknown")
      printf 'Final tree leaves: %s\n' "\$final_leaves"
    else
      printf 'Status: FAILED (no successful grafts)\n'
      printf 'Result: Backbone tree used as fallback\n'
    fi
    printf '\\n'
    printf 'Method: gotree graft; labels normalized & quotes stripped; representative renamed in cluster trees when needed.\n'
  } > grafting_report.txt

  printf 'Tree grafting process completed\n' >> grafting_log.txt

  # Versions (single printf — no unclosed blocks)
  gotree_version=\$(gotree version 2>/dev/null | sed 's/^/    /' || echo '    unknown')
  printf '\"%s\":\n    %s\n' "\$PROC_NAME" "gotree: \${gotree_version}" > versions.yml
  """
}
