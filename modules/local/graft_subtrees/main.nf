process GRAFT_SUBTREES {
    tag "tree_grafting"
    label 'process_medium'
    container "quay.io/biocontainers/gotree:0.4.4--h21ec9f0_0"
    
    publishDir "${params.outdir}/Final_Results", mode: params.publish_dir_mode, pattern: "*.{tre,txt,pdf}"

    input:
    path backbone_tree
    path cluster_trees
    path cluster_representatives

    output:
    path "global_grafted.treefile", emit: grafted_tree
    path "grafting_report.txt", emit: report
    path "grafting_log.txt", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "Starting tree grafting process" > grafting_log.txt
    echo "Backbone tree: ${backbone_tree}" >> grafting_log.txt
    
    # Count cluster trees
    cluster_count=\$(ls *.treefile 2>/dev/null | wc -l)
    echo "Found \$cluster_count cluster trees to graft" >> grafting_log.txt
    
    if [ \$cluster_count -eq 0 ]; then
        echo "WARNING: No cluster trees found for grafting" >> grafting_log.txt
        echo "Copying backbone tree as final result"
        cp ${backbone_tree} global_grafted.treefile
        
        echo "TREE GRAFTING REPORT" > grafting_report.txt
        echo "===================" >> grafting_report.txt
        echo "Status: No cluster trees to graft" >> grafting_report.txt
        echo "Result: Backbone tree used as final tree" >> grafting_report.txt
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gotree: \$(gotree version | sed 's/^/    /')
        END_VERSIONS
        
        exit 0
    fi
    
    # Start with backbone tree
    cp ${backbone_tree} current_tree.tre
    
    # Read cluster representatives mapping
    echo "Processing cluster representatives..." >> grafting_log.txt
    
    # Process each cluster tree for grafting
    grafted_count=0
    failed_count=0
    
    for cluster_tree in *.treefile; do
        if [ "\$cluster_tree" = "*.treefile" ]; then
            continue  # No files matched
        fi
        
        # Extract cluster ID from filename
        cluster_id=\$(basename \$cluster_tree .treefile)
        echo "Processing cluster: \$cluster_id" >> grafting_log.txt
        
        # Find representative for this cluster
        representative_id=""
        
        # Try to extract representative from cluster representatives file
        if [ -f "${cluster_representatives}" ]; then
            representative_id=\$(grep "^\$cluster_id" ${cluster_representatives} | cut -f2 2>/dev/null || echo "")
        fi
        
        # If not found, try to infer from tree file content
        if [ -z "\$representative_id" ]; then
            # Get first leaf name from cluster tree as representative
            representative_id=\$(gotree labels -i \$cluster_tree | head -n1 2>/dev/null || echo "")
        fi
        
        if [ -z "\$representative_id" ]; then
            echo "WARNING: Could not determine representative for cluster \$cluster_id" >> grafting_log.txt
            failed_count=\$((failed_count + 1))
            continue
        fi
        
        echo "Representative for \$cluster_id: \$representative_id" >> grafting_log.txt
        
        # Check if representative exists in current tree
        if ! gotree labels -i current_tree.tre | grep -q "^\$representative_id\$"; then
            echo "WARNING: Representative \$representative_id not found in backbone tree" >> grafting_log.txt
            failed_count=\$((failed_count + 1))
            continue
        fi
        
        # Extract subtree rooted at representative from cluster tree
        gotree subtree -i \$cluster_tree -l \$representative_id > \${cluster_id}_subtree.tre 2>/dev/null || {
            echo "WARNING: Could not extract subtree for cluster \$cluster_id" >> grafting_log.txt
            failed_count=\$((failed_count + 1))
            continue
        }
        
        # Graft subtree onto backbone at representative position
        gotree graft \\
            -i current_tree.tre \\
            -g \${cluster_id}_subtree.tre \\
            -l \$representative_id \\
            > temp_grafted.tre 2>/dev/null || {
            echo "WARNING: Grafting failed for cluster \$cluster_id" >> grafting_log.txt
            failed_count=\$((failed_count + 1))
            continue
        }
        
        # Update current tree
        mv temp_grafted.tre current_tree.tre
        grafted_count=\$((grafted_count + 1))
        echo "Successfully grafted cluster \$cluster_id" >> grafting_log.txt
    done
    
    # Final result
    cp current_tree.tre global_grafted.treefile
    
    echo "Grafting completed: \$grafted_count successful, \$failed_count failed" >> grafting_log.txt
    
    # Create detailed report
    echo "TREE GRAFTING REPORT" > grafting_report.txt
    echo "===================" >> grafting_report.txt
    echo "Backbone tree: ${backbone_tree}" >> grafting_report.txt
    echo "Total cluster trees: \$cluster_count" >> grafting_report.txt
    echo "Successfully grafted: \$grafted_count" >> grafting_report.txt
    echo "Failed grafting: \$failed_count" >> grafting_report.txt
    echo "" >> grafting_report.txt
    
    if [ \$grafted_count -gt 0 ]; then
        echo "Status: SUCCESS (partial or complete)" >> grafting_report.txt
        echo "Final tree: global_grafted.treefile" >> grafting_report.txt
        
        # Count final tree leaves
        final_leaves=\$(gotree labels -i global_grafted.treefile | wc -l 2>/dev/null || echo "unknown")
        echo "Final tree leaves: \$final_leaves" >> grafting_report.txt
    else
        echo "Status: FAILED (no successful grafts)" >> grafting_report.txt
        echo "Result: Backbone tree used as fallback" >> grafting_report.txt
    fi
    
    echo "" >> grafting_report.txt
    echo "Grafting method: gotree graft with subtree extraction" >> grafting_report.txt
    echo "Note: Each cluster forms a monophyletic group grafted at its representative" >> grafting_report.txt
    
    # Optional: Rescale branches if requested
    if [ "${params.rescale_grafted_branches}" = "true" ]; then
        echo "Rescaling grafted branches..." >> grafting_log.txt
        
        # This would require more sophisticated branch length adjustment
        # For now, just note that rescaling was requested
        echo "Branch rescaling: Requested but not implemented" >> grafting_report.txt
    fi
    
    echo "Tree grafting process completed"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gotree: \$(gotree version | sed 's/^/    /')
    END_VERSIONS
    """
}