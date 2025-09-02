process COLLECT_REPRESENTATIVES {
    tag "collect_reps"
    label 'process_low'
    container "ubuntu:jammy"

    input:
    path representative_files
    path cluster_info

    output:
    path "representatives.fa", emit: representatives_fasta
    path "cluster_representatives.tsv", emit: representatives_mapping
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "Collecting cluster representatives into single FASTA file"
    
    # Initialize output files
    > representatives.fa
    echo -e "cluster_id\\trepresentative_id\\tfile_path" > cluster_representatives.tsv
    
    # Count input files
    rep_count=\$(ls *.fa 2>/dev/null | wc -l)
    echo "Found \$rep_count representative files"
    
    if [ \$rep_count -eq 0 ]; then
        echo "WARNING: No representative files found"
        echo "Creating empty output files"
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ubuntu: \$(awk -F ' ' '{print \$2,\$3}' /etc/issue | tr -d '\n')
END_VERSIONS
        
        exit 0
    fi
    
    # Process each representative file
    for rep_file in *.fa; do
        if [ "\$rep_file" = "*.fa" ] || [ "\$rep_file" = "representatives.fa" ]; then
            continue  # No files matched or skip output file
        fi

        echo "Processing representative file: \$rep_file"

        # Extract representative ID from filename (remove .fa extension)
        rep_id=\$(basename \$rep_file .fa)

        # Try to determine cluster ID
        # This assumes the file naming convention includes cluster information
        cluster_id=""

        # Method 1: Look for cluster info in filename
        if [[ "\$rep_file" == *"cluster_"* ]]; then
            cluster_id=\$(echo \$rep_file | sed 's/.*cluster_\\([^_]*\\).*/\\1/')
        fi

        # Method 2: Use representative ID as cluster ID if not found
        if [ -z "\$cluster_id" ]; then
            cluster_id=\$rep_id
        fi

        echo "Representative: \$rep_id, Cluster: \$cluster_id"

        # Check if file has content
        if [ -s "\$rep_file" ]; then
            # Add to combined FASTA
            cat \$rep_file >> representatives.fa

            # Add to mapping
            echo -e "\$cluster_id\t\$rep_id\t\$rep_file" >> cluster_representatives.tsv
        else
            echo "WARNING: Empty representative file: \$rep_file"
        fi
    done
    
    # Verify output
    final_count=\$(grep -c "^>" representatives.fa 2>/dev/null || echo "0")
    echo "Collected \$final_count representatives in final FASTA"
    
    if [ "\$final_count" -eq 0 ]; then
        echo "WARNING: No sequences collected. Creating minimal representative."
        echo ">dummy_representative" > representatives.fa
        echo "ATCG" >> representatives.fa
        echo -e "dummy\\tdummy_representative\\tdummy.fa" >> cluster_representatives.tsv
    fi
    
    # Summary
    echo "Representative collection summary:"
    echo "- Input files: \$rep_count"
    echo "- Output sequences: \$final_count"
    echo "- Output file: representatives.fa"
    echo "- Mapping file: cluster_representatives.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ubuntu: \$(awk -F ' ' '{print \$2,\$3}' /etc/issue | tr -d '\n')
END_VERSIONS
    """
}