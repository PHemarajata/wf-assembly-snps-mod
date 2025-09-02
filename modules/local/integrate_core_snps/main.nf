process INTEGRATE_CORE_SNPS {
    tag "global_integration"
    label 'process_medium'
    container "quay.io/biocontainers/python:3.9--1"
    
    publishDir "${params.outdir}/Integrated_Results", mode: params.publish_dir_mode, pattern: "*.{fa,tsv,txt}"

    input:
    path core_snp_files
    path snp_position_files
    path clusters_file

    output:
    path "integrated_core_snps.fa", emit: integrated_alignment
    path "integrated_snp_positions.tsv", emit: integrated_positions
    path "core_snp_summary.txt", emit: summary
    path "sample_cluster_mapping.tsv", emit: sample_mapping
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Install compatible versions of numpy and pandas
    pip install --upgrade numpy>=1.15.4
    pip install pandas biopython

    echo "Starting core SNP integration..."
    echo "Input files received:"
    ls -la
    echo ""
    echo "Core SNP files:"
    ls -la *_core_snps.fa 2>/dev/null || echo "No core SNP files found"
    echo ""
    echo "SNP position files:"
    ls -la *_snp_positions.tsv 2>/dev/null || echo "No SNP position files found"
    echo ""

    python3 << 'EOF'
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict
import glob

def integrate_core_snps():
    # Integrate core SNPs from all clusters into a global alignment
    
    print("Starting core SNP integration process...")
    
    # Read cluster assignments
    try:
        clusters_df = pd.read_csv("${clusters_file}", sep='\\t')
        print(f"Read cluster assignments: {len(clusters_df)} entries")
        sample_to_cluster = dict(zip(clusters_df['sample_id'], clusters_df['cluster_id']))
    except Exception as e:
        print(f"Error reading cluster file: {e}")
        clusters_df = pd.DataFrame()
        sample_to_cluster = {}
    
    # Collect all core SNP sequences
    all_sequences = {}
    cluster_snp_counts = {}
    
    # Process each core SNP file
    core_snp_files = glob.glob("*_core_snps.fa")
    print(f"Found {len(core_snp_files)} core SNP files: {core_snp_files}")
    
    if not core_snp_files:
        print("WARNING: No core SNP files found!")
        print("This usually means:")
        print("1. EXTRACT_CORE_SNPS processes failed")
        print("2. No variable sites found in cluster alignments")
        print("3. All cluster alignments were empty")
    
    for snp_file in core_snp_files:
        cluster_id = snp_file.replace("_core_snps.fa", "")
        snp_count = 0
        
        print(f"Processing {snp_file} for cluster {cluster_id}")
        
        # Check file size
        file_size = os.path.getsize(snp_file)
        print(f"  File size: {file_size} bytes")
        
        if file_size == 0:
            print(f"  WARNING: {snp_file} is empty")
            cluster_snp_counts[cluster_id] = 0
            continue
        
        try:
            sequences_in_file = 0
            for record in SeqIO.parse(snp_file, "fasta"):
                sample_id = record.id
                snp_sequence = str(record.seq)
                sequences_in_file += 1
                
                if sample_id not in all_sequences:
                    all_sequences[sample_id] = []
                
                all_sequences[sample_id].append({
                    'cluster': cluster_id,
                    'sequence': snp_sequence,
                    'length': len(snp_sequence)
                })
                snp_count = len(snp_sequence)
            
            cluster_snp_counts[cluster_id] = snp_count
            print(f"  Processed {sequences_in_file} sequences, {snp_count} SNPs per sequence")
            
        except Exception as e:
            print(f"  ERROR processing {snp_file}: {e}")
            cluster_snp_counts[cluster_id] = 0
    
    print(f"\\nTotal samples with SNP data: {len(all_sequences)}")
    print(f"Cluster SNP counts: {cluster_snp_counts}")
    
    # Create integrated alignment by concatenating SNPs from all clusters
    integrated_sequences = {}
    
    for sample_id, cluster_data in all_sequences.items():
        # Sort by cluster to ensure consistent order
        cluster_data.sort(key=lambda x: x['cluster'])
        
        # Concatenate sequences from all clusters
        concatenated_seq = ''.join([data['sequence'] for data in cluster_data])
        integrated_sequences[sample_id] = concatenated_seq
        print(f"Sample {sample_id}: {len(concatenated_seq)} total SNPs")
    
    # Write integrated core SNPs alignment
    sequences_written = 0
    with open("integrated_core_snps.fa", 'w') as f:
        for sample_id, sequence in integrated_sequences.items():
            if sequence:  # Only write non-empty sequences
                f.write(f">{sample_id}\\n{sequence}\\n")
                sequences_written += 1
    
    print(f"Wrote {sequences_written} sequences to integrated_core_snps.fa")
    
    # Check output file
    output_size = os.path.getsize("integrated_core_snps.fa")
    print(f"Output file size: {output_size} bytes")
    
    # Integrate SNP position information
    all_positions = []
    position_offset = 0
    
    snp_position_files = glob.glob("*_snp_positions.tsv")
    print(f"Found {len(snp_position_files)} SNP position files")
    
    for pos_file in snp_position_files:
        try:
            cluster_positions = pd.read_csv(pos_file, sep='\\t')
            if not cluster_positions.empty:
                # Adjust positions by offset
                cluster_positions['global_position'] = cluster_positions['position'] + position_offset
                cluster_positions['original_position'] = cluster_positions['position']
                all_positions.append(cluster_positions)
                
                # Update offset for next cluster
                position_offset += cluster_positions['position'].max()
                print(f"  Processed {len(cluster_positions)} positions from {pos_file}")
            else:
                print(f"  WARNING: {pos_file} is empty")
                
        except Exception as e:
            print(f"  ERROR processing {pos_file}: {e}")
    
    # Combine all position data
    if all_positions:
        integrated_positions = pd.concat(all_positions, ignore_index=True)
        integrated_positions.to_csv("integrated_snp_positions.tsv", sep='\\t', index=False)
        print(f"Integrated {len(integrated_positions)} SNP positions")
    else:
        # Create empty file
        pd.DataFrame(columns=['position', 'ref_base', 'alt_bases', 'cluster_id', 'global_position', 'original_position']).to_csv("integrated_snp_positions.tsv", sep='\\t', index=False)
        print("No SNP positions to integrate - created empty file")
    
    # Create sample-cluster mapping
    sample_mapping = []
    for sample_id in integrated_sequences.keys():
        cluster_id = sample_to_cluster.get(sample_id, 'unknown')
        sample_mapping.append({
            'sample_id': sample_id,
            'cluster_id': cluster_id,
            'total_snps': len(integrated_sequences[sample_id])
        })
    
    mapping_df = pd.DataFrame(sample_mapping)
    mapping_df.to_csv("sample_cluster_mapping.tsv", sep='\\t', index=False)
    print(f"Created sample mapping for {len(mapping_df)} samples")
    
    # Create summary report
    with open("core_snp_summary.txt", 'w') as f:
        f.write("INTEGRATED CORE SNP ANALYSIS SUMMARY\\n")
        f.write("=" * 50 + "\\n\\n")
        
        f.write(f"Total samples: {len(integrated_sequences)}\\n")
        f.write(f"Total clusters processed: {len(cluster_snp_counts)}\\n")
        f.write(f"Total integrated SNP positions: {sum(cluster_snp_counts.values())}\\n\\n")
        
        f.write("SNPs per cluster:\\n")
        f.write("-" * 20 + "\\n")
        for cluster_id, count in sorted(cluster_snp_counts.items()):
            f.write(f"Cluster {cluster_id}: {count} SNPs\\n")
        
        f.write(f"\\nIntegrated alignment length: {max([len(seq) for seq in integrated_sequences.values()]) if integrated_sequences else 0}\\n")
        
        if integrated_sequences:
            avg_snps = sum([len(seq) for seq in integrated_sequences.values()]) / len(integrated_sequences)
            f.write(f"Average SNPs per sample: {avg_snps:.2f}\\n")
        else:
            f.write("\\nWARNING: No sequences integrated!\\n")
            f.write("Possible causes:\\n")
            f.write("- All cluster alignments were empty\\n")
            f.write("- No variable sites found in any cluster\\n")
            f.write("- EXTRACT_CORE_SNPS processes failed\\n")
    
    print(f"Integration complete: {len(integrated_sequences)} samples, {sum(cluster_snp_counts.values())} total SNP positions")
    
    if len(integrated_sequences) == 0:
        print("\\nWARNING: No sequences were integrated!")
        print("This will cause downstream processes to fail.")
        print("Check the EXTRACT_CORE_SNPS and cluster alignment processes.")

# Run integration
integrate_core_snps()
EOF

    echo ""
    echo "Integration process completed. Final output files:"
    ls -la integrated_core_snps.fa sample_cluster_mapping.tsv integrated_snp_positions.tsv core_snp_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}