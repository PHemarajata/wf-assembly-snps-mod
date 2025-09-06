#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process EXTRACT_CORE_SNPS {
    tag "cluster_${cluster_id}"
    label 'process_low'
    container "quay.io/biocontainers/python:3.9--1"
    
    publishDir "${params.outdir}/Core_SNPs", mode: params.publish_dir_mode, pattern: "*.{fa,tsv}"

    input:
    tuple val(cluster_id), path(alignment)
    path clusters_file

    output:
    tuple val(cluster_id), path("${cluster_id}_core_snps.fa"), emit: core_snps
    tuple val(cluster_id), path("${cluster_id}_snp_positions.tsv"), emit: snp_positions
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Install compatible versions of numpy and pandas
    pip install --upgrade numpy>=1.15.4
    pip install pandas biopython

    echo "Starting core SNP extraction for cluster ${cluster_id}"
    echo "Input alignment file: ${alignment}"
    echo "Alignment file size: \$(wc -c < "${alignment}") bytes"
    echo "Number of sequences: \$(grep -c "^>" "${alignment}" 2>/dev/null || echo "0")"
    echo ""

    python3 << 'EOF'
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

def extract_core_snps(alignment_file, cluster_id, clusters_file):
    # Extract core SNP positions from cluster alignment
    
    print(f"Processing cluster {cluster_id}")
    print(f"Alignment file: {alignment_file}")
    
    # Check if alignment file exists and has content
    if not os.path.exists(alignment_file):
        print(f"ERROR: Alignment file {alignment_file} does not exist")
        # Create empty files
        with open(f"{cluster_id}_core_snps.fa", 'w') as f:
            pass
        with open(f"{cluster_id}_snp_positions.tsv", 'w') as f:
            f.write("position\\tref_base\\talt_bases\\tcluster_id\\n")
        return
    
    file_size = os.path.getsize(alignment_file)
    print(f"Alignment file size: {file_size} bytes")
    
    if file_size == 0:
        print(f"WARNING: Alignment file {alignment_file} is empty")
        # Create empty files
        with open(f"{cluster_id}_core_snps.fa", 'w') as f:
            pass
        with open(f"{cluster_id}_snp_positions.tsv", 'w') as f:
            f.write("position\\tref_base\\talt_bases\\tcluster_id\\n")
        return
    
    # Read cluster assignments to get sample names
    try:
        clusters_df = pd.read_csv(clusters_file, sep='\\t')
        cluster_samples = clusters_df[clusters_df['cluster_id'] == cluster_id]['sample_id'].tolist()
        print(f"Expected samples in cluster {cluster_id}: {cluster_samples}")
    except Exception as e:
        print(f"Warning: Could not read cluster assignments: {e}")
        cluster_samples = []
    
    # Read alignment
    sequences = {}
    try:
        for record in SeqIO.parse(alignment_file, "fasta"):
            sequences[record.id] = str(record.seq)
        print(f"Read {len(sequences)} sequences from alignment")
        
        if sequences:
            seq_names = list(sequences.keys())
            print(f"Sequence names: {seq_names}")
            alignment_length = len(sequences[seq_names[0]])
            print(f"Alignment length: {alignment_length} bp")
        
    except Exception as e:
        print(f"ERROR reading alignment file: {e}")
        sequences = {}
    
    if not sequences:
        print(f"WARNING: No sequences found in alignment for cluster {cluster_id}")
        # Create empty files
        with open(f"{cluster_id}_core_snps.fa", 'w') as f:
            pass
        with open(f"{cluster_id}_snp_positions.tsv", 'w') as f:
            f.write("position\\tref_base\\talt_bases\\tcluster_id\\n")
        return
    
    # Get alignment length
    seq_names = list(sequences.keys())
    alignment_length = len(sequences[seq_names[0]])
    
    # Find variable positions (SNPs)
    snp_positions = []
    core_snps = {name: [] for name in seq_names}
    
    print(f"Scanning {alignment_length} positions for variable sites...")
    
    for pos in range(alignment_length):
        bases_at_pos = set()
        for name in seq_names:
            base = sequences[name][pos]
            if base not in ['-', 'N', 'n']:  # Exclude gaps and Ns
                bases_at_pos.add(base.upper())
        
        # If more than one base type at this position, it's a SNP
        if len(bases_at_pos) > 1:
            snp_positions.append({
                'position': pos + 1,  # 1-based position
                'ref_base': sorted(bases_at_pos)[0],  # Use first alphabetically as ref
                'alt_bases': ','.join(sorted(bases_at_pos)[1:]),
                'cluster_id': cluster_id
            })
            
            # Add this position to core SNPs for each sequence
            for name in seq_names:
                base = sequences[name][pos]
                core_snps[name].append(base.upper() if base not in ['-', 'N', 'n'] else 'N')
    
    print(f"Found {len(snp_positions)} variable positions")
    
    # Write core SNPs FASTA
    sequences_written = 0
    with open(f"{cluster_id}_core_snps.fa", 'w') as f:
        for name in seq_names:
            if core_snps[name]:  # Only write if there are SNPs
                snp_seq = ''.join(core_snps[name])
                f.write(f">{name}\\n{snp_seq}\\n")
                sequences_written += 1
    
    print(f"Wrote {sequences_written} sequences with SNPs to {cluster_id}_core_snps.fa")
    
    # Write SNP positions table
    snp_df = pd.DataFrame(snp_positions)
    snp_df.to_csv(f"{cluster_id}_snp_positions.tsv", sep='\\t', index=False)
    
    print(f"Cluster {cluster_id}: Found {len(snp_positions)} core SNP positions")
    
    if len(snp_positions) == 0:
        print(f"WARNING: No variable sites found in cluster {cluster_id}")
        print("This could be because:")
        print("- All sequences in the cluster are identical")
        print("- The alignment contains only gaps or Ns")
        print("- The cluster has insufficient diversity")

# Run extraction
extract_core_snps("${alignment}", "${cluster_id}", "${clusters_file}")
EOF

    echo ""
    echo "Core SNP extraction completed for cluster ${cluster_id}"
    echo "Output files:"
    ls -la ${cluster_id}_core_snps.fa ${cluster_id}_snp_positions.tsv
    echo "Core SNPs file size: \$(wc -c < "${cluster_id}_core_snps.fa") bytes"
    echo "SNP positions file size: \$(wc -c < "${cluster_id}_snp_positions.tsv") bytes"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
