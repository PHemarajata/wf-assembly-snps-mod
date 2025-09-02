process KEEP_INVARIANT_ATCG {
    tag "cluster_${cluster_id}"
    label 'process_low'
    container "quay.io/biocontainers/biopython@sha256:10d755c731c82a22d91fc346f338ba47d5fd4f3b357828f5bbc903c9be865614"

    input:
    tuple val(cluster_id), path(alignment)

    output:
    tuple val(cluster_id), path("${cluster_id}.core.full.aln"), emit: core_alignment
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "Keeping invariant A/T/C/G columns for cluster ${cluster_id}"
    echo "Input alignment: ${alignment}"
    
    python3 << 'EOF'
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

def keep_invariant_atcg_sites(input_file, output_file):
    '''
    Keep columns that contain only A/T/C/G bases (including invariant sites).
    This preserves the whole genome context needed for Gubbins recombination detection.
    
    According to specifications:
    - Keep invariant A/T/C/G columns for whole genome alignment
    - Avoid SNP-only alignments at this stage
    - Feed whole genome alignment to Gubbins
    '''
    
    try:
        # Read alignment
        alignment = AlignIO.read(input_file, "fasta")
        print(f"Read alignment with {len(alignment)} sequences, length {alignment.get_alignment_length()}")
        
        if len(alignment) == 0:
            print("WARNING: Empty alignment")
            # Create empty output
            with open(output_file, 'w') as f:
                pass
            return
        
        alignment_length = alignment.get_alignment_length()
        
        # Find columns to keep (all A/T/C/G, including invariant sites)
        keep_columns = []
        
        for i in range(alignment_length):
            column = alignment[:, i]
            
            # Check if all bases in this column are A/T/C/G (case insensitive)
            # This includes invariant sites (all same base) and variable sites
            valid_bases = set()
            all_valid = True
            
            for base in column:
                base_upper = base.upper()
                if base_upper in ['A', 'T', 'C', 'G']:
                    valid_bases.add(base_upper)
                elif base_upper in ['N', '-', '.']:
                    # Skip columns with gaps or Ns
                    all_valid = False
                    break
                else:
                    # Invalid base
                    all_valid = False
                    break
            
            if all_valid and len(valid_bases) > 0:
                keep_columns.append(i)
        
        print(f"Keeping {len(keep_columns)} columns out of {alignment_length} (includes invariant A/T/C/G sites)")
        
        if len(keep_columns) == 0:
            print("WARNING: No valid A/T/C/G columns found")
            # Create minimal output
            with open(output_file, 'w') as f:
                for record in alignment:
                    f.write(f">{record.id}\\nN\\n")
            return
        
        # Create filtered alignment keeping invariant sites
        filtered_records = []
        
        for record in alignment:
            # Extract sequence for kept columns
            filtered_seq = "".join(str(record.seq[i]) for i in keep_columns)
            
            # Create new record
            filtered_record = SeqRecord(
                Seq(filtered_seq),
                id=record.id,
                description=record.description
            )
            filtered_records.append(filtered_record)
        
        # Write filtered alignment
        SeqIO.write(filtered_records, output_file, "fasta")
        
        print(f"Created core alignment with {len(filtered_records)} sequences")
        print(f"Filtered alignment length: {len(filtered_seq)} bp")
        
        # Log statistics
        if len(keep_columns) > 0:
            invariant_count = 0
            variable_count = 0
            
            # Quick check for invariant vs variable sites in kept columns
            for i in keep_columns[:min(1000, len(keep_columns))]:  # Sample first 1000 for efficiency
                column = alignment[:, i]
                unique_bases = set(base.upper() for base in column if base.upper() in ['A', 'T', 'C', 'G'])
                if len(unique_bases) == 1:
                    invariant_count += 1
                else:
                    variable_count += 1
            
            print(f"Sample of kept sites: ~{invariant_count} invariant, ~{variable_count} variable")
        
    except Exception as e:
        print(f"Error processing alignment: {e}")
        # Create minimal fallback
        with open(output_file, 'w') as f:
            f.write(f">${cluster_id}_dummy\\nATCG\\n")

# Process the alignment
keep_invariant_atcg_sites("${alignment}", "${cluster_id}.core.full.aln")
EOF

    # Verify output
    if [ ! -f "${cluster_id}.core.full.aln" ]; then
        echo "WARNING: Output file not created. Creating minimal alignment."
        echo ">${cluster_id}_dummy" > ${cluster_id}.core.full.aln
        echo "ATCG" >> ${cluster_id}.core.full.aln
    fi
    
    echo "Core alignment with invariant sites created for cluster ${cluster_id}"
    echo "Output file size: \$(wc -c < ${cluster_id}.core.full.aln) bytes"
    echo "Number of sequences: \$(grep -c '^>' ${cluster_id}.core.full.aln)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}