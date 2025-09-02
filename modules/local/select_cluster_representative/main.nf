process SELECT_CLUSTER_REPRESENTATIVE {
    tag "cluster_${cluster_id}"
    label 'process_low'
    container "quay.io/biocontainers/python:3.9--1"

    input:
    tuple val(cluster_id), val(sample_ids), path(assemblies)
    path mash_distances

    output:
        tuple val(cluster_id), path("representative_id.txt"), path("*.fa"), emit: representative
        path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pip install pandas numpy

python3 << 'EOF'
import pandas as pd
import numpy as np
import os
import shutil

def normalize_name(name):
    return os.path.splitext(os.path.basename(name))[0]

def select_medoid(cluster_id, sample_ids, mash_distances_file):
    print(f"Selecting representative for cluster {cluster_id}")
    print(f"Samples in cluster: {sample_ids}")

    try:
        distances_df = pd.read_csv(mash_distances_file, sep='\t', index_col=0)
        distances_df.index = [normalize_name(x) for x in distances_df.index]
        distances_df.columns = [normalize_name(x) for x in distances_df.columns]
        print(f"Loaded distance matrix with shape: {distances_df.shape}")
    except Exception as e:
        print(f"Error reading distance matrix: {e}")
        return normalize_name(sample_ids[0])

    # Normalize all sample_ids
    normalized_ids = [normalize_name(s) for s in sample_ids]
    cluster_samples = [s for s in normalized_ids if s in distances_df.index and s in distances_df.columns]

    if len(cluster_samples) < len(normalized_ids):
        print(f"Warning: Only {len(cluster_samples)} of {len(normalized_ids)} samples found in distance matrix")
        if not cluster_samples:
            print("No samples found in distance matrix, selecting first sample")
            return normalized_ids[0]

    if len(cluster_samples) == 1:
        return cluster_samples[0]

    cluster_distances = distances_df.loc[cluster_samples, cluster_samples]
    distance_sums = cluster_distances.sum(axis=1)
    medoid = distance_sums.idxmin()

    print(f"Selected medoid: {medoid} (sum of distances: {distance_sums[medoid]:.6f})")
    return medoid

cluster_id = "${cluster_id}"
sample_ids = "${sample_ids}".strip('[]').replace(' ', '').split(',')
mash_distances_file = "${mash_distances}"

representative_id = select_medoid(cluster_id, sample_ids, mash_distances_file)
print(f"Representative for cluster {cluster_id}: {representative_id}")

# Write representative_id to file for Nextflow output
with open("representative_id.txt", "w") as f:
    f.write(representative_id + "\\n")

assembly_files = [f for f in os.listdir('.') if f.endswith('.fa') or f.endswith('.fasta') or f.endswith('.fna')]
representative_file = None

for assembly_file in assembly_files:
    base_name = normalize_name(assembly_file)
    if base_name == representative_id:
        representative_file = assembly_file
        break

if representative_file:
    shutil.copy(representative_file, f"{representative_id}.fa")
    print(f"Copied {representative_file} to {representative_id}.fa")
else:
    print(f"Warning: Could not find assembly file for representative {representative_id}")
    with open(f"{representative_id}.fa", 'w') as f:
        f.write(f">{representative_id}\\nN\\n")
EOF

cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
END_VERSIONS
    """
}