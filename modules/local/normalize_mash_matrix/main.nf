process NORMALIZE_MASH_MATRIX {
    tag "normalize_mash_matrix"
    label 'process_low'
    container "python:3.9"

    input:
    path mash_matrix

    output:
    path "mash_distances.tsv", emit: distances
    path "versions.yml", emit: versions

    script:
    """
pip install pandas

python3 << 'EOF'
import pandas as pd
import os
def normalize_name(name):
    return os.path.splitext(os.path.basename(name))[0]
df = pd.read_csv("${mash_matrix}", sep='\t', index_col=0)
df.index = [normalize_name(x) for x in df.index]
df.columns = [normalize_name(x) for x in df.columns]
df.to_csv('mash_distances.tsv', sep='\t')
EOF

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed 's/Python //')
    pandas: \$(python -c "import pandas; print(pandas.__version__)")
END_VERSIONS
    """
}
