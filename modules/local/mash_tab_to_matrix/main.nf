process MASH_TAB_TO_MATRIX {
    tag "mash_tab_to_matrix"
    label 'process_low'
    container "quay.io/biocontainers/python:3.9--1"

    input:
    path mash_tabular

    output:
    path "mash_matrix.tsv", emit: matrix
    path "versions.yml", emit: versions

    script:
    """
pip install pandas numpy

python3 ${projectDir}/bin/mash_tab_to_matrix.py \
    $mash_tabular \
    mash_matrix.tsv

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed 's/Python //')
    pandas: \$(python -c "import pandas; print(pandas.__version__)")
    numpy: \$(python -c "import numpy; print(numpy.__version__)")
END_VERSIONS
    """
}
