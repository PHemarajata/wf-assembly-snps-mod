#!/usr/bin/env python3
"""
Convert Mash tabular output to a square distance matrix for downstream
clustering / medoid selection.

Behaviour-compatible drop-in for the original script.  The original filled the
matrix with

    for _, row in df.iterrows():
        matrix.at[row['ref'], row['query']] = row['distance']
        matrix.at[row['query'], row['ref']] = row['distance']

which is O(n^2) Python-level iterations with two label lookups each: at n = 2000
that is 4,000,000 rows.  The fill is now a single pair of numpy fancy-index
assignments (see bin/mash_matrix_io.read_mash_edges).

Semantics preserved exactly:
  - sample labels are basename-minus-one-extension, sorted;
  - unreported pairs remain NaN;
  - the matrix is symmetrized by assigning both (i,j) and (j,i);
  - the diagonal is forced to 0.

Usage:
    python3 mash_tab_to_matrix.py mash_tabular.tsv mash_matrix.tsv
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mash_matrix_io import read_mash_edges, write_square  # noqa: E402


def main(tabular_file, matrix_file):
    labels, matrix = read_mash_edges(tabular_file)
    write_square(labels, matrix, matrix_file)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python mash_tab_to_matrix.py mash_tabular.tsv mash_matrix.tsv")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
