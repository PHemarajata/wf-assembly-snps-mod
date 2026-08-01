#!/usr/bin/env python3
"""
Convert `mash triangle` lower-triangular relaxed Phylip to a labelled square TSV.

Emitted for backwards compatibility with consumers that expect the historical
mash_matrix.tsv shape produced by bin/mash_tab_to_matrix.py.  cluster_mash.py
reads the Phylip directly and does not need this.

Usage:
    python3 mash_phylip_to_matrix.py mash_distances.phylip mash_matrix.tsv
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mash_matrix_io import read_phylip_lower, write_square  # noqa: E402


def main(phylip_file, matrix_file):
    labels, matrix = read_phylip_lower(phylip_file)
    write_square(labels, matrix, matrix_file)
    print("Wrote %d x %d square matrix to %s"
          % (matrix.shape[0], matrix.shape[1], matrix_file))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
