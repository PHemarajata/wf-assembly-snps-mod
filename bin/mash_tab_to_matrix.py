#!/usr/bin/env python3
"""
Convert Mash tabular output to a square distance matrix for downstream clustering/medoid selection.
Usage:
    python mash_tab_to_matrix.py mash_tabular.tsv mash_matrix.tsv
"""
import sys
import pandas as pd
import numpy as np
import os

def normalize_name(name):
    return os.path.splitext(os.path.basename(str(name)))[0]

def main(tabular_file, matrix_file):
    # Mash tabular columns: ref, query, distance, p_value, shared_hashes
    df = pd.read_csv(tabular_file, sep='\t', header=None,
                     names=['ref', 'query', 'distance', 'p_value', 'shared_hashes'])
    # Filter out rows where ref or query is not a string
    df = df[df['ref'].apply(lambda x: isinstance(x, str))]
    df = df[df['query'].apply(lambda x: isinstance(x, str))]
    # Normalize sample names
    df['ref'] = df['ref'].apply(normalize_name)
    df['query'] = df['query'].apply(normalize_name)
    # Get all unique sample names
    samples = sorted(set(df['ref'].tolist() + df['query'].tolist()))
    # Initialize square matrix
    matrix = pd.DataFrame(np.nan, index=samples, columns=samples)
    # Fill in distances
    for _, row in df.iterrows():
        matrix.at[row['ref'], row['query']] = row['distance']
        matrix.at[row['query'], row['ref']] = row['distance']  # symmetric
    # Fill diagonal with zeros
    np.fill_diagonal(matrix.values, 0)
    # Normalize matrix row/col names
    matrix.index = [normalize_name(x) for x in matrix.index]
    matrix.columns = [normalize_name(x) for x in matrix.columns]
    # Save as TSV
    matrix.to_csv(matrix_file, sep='\t')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python mash_tab_to_matrix.py mash_tabular.tsv mash_matrix.tsv")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
