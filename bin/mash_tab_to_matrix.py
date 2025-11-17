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
    print(f"Reading tabular file: {tabular_file}")

    # Mash tabular columns: ref, query, distance, p_value, shared_hashes
    df = pd.read_csv(tabular_file, sep='\t', header=None,
                     names=['ref', 'query', 'distance', 'p_value', 'shared_hashes'])

    print(f"Loaded {len(df)} pairwise distances")

    # Filter out rows where ref or query is not a string
    df = df[df['ref'].apply(lambda x: isinstance(x, str))]
    df = df[df['query'].apply(lambda x: isinstance(x, str))]

    # Normalize sample names
    df['ref'] = df['ref'].apply(normalize_name)
    df['query'] = df['query'].apply(normalize_name)

    # Get all unique sample names
    samples = sorted(set(df['ref'].tolist() + df['query'].tolist()))
    print(f"Found {len(samples)} unique samples")

    # Create sample to index mapping for fast lookup
    sample_to_idx = {sample: idx for idx, sample in enumerate(samples)}

    # Initialize square matrix with zeros (faster than NaN)
    n_samples = len(samples)
    matrix_array = np.zeros((n_samples, n_samples), dtype=np.float64)

    print("Building distance matrix using vectorized operations...")

    # Vectorized approach - much faster than iterrows()
    ref_indices = df['ref'].map(sample_to_idx).values
    query_indices = df['query'].map(sample_to_idx).values
    distances = df['distance'].values

    # Fill matrix in both directions (symmetric)
    matrix_array[ref_indices, query_indices] = distances
    matrix_array[query_indices, ref_indices] = distances

    # Diagonal is already zeros from initialization

    print("Converting to DataFrame...")
    # Convert to pandas DataFrame for output
    matrix = pd.DataFrame(matrix_array, index=samples, columns=samples)

    print(f"Saving matrix to: {matrix_file}")
    # Save as TSV
    matrix.to_csv(matrix_file, sep='\t')

    print("Done!")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python mash_tab_to_matrix.py mash_tabular.tsv mash_matrix.tsv")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
