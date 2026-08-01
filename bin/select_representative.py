#!/usr/bin/env python3
"""
Pick the medoid of one cluster from a Mash distance matrix and stage it as the
per-cluster reference assembly.

This is the inline heredoc that used to live in
modules/local/select_cluster_representative/main.nf, lifted into bin/ so that it
can be unit-tested outside a container and so the module no longer has to
`pip install pandas numpy` at runtime.

Behaviour preserved from the original, including the fallbacks:
  - labels are basename-minus-one-extension;
  - a version-tolerant second matching pass (a matrix label 'GCA_000259775' must
    match a staged file 'GCA_000259775.1.fasta');
  - if the representative cannot be matched to a staged file, a real assembly
    from the cluster is copied instead of emitting a 1-bp placeholder, because a
    degenerate reference collapses the per-cluster Snippy alignment and makes
    Gubbins fail.

Performance change: the matrix argument is normally the CLUSTER'S OWN
submatrix (k x k), written once by cluster_mash.py --emit-submatrices, instead of
the full n x n matrix re-parsed in every per-cluster task.  A full matrix is
still accepted; only the cluster's rows/columns are read from it, via pandas
usecols so the unused n-k columns are never materialized.  The medoid search
itself is O(k^2), which is negligible at k <= max_cluster_size.
"""

import argparse
import os
import shutil
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mash_matrix_io import normalize_name, sniff_format, read_phylip_lower  # noqa: E402


def version_stripped(name):
    """Drop a trailing assembly version ('.1' or '_1', up to 3 digits)."""
    b = normalize_name(name)
    if "." in b:
        head, tail = b.rsplit(".", 1)
        if tail.isdigit() and len(tail) <= 3:
            return head
    if "_" in b:
        head, tail = b.rsplit("_", 1)
        if tail.isdigit() and len(tail) <= 3:
            return head
    return b


def load_cluster_block(matrix_file, wanted):
    """Return (labels, square float64 block) restricted to `wanted`.

    For a square TSV only the needed columns are parsed (usecols), so handing
    this the full n x n matrix costs O(n * k) instead of O(n^2).
    """
    import pandas as pd

    fmt = sniff_format(matrix_file)
    if fmt == "phylip":
        labels, matrix = read_phylip_lower(matrix_file)
        idx = {s: i for i, s in enumerate(labels)}
        sel = [idx[s] for s in wanted if s in idx]
        return [labels[i] for i in sel], matrix[np.ix_(sel, sel)]

    header = pd.read_csv(matrix_file, sep="\t", nrows=0)
    all_cols = list(header.columns)
    index_col = all_cols[0]
    norm = {normalize_name(c): c for c in all_cols[1:]}
    keep = [norm[s] for s in wanted if s in norm]
    if not keep:
        return [], np.zeros((0, 0))

    df = pd.read_csv(matrix_file, sep="\t", usecols=[index_col] + keep,
                     index_col=0)
    df.index = [normalize_name(x) for x in df.index]
    df.columns = [normalize_name(c) for c in df.columns]
    rows = [s for s in wanted if s in df.index]
    df = df.loc[rows, rows].apply(pd.to_numeric, errors="coerce")
    return list(df.index), df.to_numpy(dtype=np.float64)


def select_medoid(cluster_id, sample_ids, matrix_file):
    print("Selecting representative for cluster %s" % cluster_id)
    print("Samples in cluster: %s" % (sample_ids,))

    normalized_ids = [normalize_name(s) for s in sample_ids]
    if len(normalized_ids) == 1:
        return normalized_ids[0]

    try:
        labels, block = load_cluster_block(matrix_file, normalized_ids)
        print("Loaded %d x %d distance block for this cluster" % block.shape)
    except Exception as exc:                                # noqa: BLE001
        print("Error reading distance matrix: %s" % exc)
        return normalized_ids[0]

    if len(labels) < len(normalized_ids):
        print("Warning: only %d of %d samples found in distance matrix"
              % (len(labels), len(normalized_ids)))
    if not labels:
        print("No samples found in distance matrix, selecting first sample")
        return normalized_ids[0]
    if len(labels) == 1:
        return labels[0]

    # NaN-safe row sums; identical to pandas .sum(axis=1) which skips NaN.
    sums = np.nansum(np.where(np.isfinite(block), block, np.nan), axis=1)
    best = int(np.argmin(sums))
    print("Selected medoid: %s (sum of distances: %.6f)" % (labels[best], sums[best]))
    return labels[best]


def stage_representative(representative_id, cluster_id):
    exts = (".fa", ".fasta", ".fna")
    assembly_files = sorted(f for f in os.listdir(".") if f.endswith(exts))
    # Do not let a previously written representative_id.fa be picked as input.
    assembly_files = [f for f in assembly_files
                      if f != "%s.fa" % representative_id]

    representative_file = None
    for f in assembly_files:                       # 1) exact normalized match
        if normalize_name(f) == representative_id:
            representative_file = f
            break
    if representative_file is None:                # 2) version-tolerant match
        for f in assembly_files:
            if (version_stripped(f) == representative_id
                    or version_stripped(f) == version_stripped(representative_id)):
                representative_file = f
                break

    out = "%s.fa" % representative_id
    if representative_file:
        shutil.copy(representative_file, out)
        print("Copied %s to %s" % (representative_file, out))
    elif assembly_files:
        shutil.copy(assembly_files[0], out)
        print("WARNING: could not match representative %s; using %s instead"
              % (representative_id, assembly_files[0]))
    else:
        print("ERROR: no assembly files staged for cluster %s" % cluster_id)
        with open(out, "w") as fh:
            fh.write(">%s\nN\n" % representative_id)
    return out


def main():
    ap = argparse.ArgumentParser(description="Select a cluster medoid representative")
    ap.add_argument("--cluster-id", required=True)
    ap.add_argument("--sample-ids", required=True,
                    help="Nextflow list literal, e.g. \"[a, b, c]\"")
    ap.add_argument("--matrix", required=True,
                    help="Cluster submatrix (preferred) or full square matrix "
                         "or lower-triangular Phylip")
    ap.add_argument("--out-id", default="representative_id.txt")
    args = ap.parse_args()

    sample_ids = [s for s in args.sample_ids.strip("[]").replace(" ", "").split(",") if s]

    representative_id = select_medoid(args.cluster_id, sample_ids, args.matrix)
    print("Representative for cluster %s: %s" % (args.cluster_id, representative_id))

    with open(args.out_id, "w") as fh:
        fh.write(representative_id + "\n")

    stage_representative(representative_id, args.cluster_id)


if __name__ == "__main__":
    main()
