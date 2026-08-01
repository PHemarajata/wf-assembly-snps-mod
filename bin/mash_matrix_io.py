#!/usr/bin/env python3
"""
Shared, vectorized readers/writers for Mash distance data.

Three on-disk shapes are supported so that the clustering front end can be
switched between the legacy `mash dist` path and the new `mash triangle` path
without touching downstream code:

  1. "phylip"  - lower-triangular relaxed Phylip, i.e. the output of
                 `mash triangle`.  First line is the taxon count, then one line
                 per taxon: <label> <d_i0> ... <d_i,i-1>.
  2. "square"  - a labelled square TSV (row 0 = column names, column 0 = row
                 names).  This is what bin/mash_tab_to_matrix.py has always
                 produced and what cluster_mash.py / SELECT_CLUSTER_REPRESENTATIVE
                 have always consumed.
  3. "edges"   - raw 5-column `mash dist` tabular output
                 (ref, query, distance, p_value, shared_hashes).

All three are read into the same (labels, numpy float64 square matrix) pair, so
callers do not need to know which one they were handed.  Nothing here loops in
Python over O(n^2) cells: the legacy `for _, row in df.iterrows(): matrix.at[...]`
pattern is replaced with numpy fancy indexing / pivot.
"""

import os
import numpy as np


def normalize_name(name):
    """Strip directory and a single extension - matches the historical helper."""
    return os.path.splitext(os.path.basename(str(name)))[0]


def sniff_format(path):
    """Return 'phylip', 'edges' or 'square' by inspecting the first two lines."""
    with open(path) as fh:
        first = fh.readline()
        second = fh.readline()

    # `mash triangle` writes the taxon count alone on line 1 (often tab-indented).
    if first.strip().isdigit():
        return "phylip"

    # `mash dist` writes 5 tab-separated fields with no header.
    fields = first.rstrip("\n").split("\t")
    if len(fields) == 5:
        try:
            float(fields[2])
            return "edges"
        except ValueError:
            pass

    # Otherwise assume a labelled square matrix. Guard against a headerless file.
    if second and len(second.rstrip("\n").split("\t")) == len(fields):
        return "square"
    return "square"


def read_phylip_lower(path, normalize=True):
    """Read lower-triangular relaxed Phylip into (labels, square float64 matrix).

    `mash triangle` labels rows by the *input file path*, so name normalization
    still has to happen here exactly as it did in mash_tab_to_matrix.py.
    """
    with open(path) as fh:
        lines = fh.read().splitlines()

    # Drop blank lines but keep order; line 0 is the count.
    lines = [ln for ln in lines if ln.strip() != ""]
    n = int(lines[0].strip())
    body = lines[1:]
    if len(body) != n:
        raise ValueError(
            "Phylip header declares %d taxa but file has %d data lines (%s)"
            % (n, len(body), path)
        )

    labels = []
    matrix = np.zeros((n, n), dtype=np.float64)
    for i, ln in enumerate(body):
        parts = ln.split("\t") if "\t" in ln else ln.split()
        labels.append(normalize_name(parts[0]) if normalize else parts[0])
        if i:
            # Row i carries exactly i distances: d(i,0) .. d(i,i-1).
            vals = np.array(parts[1:1 + i], dtype=np.float64)
            if vals.size != i:
                raise ValueError(
                    "Phylip row %d has %d values, expected %d (%s)"
                    % (i, vals.size, i, path)
                )
            matrix[i, :i] = vals
            matrix[:i, i] = vals
    np.fill_diagonal(matrix, 0.0)
    return labels, matrix


def read_square(path, normalize=True):
    """Read a labelled square TSV into (labels, square float64 matrix)."""
    import pandas as pd

    df = pd.read_csv(path, sep="\t", index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce")
    labels = [normalize_name(x) for x in df.index] if normalize else list(df.index)
    return labels, df.to_numpy(dtype=np.float64)


def read_mash_edges(path, normalize=True):
    """Read 5-column `mash dist` output into (labels, square float64 matrix).

    Vectorized replacement for the row-at-a-time `.at[]` fill in
    bin/mash_tab_to_matrix.py.  Missing pairs stay NaN, exactly as before;
    the diagonal is forced to 0.
    """
    import pandas as pd

    df = pd.read_csv(
        path, sep="\t", header=None,
        names=["ref", "query", "distance", "p_value", "shared_hashes"],
        dtype={"ref": str, "query": str},
    )
    df = df[df["ref"].notna() & df["query"].notna()]
    df["distance"] = pd.to_numeric(df["distance"], errors="coerce")

    if normalize:
        # .map over the unique values only - normalize_name on 4M rows one at a
        # time is what made the original slow.
        ref_u = {v: normalize_name(v) for v in df["ref"].unique()}
        qry_u = {v: normalize_name(v) for v in df["query"].unique()}
        ref = df["ref"].map(ref_u).to_numpy()
        qry = df["query"].map(qry_u).to_numpy()
    else:
        ref = df["ref"].to_numpy()
        qry = df["query"].to_numpy()

    labels = sorted(set(ref.tolist()) | set(qry.tolist()))
    idx = {s: i for i, s in enumerate(labels)}
    n = len(labels)

    ri = np.fromiter((idx[s] for s in ref), dtype=np.int64, count=ref.size)
    qi = np.fromiter((idx[s] for s in qry), dtype=np.int64, count=qry.size)
    vals = df["distance"].to_numpy(dtype=np.float64)

    matrix = np.full((n, n), np.nan, dtype=np.float64)
    matrix[ri, qi] = vals
    matrix[qi, ri] = vals          # symmetric, same as the legacy double assignment
    np.fill_diagonal(matrix, 0.0)
    return labels, matrix


def read_any(path, normalize=True):
    """Dispatch on sniff_format and return (labels, square float64 matrix)."""
    fmt = sniff_format(path)
    if fmt == "phylip":
        return read_phylip_lower(path, normalize=normalize)
    if fmt == "edges":
        return read_mash_edges(path, normalize=normalize)
    return read_square(path, normalize=normalize)


def write_square(labels, matrix, path, float_format="%.6g"):
    """Write a labelled square TSV byte-compatible with mash_tab_to_matrix.py."""
    import pandas as pd

    df = pd.DataFrame(matrix, index=labels, columns=labels)
    df.to_csv(path, sep="\t", float_format=float_format)


def write_submatrix(labels, matrix, keep, path):
    """Write the square submatrix restricted to `keep` (list of labels)."""
    idx = {s: i for i, s in enumerate(labels)}
    sel = np.array([idx[s] for s in keep if s in idx], dtype=np.int64)
    sub = matrix[np.ix_(sel, sel)]
    write_square([labels[i] for i in sel], sub, path)
    return sel.size
