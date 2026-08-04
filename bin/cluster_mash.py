#!/usr/bin/env python3
"""
Cluster genomes based on Mash distances using single-linkage clustering.

Rewritten for large collections (2000+ Burkholderia pseudomallei assemblies)
while keeping the CLI surface and the default clustering semantics intact.

What changed
------------
1. Adjacency construction is a vectorized numpy mask instead of an
   O(n^2) Python double loop with `filtered_df.iloc[i, j]` per cell.
   Cluster labels are unchanged: run `cluster_mash.py <matrix> --self-test`
   to verify that on your own data.  The self-test re-implements the original
   upstream double loop inline and asserts that the connected-component labels,
   the component partition, and the order-split output all match.

2. Input can now be lower-triangular relaxed Phylip (`mash triangle`), the
   legacy labelled square TSV, or raw `mash dist` tabular output.  The format
   is sniffed, so the same script serves the old and the new front end.

3. Oversized connected components are split by phylogenetic similarity
   (average-linkage on the component's Mash submatrix) rather than by position
   in a Python list.  This is a *scientific* behaviour change and is gated on
   --split-method; pass `order` to reproduce the historical output exactly.

Memory note: the matrix is held dense as float64.  For n = 2000 that is
2000 * 2000 * 8 B = 32 MB (16 MB as float32), so there is no reason to
complicate this with a sparse or out-of-core representation at the scales this
workflow targets.  A dense n = 10000 matrix would still only be 800 MB.
"""

import argparse
import os
import time
import sys

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

# mash_matrix_io.py is shipped alongside this script in bin/
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mash_matrix_io import read_any, write_square, normalize_name  # noqa: E402


# --------------------------------------------------------------------------- #
# Connected components (single linkage at the threshold)
# --------------------------------------------------------------------------- #

def build_adjacency(matrix, threshold):
    """Boolean adjacency: edge iff distance <= threshold, off-diagonal, non-NaN.

    Vectorized equivalent of the legacy double loop.  NaN cells (pairs Mash did
    not report) are treated as "no edge", which is what
    `not np.isnan(filtered_df.iloc[i, j])` did.
    """
    with np.errstate(invalid="ignore"):
        mask = (matrix <= threshold) & np.isfinite(matrix)
    np.fill_diagonal(mask, False)
    rows, cols = np.nonzero(mask)
    adj = csr_matrix(
        (np.ones(rows.size, dtype=np.int8), (rows, cols)),
        shape=matrix.shape,
    )
    return adj


def components_in_label_order(adj, labels):
    """Return {component_label: [sample names]} preserving legacy ordering.

    The legacy code iterated `enumerate(labels)` and appended into a dict, so
    (a) members of a component come out in ascending matrix-index order and
    (b) components are keyed/ordered by first appearance.  Downstream cluster
    numbering depends on both, so both are reproduced here.
    """
    _, comp = connected_components(adj, directed=False)
    clusters = {}
    for idx, lab in enumerate(comp):
        clusters.setdefault(lab, []).append(labels[idx])
    return clusters, comp


# --------------------------------------------------------------------------- #
# Splitting oversized components
# --------------------------------------------------------------------------- #

def split_by_order(members, max_size):
    """Legacy behaviour: contiguous chunks of the member list.

    Retained only for reproducibility.  The chunk boundaries are arbitrary with
    respect to phylogeny: every chunk tends to contain a sample of *all* the
    lineages inside the component, which inflates within-cluster divergence and
    therefore both the Gubbins/IQ-TREE cost and the per-cluster recombination
    inference.
    """
    return [members[i:i + max_size] for i in range(0, len(members), max_size)]


def split_by_similarity(members, sub, max_size, consolidate=True,
                        consolidate_tolerance=2.0):
    """Average-linkage partition of one oversized component.

    `sub` is the component's square Mash submatrix, aligned with `members`.
    We cut the average-linkage dendrogram at the smallest number of flat
    clusters (k) for which no cluster exceeds max_size; if a k-cut still leaves
    an oversize cluster (possible with pathological ties), that cluster is
    recursively re-partitioned, and any residual oversize group falls back to
    order-based chunking so the max_size contract is never violated.

    With consolidate=True (default) a greedy merge pass then recombines the
    closest groups whose union still fits max_size, so the number of clusters
    stays near the ceil(n / max_size) floor instead of overshooting it - see
    consolidate_groups().
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    n = len(members)
    if n <= max_size:
        return [list(members)]

    # NaN / non-finite cells would poison the linkage. Mash reports every pair
    # inside a connected component, but be defensive: substitute the largest
    # observed finite distance, i.e. "maximally dissimilar".
    d = np.array(sub, dtype=np.float64, copy=True)
    d = 0.5 * (d + d.T)                       # enforce exact symmetry
    np.fill_diagonal(d, 0.0)
    bad = ~np.isfinite(d)
    if bad.any():
        finite_max = np.nanmax(d[np.isfinite(d)]) if np.isfinite(d).any() else 1.0
        d[bad] = finite_max
        np.fill_diagonal(d, 0.0)
    d[d < 0] = 0.0

    Z = linkage(squareform(d, checks=False), method="average")

    # Smallest k with all parts <= max_size. Monotone in k, so bisect.
    lo = int(np.ceil(n / float(max_size)))
    hi = n
    best = None
    while lo <= hi:
        mid = (lo + hi) // 2
        lab = fcluster(Z, t=mid, criterion="maxclust")
        sizes = np.bincount(lab)[1:]
        if sizes.max() <= max_size:
            best = lab
            hi = mid - 1
        else:
            lo = mid + 1

    if best is None:
        # fcluster could not honour max_size (heavy ties): recurse then fall back.
        lab = fcluster(Z, t=2, criterion="maxclust")
        out = []
        for c in np.unique(lab):
            sel = np.nonzero(lab == c)[0]
            part = [members[i] for i in sel]
            if len(part) == n:                        # no progress -> chunk it
                return split_by_order(list(members), max_size)
            if len(part) > max_size:
                out.extend(split_by_similarity(
                    part, d[np.ix_(sel, sel)], max_size,
                    consolidate=consolidate,
                    consolidate_tolerance=consolidate_tolerance))
            else:
                out.append(part)
        return out

    groups_idx = [np.nonzero(best == c)[0] for c in np.unique(best)]
    if consolidate:
        groups_idx = consolidate_groups(groups_idx, d, max_size,
                                        tolerance=consolidate_tolerance)

    groups = [[members[i] for i in sel] for sel in groups_idx]
    # Deterministic output order: by first member index within the component.
    order = {m: i for i, m in enumerate(members)}
    groups.sort(key=lambda g: min(order[m] for m in g))
    return groups


def _mean_within(idx, d):
    """Mean pairwise distance inside one index group (0.0 for singletons)."""
    if idx.size < 2:
        return 0.0
    sub = d[np.ix_(idx, idx)]
    iu = np.triu_indices(idx.size, k=1)
    return float(sub[iu].mean())


def consolidate_groups(groups_idx, d, max_size, tolerance=2.0):
    """Merge sibling sub-groups that are as tight together as they are inside.

    Cutting an average-linkage dendrogram at the smallest maxclust that honours
    max_size overshoots when the tree is unbalanced: on synthetic data with 8
    lineages of 250 and max_cluster_size=50 the cut produced 98 groups where 40
    would suffice.  Each extra group is another Snippy + Gubbins + IQ-TREE
    invocation, so on a modest workstation the group COUNT matters too.

    Merging purely by "closest pair that fits", however, buys a lower count with
    coherence: measured on the n=500 synthetic set (8 lineages of 62,
    max_cluster_size=50, so every lineage must span 2 groups) unconditional
    greedy merging cut 35 groups to 12 but dropped lineage purity from 1.000 to
    0.888 and raised mean within-cluster Mash distance from 0.0027 to 0.0055.
    That is the same failure mode as order-based chunking, just milder.

    So the merge is bounded instead of unconditional: two groups are merged only
    when their mean INTER-group distance is no more than `tolerance` times the
    larger of their two mean WITHIN-group distances, i.e. only when they look
    like siblings split by the max_size ceiling rather than distinct lineages.
    With the default tolerance of 2.0 and B. pseudomallei-like structure
    (within-lineage ~0.003 vs between-lineage ~0.022, a 7x gap) cross-lineage
    merges are refused, so the count drops without a purity cost.  Raise
    --consolidate-tolerance to trade coherence for fewer clusters; set
    --no-consolidate to skip the pass entirely.

    Performance note. The obvious implementation recomputes every inter-group
    mean, `d[np.ix_(groups[a], groups[b])].mean()`, on every pass of the merge
    loop, which is O(k^3) submatrix means for k initial groups. Profiled on the
    real 2795-genome set that was 290 million np.ix_/mean calls and 52 MINUTES
    for one clustering run, with `consolidate_groups` accounting for essentially
    all of it.

    Block means compose exactly, so nothing needs recomputing. Carrying SUMS
    instead of means:

        sum(A|B, C)  = sum(A, C) + sum(B, C)
        within(A|B)  = within_sum(A) + within_sum(B) + sum(A, B)

    a merge updates one row/column in O(k) and the pair scan becomes vectorized
    array lookups. The arithmetic is identical, not an approximation, so cluster
    membership is unchanged -- verified byte-for-byte against the previous
    implementation's output on the 2795-genome set at threshold 0.003.
    """
    groups = [np.asarray(g, dtype=np.int64) for g in groups_idx]
    k = len(groups)
    if k < 2:
        return groups

    sizes = np.array([g.size for g in groups], dtype=np.float64)

    # Sum of the within-group upper triangle, not the mean, so merges can add.
    within_sum = np.zeros(k, dtype=np.float64)
    for i, g in enumerate(groups):
        if g.size >= 2:
            sub = d[np.ix_(g, g)]
            within_sum[i] = float(sub[np.triu_indices(g.size, k=1)].sum())

    # Sum over the full a-by-b cross block. Computed once; O(k^2) submatrix
    # sums, i.e. what the old inner loop cost on a SINGLE pass.
    inter_sum = np.zeros((k, k), dtype=np.float64)
    for a in range(k):
        for b in range(a + 1, k):
            s = float(d[np.ix_(groups[a], groups[b])].sum())
            inter_sum[a, b] = s
            inter_sum[b, a] = s

    # `slots` preserves the ORDER the old list-based loop iterated in, so ties
    # resolve to the same pair: merged into position a, position b removed.
    slots = list(range(k))

    while len(slots) > 1:
        idx = np.asarray(slots, dtype=np.int64)
        s = sizes[idx]
        pairs = s * (s - 1.0) / 2.0
        within_mean = np.divide(within_sum[idx], pairs,
                                out=np.zeros_like(s), where=pairs > 0)
        inter_mean = inter_sum[np.ix_(idx, idx)] / np.outer(s, s)

        iu = np.triu_indices(idx.size, k=1)
        cand = inter_mean[iu]
        sa, sb = s[iu[0]], s[iu[1]]
        # Coherence bound. A pair of singletons (within == 0) has no internal
        # scale to compare against, so allow it: two lone genomes cannot be
        # "less coherent" than themselves.
        scale = np.maximum(within_mean[iu[0]], within_mean[iu[1]])
        ok = (sa + sb <= max_size) & ~((scale > 0.0) & (cand > tolerance * scale))
        if not ok.any():
            break

        # argmin takes the FIRST minimum, and triu_indices is row-major, which
        # together reproduce the old "a ascending, then b ascending, strict <".
        p = int(np.argmin(np.where(ok, cand, np.inf)))
        ai, bi = int(iu[0][p]), int(iu[1][p])
        A, B = int(idx[ai]), int(idx[bi])

        groups[A] = np.sort(np.concatenate([groups[A], groups[B]]))
        within_sum[A] += within_sum[B] + inter_sum[A, B]
        sizes[A] += sizes[B]
        merged_row = inter_sum[A, :] + inter_sum[B, :]
        inter_sum[A, :] = merged_row
        inter_sum[:, A] = merged_row
        inter_sum[A, A] = 0.0
        slots.pop(bi)

    return [groups[i] for i in slots]


def _consolidate_groups_reference(groups_idx, d, max_size, tolerance=2.0):
    """Original O(k^3) implementation, kept only as the correctness oracle.

    NOT called by the pipeline and deliberately unreferenced. It exists so the
    fast path above can be re-verified against it after any change:

        from cluster_mash import consolidate_groups, _consolidate_groups_reference
        assert [list(g) for g in consolidate_groups(gs, d, 50)] == \\
               [list(g) for g in _consolidate_groups_reference(gs, d, 50)]

    Equivalence was confirmed end to end on the real 2795-genome set: both
    implementations produce clusters_0.003.tsv byte-for-byte identically
    (282 clusters, 2795 genomes), in 10 s versus 3126 s.
    """
    groups = [np.asarray(g, dtype=np.int64) for g in groups_idx]
    within = [_mean_within(g, d) for g in groups]

    while len(groups) > 1:
        best_pair, best_dist = None, np.inf
        for a in range(len(groups)):
            for b in range(a + 1, len(groups)):
                if groups[a].size + groups[b].size > max_size:
                    continue
                inter = float(d[np.ix_(groups[a], groups[b])].mean())
                scale = max(within[a], within[b])
                if scale > 0.0 and inter > tolerance * scale:
                    continue
                if inter < best_dist:
                    best_dist, best_pair = inter, (a, b)
        if best_pair is None:
            return groups
        a, b = best_pair
        groups[a] = np.sort(np.concatenate([groups[a], groups[b]]))
        within[a] = _mean_within(groups[a], d)
        groups.pop(b)
        within.pop(b)
    return groups


# --------------------------------------------------------------------------- #
# Reporting helpers
# --------------------------------------------------------------------------- #

def mean_within_distance(matrix, idx_of, members):
    """Mean pairwise Mash distance inside one cluster (NaN-safe)."""
    if len(members) < 2:
        return 0.0
    sel = np.array([idx_of[m] for m in members if m in idx_of], dtype=np.int64)
    if sel.size < 2:
        return 0.0
    sub = matrix[np.ix_(sel, sel)]
    iu = np.triu_indices(sel.size, k=1)
    vals = sub[iu]
    vals = vals[np.isfinite(vals)]
    return float(vals.mean()) if vals.size else 0.0


def max_within_distance(matrix, idx_of, members):
    """Largest pairwise Mash distance inside one cluster (NaN-safe)."""
    if len(members) < 2:
        return 0.0
    sel = np.array([idx_of[m] for m in members if m in idx_of], dtype=np.int64)
    if sel.size < 2:
        return 0.0
    sub = matrix[np.ix_(sel, sel)]
    vals = sub[np.triu_indices(sel.size, k=1)]
    vals = vals[np.isfinite(vals)]
    return float(vals.max()) if vals.size else 0.0


def cluster_at(matrix, labels, idx_of, threshold, max_cluster_size,
               split_method, consolidate, consolidate_tolerance, verbose=True):
    """Connected components at `threshold`, oversized ones split. Returns dict."""
    adj = build_adjacency(matrix, threshold)
    components, _ = components_in_label_order(adj, labels)
    if verbose:
        print("Found %d connected components at threshold %s"
              % (len(components), threshold))

    final_clusters = {}
    counter = 0
    for _, members in components.items():
        if len(members) <= max_cluster_size:
            final_clusters[counter] = members
            counter += 1
            continue
        before = mean_within_distance(matrix, idx_of, members)
        if split_method == "order":
            parts = split_by_order(members, max_cluster_size)
        else:
            sel = np.array([idx_of[m] for m in members], dtype=np.int64)
            parts = split_by_similarity(members, matrix[np.ix_(sel, sel)],
                                        max_cluster_size,
                                        consolidate=consolidate,
                                        consolidate_tolerance=consolidate_tolerance)
        if verbose:
            after = [mean_within_distance(matrix, idx_of, p) for p in parts]
            print("  component of %d samples split into %d parts by '%s' "
                  "(mean within-component distance %.6f -> mean within-part %.6f)"
                  % (len(members), len(parts), split_method, before,
                     float(np.mean(after)) if after else 0.0))
        for part in parts:
            final_clusters[counter] = part
            counter += 1
    return final_clusters


def threshold_stats(final_clusters, matrix, idx_of):
    """Coverage and coherence of one clustering.

    analysable = genomes in clusters of >= 3 taxa; anything smaller cannot
    produce a tree and is dropped by the workflow.
    """
    analysable = dropped = 0
    w_max, w_means = 0.0, []
    n_clusters = 0
    for members in final_clusters.values():
        if len(members) >= 3:
            n_clusters += 1
            analysable += len(members)
            w_max = max(w_max, max_within_distance(matrix, idx_of, members))
            w_means.append(mean_within_distance(matrix, idx_of, members))
        else:
            dropped += len(members)
    return {
        "clusters": n_clusters,
        "analysable": analysable,
        "dropped": dropped,
        "within_max": w_max,
        "within_mean": float(np.mean(w_means)) if w_means else 0.0,
    }


def default_threshold_grid(matrix, n=15):
    """Candidate thresholds as quantiles of the observed pairwise distances.

    Deriving the grid from the data rather than hardcoding values is what makes
    this transferable: AUDIT_REPORT.md's 0.003 was calibrated on 112 genomes and
    silently discards 203 of 2795 on a wider collection.
    """
    iu = np.triu_indices(matrix.shape[0], k=1)
    vals = matrix[iu]
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return [0.03]
    qs = np.linspace(5, 95, n)
    grid = sorted(set(round(float(v), 6) for v in np.percentile(vals, qs)))
    return [g for g in grid if g > 0]


def select_threshold(matrix, labels, idx_of, grid, max_cluster_size,
                     split_method, consolidate, consolidate_tolerance,
                     coherence=0.90, report_path=None, refine=12):
    """Pick the threshold that analyses the most genomes without incoherent clusters.

    Two things are in tension. A LOW threshold leaves genomes in components of
    <3 taxa, which cannot produce a tree and are dropped outright. A HIGH one
    fuses distinct lineages into a single cluster, and Gubbins assumes samples of
    limited diversity sharing a recent common ancestor -- its recombination calls
    inside such a bin are not interpretable (AUDIT_REPORT.md G.3).

    The coherence bound is expressed RELATIVE to the collection's own diversity
    rather than as an absolute Mash distance, so it transfers to any dataset: a
    cluster may not span more than `coherence` of the largest pairwise distance
    in the whole collection. Measured on the 2795-genome B. pseudomallei set,
    whose total span is 0.01088:

        threshold  clusters>=3  dropped  within-MAX  span of collection
          0.005         67         58     0.00946        87%   accepted
          0.006         61          5     0.00946        87%   accepted <- chosen
          0.007         60          0     0.01061        98%   REJECTED
          0.008         60          0     0.01061        98%   REJECTED

    0.007 gains the last 5 genomes by creating one cluster that spans 98% of the
    entire collection. Among accepted thresholds the most genomes analysed wins,
    which is 0.006.
    """
    iu = np.triu_indices(matrix.shape[0], k=1)
    allv = matrix[iu]
    allv = allv[np.isfinite(allv)]
    span = float(allv.max()) if allv.size else 0.0
    limit = coherence * span

    print("Auto-threshold: collection span (max pairwise Mash) = %.6f, "
          "per-cluster limit = %.6f (%.0f%%)" % (span, limit, 100 * coherence))

    rows = []

    def evaluate(t, tag=""):
        fc = cluster_at(matrix, labels, idx_of, t, max_cluster_size,
                        split_method, consolidate, consolidate_tolerance,
                        verbose=False)
        st = threshold_stats(fc, matrix, idx_of)
        st["threshold"] = t
        st["ok"] = st["within_max"] <= limit and st["clusters"] > 0
        rows.append(st)
        print("  threshold %-9.6f clusters>=3=%-5d analysable=%-6d dropped=%-6d "
              "within_max=%.6f  %s%s"
              % (t, st["clusters"], st["analysable"], st["dropped"],
                 st["within_max"], "ok" if st["ok"] else "REJECTED (incoherent)",
                 tag))
        return st

    for t in grid:
        evaluate(t)

    # Refine above the best accepted threshold, where the coverage/coherence
    # trade-off actually lives. A UNIFORM scan, not a bisection: acceptance is
    # NOT monotonic in the threshold. Measured on the 2795-genome set, 0.005899
    # and 0.006180 both exceed the limit (within_max 0.0103) while 0.006000
    # between them does not (0.00946, 5 genomes dropped) -- one borderline genome
    # joins a different cluster there and swings that cluster's span across the
    # bound. A bisection assumes a single boundary and walks straight past such a
    # window; it converged on 0.005162 (38 dropped) and never sampled 0.006.
    accepted = [r["threshold"] for r in rows if r["ok"]]
    if accepted and refine > 0:
        lo = max(accepted)
        hi = max(r["threshold"] for r in rows)
        if hi > lo:
            step = (hi - lo) / float(refine + 1)
            print("  refining %.6f..%.6f in %d uniform steps "
                  "(acceptance is not monotonic, so this scans rather than bisects)"
                  % (lo, hi, refine))
            for i in range(1, refine + 1):
                evaluate(lo + i * step, tag="  [refine]")

    # Chosen over every threshold evaluated, coarse or refined. Acceptance is not
    # strictly monotonic (within_max can dip as a component re-splits), so this
    # scans all results rather than trusting the bisection endpoint.
    best = None
    for st in rows:
        if not st["ok"]:
            continue
        key = (-st["analysable"], st["within_max"], st["threshold"])
        if best is None or key < best[0]:
            best = (key, st)

    if report_path:
        with open(report_path, "w") as fh:
            fh.write("threshold\tclusters_ge3\tanalysable\tdropped\t"
                     "within_mean\twithin_max\tcollection_span\tlimit\taccepted\n")
            for st in rows:
                fh.write("%.6f\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%s\n"
                         % (st["threshold"], st["clusters"], st["analysable"],
                            st["dropped"], st["within_mean"], st["within_max"],
                            span, limit, "yes" if st["ok"] else "no"))
        print("Auto-threshold sweep written to %s" % report_path)

    if best is None:
        raise SystemExit(
            "ERROR: no candidate threshold produced coherent clusters "
            "(every option exceeded %.6f, %.0f%% of the collection span %.6f). "
            "Widen --auto-grid, raise --auto-coherence, or set --threshold "
            "explicitly." % (limit, 100 * coherence, span))

    chosen = best[1]
    print("Auto-threshold: selected %.6f (%d clusters >=3 taxa, %d genomes "
          "analysable, %d dropped, within_max %.6f = %.0f%% of collection span)"
          % (chosen["threshold"], chosen["clusters"], chosen["analysable"],
             chosen["dropped"], chosen["within_max"],
             100 * chosen["within_max"] / span if span else 0.0))
    return chosen["threshold"]


def write_clusters(final_clusters, output_file):
    with open(output_file, "w") as fh:
        fh.write("cluster_id\tsample_id\n")
        for cluster_id, members in final_clusters.items():
            for sample in members:
                fh.write("cluster_%s\t%s\n" % (cluster_id, normalize_name(sample)))


def write_submatrices(final_clusters, labels, matrix, outdir):
    """Emit one small square TSV per cluster.

    SELECT_CLUSTER_REPRESENTATIVE historically re-read the full n x n matrix in
    every per-cluster task (2000 x 2000 TSV parsed once per cluster).  Writing
    the submatrix once here makes each medoid task read only k x k.
    """
    os.makedirs(outdir, exist_ok=True)
    idx_of = {s: i for i, s in enumerate(labels)}
    written = []
    for cluster_id, members in final_clusters.items():
        sel = np.array([idx_of[m] for m in members if m in idx_of], dtype=np.int64)
        if sel.size == 0:
            continue
        path = os.path.join(outdir, "cluster_%s.matrix.tsv" % cluster_id)
        write_square([labels[i] for i in sel], matrix[np.ix_(sel, sel)], path)
        written.append(path)
    return written


# --------------------------------------------------------------------------- #

def _legacy_build_adjacency(df, threshold):
    """The ORIGINAL upstream build_distance_matrix(), copied verbatim.

    Kept only so --self-test can prove the vectorized path agrees with it on
    real input. Do not use in the pipeline: this is the O(n^2) Python double
    loop the rewrite exists to remove.
    """
    samples = list(df.index)
    filtered_df = df.copy()
    filtered_df[filtered_df > threshold] = np.nan
    rows, cols = [], []
    for i, _sample_i in enumerate(samples):
        for j, _sample_j in enumerate(samples):
            if i != j and not np.isnan(filtered_df.iloc[i, j]):
                rows.append(i)
                cols.append(j)
    data = [1] * len(rows)
    adj = csr_matrix((data, (rows, cols)), shape=(len(samples), len(samples)))
    return adj, samples


def self_test(matrix_file, threshold, max_cluster_size):
    """Assert the vectorized path reproduces the original implementation.

    Runs both the legacy double loop and the numpy mask on the SAME input and
    compares (a) connected-component labels elementwise, (b) the component
    partition, and (c) the final order-split cluster assignment. Exits non-zero
    on any mismatch, so it can be used as a CI check or a one-off sanity check
    on real data.
    """
    import pandas as pd

    labels, matrix = read_any(matrix_file)
    print("Self-test on %s: %d x %d matrix, threshold %s, max_cluster_size %s"
          % (matrix_file, matrix.shape[0], matrix.shape[1],
             threshold, max_cluster_size))

    # Legacy path wants a labelled DataFrame, as it got from the square TSV.
    df = pd.DataFrame(matrix, index=labels, columns=labels)

    t0 = time.time()
    adj_legacy, samples_legacy = _legacy_build_adjacency(df, threshold)
    t_legacy = time.time() - t0

    t0 = time.time()
    adj_new = build_adjacency(matrix, threshold)
    t_new = time.time() - t0

    ok = True

    if samples_legacy != labels:
        print("FAIL: label order differs")
        ok = False

    _, comp_legacy = connected_components(adj_legacy, directed=False)
    _, comp_new = connected_components(adj_new, directed=False)
    if np.array_equal(comp_legacy, comp_new):
        n_comp = int(comp_new.max()) + 1 if comp_new.size else 0
        print("PASS: connected-component labels identical (%d components)"
              % n_comp)
    else:
        n_diff = int((comp_legacy != comp_new).sum())
        print("FAIL: component labels differ in %d of %d positions"
              % (n_diff, comp_new.size))
        ok = False

    clusters_legacy, _ = components_in_label_order(adj_legacy, labels)
    clusters_new, _ = components_in_label_order(adj_new, labels)
    if ([list(v) for v in clusters_legacy.values()]
            == [list(v) for v in clusters_new.values()]):
        print("PASS: component partition and member ordering identical")
    else:
        print("FAIL: component partition or member ordering differs")
        ok = False

    # Order-split must reproduce the legacy output exactly.
    def split_all(clusters):
        out, counter = {}, 0
        for _, members in clusters.items():
            parts = (split_by_order(members, max_cluster_size)
                     if len(members) > max_cluster_size else [members])
            for part in parts:
                out[counter] = part
                counter += 1
        return out

    if ({k: list(v) for k, v in split_all(clusters_legacy).items()}
            == {k: list(v) for k, v in split_all(clusters_new).items()}):
        print("PASS: order-split cluster assignment identical "
              "(--split-method order reproduces legacy output)")
    else:
        print("FAIL: order-split cluster assignment differs")
        ok = False

    print("Timing: legacy double loop %.4f s vs numpy mask %.4f s (%.0fx)"
          % (t_legacy, t_new, (t_legacy / t_new) if t_new > 0 else float("nan")))
    print("SELF-TEST %s" % ("PASSED" if ok else "FAILED"))
    return 0 if ok else 1


def main():
    parser = argparse.ArgumentParser(
        description="Cluster genomes based on Mash distances"
    )
    parser.add_argument("mash_file",
                        help="Mash distances: lower-triangular Phylip "
                             "(mash triangle), labelled square TSV, or raw "
                             "5-column `mash dist` output. Auto-detected.")
    parser.add_argument("output_file", nargs="?", default=None,
                        help="Output cluster assignments file (not required "
                             "with --self-test)")
    parser.add_argument("--threshold", type=float, default=0.03,
                        help="Distance threshold for clustering (default: 0.03)")
    parser.add_argument("--max-cluster-size", type=int, default=100,
                        help="Maximum cluster size (default: 100)")
    parser.add_argument("--split-method", choices=["similarity", "order"],
                        default="similarity",
                        help="How to split components larger than "
                             "--max-cluster-size. 'similarity' (default) uses "
                             "average linkage on the Mash submatrix; 'order' "
                             "reproduces the legacy list-order chunking.")
    parser.add_argument("--self-test", action="store_true",
                        help="Verify that the vectorized adjacency reproduces "
                             "the original upstream double-loop implementation "
                             "on this input, then exit. Compares component "
                             "labels, the component partition and the "
                             "order-split assignment.")
    parser.add_argument("--no-consolidate", action="store_true",
                        help="With --split-method similarity, do NOT greedily "
                             "recombine adjacent sub-clusters that still fit "
                             "--max-cluster-size. Yields more, smaller clusters.")
    parser.add_argument("--consolidate-tolerance", type=float, default=2.0,
                        help="Merge two sub-clusters only if their mean "
                             "inter-cluster distance is <= this multiple of "
                             "their mean within-cluster distance (default: 2.0).")
    parser.add_argument("--emit-submatrices", metavar="DIR", default=None,
                        help="Also write one small square distance TSV per "
                             "cluster into DIR, for medoid selection.")
    parser.add_argument("--matrix-out", metavar="TSV", default=None,
                        help="Also write the full square matrix as TSV "
                             "(compatibility with consumers that still expect "
                             "mash_matrix.tsv).")

    parser.add_argument("--auto-threshold", action="store_true",
                        help="Sweep candidate thresholds and pick the one that "
                             "analyses the most genomes without producing a "
                             "cluster that spans more than --auto-coherence of "
                             "the collection's own diversity. Overrides "
                             "--threshold.")
    parser.add_argument("--auto-grid", default=None,
                        help="Comma-separated candidate thresholds for "
                             "--auto-threshold. Default: 15 quantiles of the "
                             "observed pairwise distances, so the grid adapts "
                             "to the collection.")
    parser.add_argument("--auto-coherence", type=float, default=0.90,
                        help="A cluster may not span more than this fraction "
                             "of the largest pairwise distance in the whole "
                             "collection (default: 0.90). Relative, not an "
                             "absolute Mash distance, so it transfers between "
                             "datasets.")
    parser.add_argument("--auto-refine", type=int, default=12,
                        help="Bisection steps between the best accepted and "
                             "first rejected threshold (default: 6). The "
                             "optimum sits at the coherence boundary, which a "
                             "coarse grid straddles. 0 disables refinement.")
    parser.add_argument("--auto-report", default=None,
                        help="Write the full threshold sweep to this TSV.")
    parser.add_argument("--threshold-out", default=None,
                        help="Write the chosen threshold to this file.")

    args = parser.parse_args()

    if args.self_test:
        sys.exit(self_test(args.mash_file, args.threshold,
                           args.max_cluster_size))

    if args.output_file is None:
        parser.error("output_file is required unless --self-test is given")

    print("Reading Mash distances from %s" % args.mash_file)
    labels, matrix = read_any(args.mash_file)
    print("Loaded %d x %d distance matrix (%.1f MB dense float64)"
          % (matrix.shape[0], matrix.shape[1], matrix.nbytes / 1e6))

    if args.matrix_out:
        write_square(labels, matrix, args.matrix_out)
        print("Wrote square matrix to %s" % args.matrix_out)

    idx_of = {s: i for i, s in enumerate(labels)}

    threshold = args.threshold
    if args.auto_threshold:
        grid = ([float(x) for x in args.auto_grid.split(",")]
                if args.auto_grid else default_threshold_grid(matrix))
        threshold = select_threshold(
            matrix, labels, idx_of, grid, args.max_cluster_size,
            args.split_method, not args.no_consolidate,
            args.consolidate_tolerance,
            coherence=args.auto_coherence,
            report_path=args.auto_report,
            refine=args.auto_refine)
        if args.threshold_out:
            with open(args.threshold_out, "w") as fh:
                fh.write("%.6f\n" % threshold)
            print("Chosen threshold written to %s" % args.threshold_out)

    print("Clustering with threshold %s" % threshold)
    final_clusters = cluster_at(matrix, labels, idx_of, threshold,
                                args.max_cluster_size, args.split_method,
                                not args.no_consolidate,
                                args.consolidate_tolerance, verbose=True)

    print("Found %d clusters" % len(final_clusters))
    for cluster_id, members in final_clusters.items():
        print("  Cluster %d: %d samples (mean within-cluster Mash distance %.6f)"
              % (cluster_id, len(members),
                 mean_within_distance(matrix, idx_of, members)))

    write_clusters(final_clusters, args.output_file)
    print("Cluster assignments written to %s" % args.output_file)

    if args.emit_submatrices:
        paths = write_submatrices(final_clusters, labels, matrix,
                                  args.emit_submatrices)
        print("Wrote %d per-cluster submatrices to %s"
              % (len(paths), args.emit_submatrices))


if __name__ == "__main__":
    main()
