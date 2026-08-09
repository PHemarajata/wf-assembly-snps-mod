#!/usr/bin/env python3
"""
Root-to-tip regression against collection date, whole-tree and per-cluster.

Run this BEFORE any dating analysis. A tip-dated clock model will happily return
a posterior with tight credible intervals on data that contain no temporal
signal at all; the regression is the cheap check that says whether the data can
support the question.

Two levels are reported because they answer different questions on a GRAFTED
tree:

  whole tree   Meaningless if the tree was assembled by grafting per-cluster
               trees, each inferred from its own alignment: root-to-tip distance
               then mixes incomparable branch-length scales. Reported anyway,
               because the number is the argument for not dating the whole tree.

  per cluster  Valid. Within one grafted cluster all branch lengths come from
               one alignment, and the shared backbone path from the root to the
               cluster is a constant that cancels out of both the correlation
               and the slope. A cluster with a decent positive slope is a
               candidate for its own dated analysis.

Root-to-tip distance is taken from the tree AS ROOTED IN THE FILE, or from a
precomputed column with `--rtt-col`. Re-root the tree first if you want a
different root; the per-cluster result is insensitive to this, the whole-tree
result is not.

    clock_signal_check.py --tree global_grafted.treefile \
        --annotation Tip_annotation/tip_annotation.tsv \
        --out Tip_annotation/clock_signal_by_cluster.tsv
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------
# A minimal rooted-tree reader. Only root-to-tip distance is needed, so the
# tree is stored as parent pointers rather than as a node object graph.
# --------------------------------------------------------------------------
def tokenize(text: str):
    text = re.sub(r"\[[^\]]*\]", "", text)
    i, n = 0, len(text)
    while i < n:
        ch = text[i]
        if ch in "(),:;":
            yield ch
            i += 1
        elif ch in "'\"":
            quote, i = ch, i + 1
            start = i
            while i < n and text[i] != quote:
                i += 1
            yield ("name", text[start:i])
            i += 1
        elif ch.isspace():
            i += 1
        else:
            start = i
            while i < n and text[i] not in "(),:;" and not text[i].isspace():
                i += 1
            yield ("name", text[start:i])


def root_to_tip_distances(path: Path) -> dict[str, float]:
    """Newick -> {leaf label: distance from the file's root}.

    Recursive descent over the token stream, accumulating depth on the way down,
    so nothing but the leaf distances is ever materialised.
    """
    tokens = list(tokenize(path.read_text()))
    pos = 0

    def peek():
        return tokens[pos] if pos < len(tokens) else None

    def parse_node():
        """-> (name, branch_length, [children]); children empty for a leaf."""
        nonlocal pos
        children = []
        if peek() == "(":
            pos += 1
            while True:
                children.append(parse_node())
                if peek() == ",":
                    pos += 1
                    continue
                break
            if peek() == ")":
                pos += 1
        name = ""
        if isinstance(peek(), tuple):
            name = peek()[1]
            pos += 1
        length = 0.0
        if peek() == ":":
            pos += 1
            token = peek()
            if isinstance(token, tuple):
                try:
                    length = float(token[1])
                except ValueError:
                    length = 0.0
                pos += 1
        return name, length, children

    root = parse_node()

    distances: dict[str, float] = {}
    duplicates: list[str] = []
    stack = [(root, 0.0)]
    while stack:
        (name, length, children), depth = stack.pop()
        here = depth + length
        if children:
            stack.extend((child, here) for child in children)
        elif name:
            if name in distances:
                duplicates.append(name)
            distances[name] = here
    if duplicates:
        print(f"WARNING: {len(duplicates)} duplicate tip labels in the tree, "
              f"e.g. {duplicates[:5]}", file=sys.stderr)
    return distances


def regress(x: np.ndarray, y: np.ndarray):
    if len(x) < 3 or len(set(x)) < 3:
        return None
    r = float(np.corrcoef(x, y)[0, 1])
    slope, intercept = np.polyfit(x, y, 1)
    tmrca = -intercept / slope if slope != 0 else float("nan")
    return r, float(slope), float(tmrca)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotation", required=True,
                    help="tip_annotation.tsv from build_tip_annotation.py")
    ap.add_argument("--tree", help="Newick tree; root-to-tip is computed from it")
    ap.add_argument("--rtt-col", default="prior_root_to_tip",
                    help="use a precomputed root-to-tip column instead of --tree")
    ap.add_argument("--date-col", default="date_decimal")
    ap.add_argument("--cluster-col", default="cluster_id")
    ap.add_argument("--min-tips", type=int, default=10)
    ap.add_argument("--min-dates", type=int, default=5,
                    help="minimum distinct collection dates in a cluster")
    ap.add_argument("--sites-col", default="cluster_n_filtered_polymorphic_sites",
                    help="per-cluster count of variable sites the tree was inferred from; "
                         "used to convert the slope to substitutions per genome per year")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    table = pd.read_csv(args.annotation, sep="\t")
    if args.tree:
        rtt = root_to_tip_distances(Path(args.tree))
        table["_rtt"] = table["tip"].map(rtt)
        missing = int(table["_rtt"].isna().sum())
        if missing:
            print(f"WARNING: {missing} annotation rows have no matching tip in the tree",
                  file=sys.stderr)
        source = f"tree:{args.tree}"
    else:
        if args.rtt_col not in table.columns:
            sys.exit(f"ERROR: no --tree and no column '{args.rtt_col}' in the annotation")
        table["_rtt"] = pd.to_numeric(table[args.rtt_col], errors="coerce")
        source = f"column:{args.rtt_col}"

    data = table.dropna(subset=["_rtt", args.date_col])
    whole = regress(data[args.date_col].to_numpy(float), data["_rtt"].to_numpy(float))

    rows = []
    for cluster, group in data.groupby(args.cluster_col):
        if len(group) < args.min_tips or group[args.date_col].nunique() < args.min_dates:
            continue
        fit = regress(group[args.date_col].to_numpy(float), group["_rtt"].to_numpy(float))
        if fit is None:
            continue
        r, slope, tmrca = fit
        # The per-cluster tree is inferred from the Gubbins filtered_polymorphic_
        # sites alignment, so a branch length is substitutions per VARIABLE site.
        # Substitutions per genome per year is therefore slope x n_variable_sites,
        # not slope x genome_length. Reported so the rate can be compared against
        # the literature -- roughly 1-10 SNPs/genome/year for most bacteria.
        sites = np.nan
        if args.sites_col in group.columns:
            sites = pd.to_numeric(group[args.sites_col], errors="coerce").dropna()
            sites = float(sites.iloc[0]) if len(sites) else np.nan
        rows.append({
            "cluster_id": cluster,
            "n_dated_tips": len(group),
            "n_distinct_dates": int(group[args.date_col].nunique()),
            "date_span_years": round(float(group[args.date_col].max() - group[args.date_col].min()), 2),
            "r": round(r, 4),
            "r_squared": round(r * r, 4),
            "slope_subs_per_alignment_site_per_year": f"{slope:.4g}",
            "n_variable_sites": "" if np.isnan(sites) else int(sites),
            "subs_per_genome_per_year": "" if np.isnan(sites) else round(slope * sites, 2),
            "implied_tmrca_year": round(tmrca, 1) if np.isfinite(tmrca) else "",
            "temporal_signal": "usable" if r >= 0.3 else ("weak" if r > 0 else "none/negative"),
        })

    out = pd.DataFrame(rows).sort_values("r", ascending=False)
    out.to_csv(args.out, sep="\t", index=False)

    print(f"root-to-tip source : {source}")
    print(f"dated tips         : {len(data)} of {len(table)}")
    if whole:
        r, slope, tmrca = whole
        print(f"\nWHOLE TREE   r = {r:+.4f}   R2 = {r * r:.4f}   "
              f"slope = {slope:.3g} subs/site/yr")
        print("  On a grafted tree this number is not interpretable as a clock rate; "
              "an |r| near zero here\n  is the expected consequence of mixing "
              "per-cluster branch-length scales, not evidence of\n  a slow clock.")
    print(f"\nPER CLUSTER  {len(out)} clusters with >= {args.min_tips} dated tips "
          f"and >= {args.min_dates} distinct dates")
    if len(out):
        print(f"  usable (r >= 0.3)  : {(out['temporal_signal'] == 'usable').sum()}")
        print(f"  weak   (0 < r < .3): {(out['temporal_signal'] == 'weak').sum()}")
        print(f"  none/negative      : {(out['temporal_signal'] == 'none/negative').sum()}")
        print(f"  median r           : {out['r'].median():+.3f}")
        print("\n" + out.head(15).to_string(index=False))
    print(f"\nwrote -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
