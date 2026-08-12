#!/usr/bin/env python3
"""
Convert a PopPUNK clusters CSV (Taxon,Cluster) into the pipeline's
`cluster_id<TAB>sample_id` contract, applying a minimum-cluster-size gate.

The gate is the point of this script as much as the format change. PopPUNK
partitions everything -- on the 2,802-genome B. pseudomallei collection its
refined fit assigned 100% of genomes to 271 strains, but 158 of those are
singletons and only 34 strains have n >= 7. A singleton is not an analysis unit:
Gubbins needs taxa to detect recombination against, and a 1- or 2-genome
"cluster" produces a Tier4 result at best. Filtering here rather than letting
those clusters run and fail keeps the failure visible (they are written to the
excluded TSV with a reason) instead of buried in per-cluster diagnostics.

Cluster ids are prefixed (default "strain_") because a bare integer is a poor
directory name and collides with the `cluster_<id>` publishDir convention.
"""

import argparse
import csv
import re
import sys
from collections import defaultdict


def sanitise(name):
    """Reproduce PopPUNK's sample-name rewriting.

    PopPUNK rewrites characters it will not carry through its database into
    underscores, so `GCF_015714675_1_Virgin_Islands_St._John` comes back as
    `..._St__John`. Exactly one genome in the 2,802-genome collection is
    affected -- which is precisely why it matters: the sample keeps its original
    name in the samplesheet, so an unreconciled join drops it from its cluster
    silently, and a one-genome loss in a 2,382-genome run is invisible. Map
    PopPUNK's names back to the originals rather than trusting them to match.
    """
    return re.sub(r"[^A-Za-z0-9_-]", "_", name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters", required=True,
                    help="PopPUNK *_clusters.csv with Taxon,Cluster columns")
    ap.add_argument("--rfile",
                    help="rfile used for the run; if given, sample names are "
                         "checked against it and any mismatch is reported "
                         "rather than silently dropped")
    ap.add_argument("--min-cluster-size", type=int, default=7)
    ap.add_argument("--prefix", default="strain_")
    ap.add_argument("--out", required=True)
    ap.add_argument("--excluded")
    a = ap.parse_args()

    with open(a.clusters) as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        sys.exit("ERROR: %s has no rows" % a.clusters)

    cols = {c.lower(): c for c in rows[0]}
    taxon_col = cols.get("taxon")
    clust_col = cols.get("cluster")
    if not taxon_col or not clust_col:
        sys.exit("ERROR: expected Taxon and Cluster columns, found %s"
                 % list(rows[0]))

    members = defaultdict(list)
    for r in rows:
        members[r[clust_col]].append(r[taxon_col])

    # Reconcile PopPUNK's rewritten names back to the originals, then cross-check
    # against the rfile. A residual mismatch means PopPUNK genuinely dropped
    # genomes, which must not pass unnoticed.
    n_reconciled = 0
    if a.rfile:
        want = set()
        with open(a.rfile) as fh:
            for line in fh:
                if line.strip():
                    want.add(line.split("\t")[0].strip())
        by_sanitised = {sanitise(w): w for w in want}
        for cid, ms in members.items():
            fixed = []
            for t in ms:
                if t not in want and t in by_sanitised:
                    fixed.append(by_sanitised[t])
                    n_reconciled += 1
                else:
                    fixed.append(t)
            members[cid] = fixed
        got = {t for ms in members.values() for t in ms}
        missing, extra = want - got, got - want
        if n_reconciled:
            print("Reconciled %d PopPUNK-rewritten sample name(s) back to the "
                  "rfile originals" % n_reconciled, file=sys.stderr)
        if missing:
            print("WARNING: %d genome(s) in the rfile got no cluster, e.g. %s"
                  % (len(missing), sorted(missing)[:5]), file=sys.stderr)
        if extra:
            print("WARNING: %d clustered name(s) absent from the rfile even "
                  "after name reconciliation, e.g. %s"
                  % (len(extra), sorted(extra)[:5]), file=sys.stderr)

    kept, dropped = {}, {}
    for cid, ms in members.items():
        (kept if len(ms) >= a.min_cluster_size else dropped)[cid] = ms

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["cluster_id", "sample_id"])
        for cid in sorted(kept, key=lambda c: (-len(kept[c]), str(c))):
            for s in sorted(kept[cid]):
                w.writerow(["%s%s" % (a.prefix, cid), s])

    if a.excluded:
        with open(a.excluded, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["cluster_id", "n", "reason"])
            for cid in sorted(dropped, key=lambda c: (-len(dropped[c]), str(c))):
                w.writerow(["%s%s" % (a.prefix, cid), len(dropped[cid]),
                            "below_min_cluster_size_%d" % a.min_cluster_size])

    n_kept_g = sum(len(v) for v in kept.values())
    n_drop_g = sum(len(v) for v in dropped.values())
    sizes = sorted((len(v) for v in kept.values()), reverse=True)
    print("PopPUNK partition")
    print("  genomes clustered        : %d" % (n_kept_g + n_drop_g))
    print("  clusters total           : %d" % len(members))
    print("  min cluster size gate    : %d" % a.min_cluster_size)
    print("  clusters kept            : %d  (%d genomes)" % (len(kept), n_kept_g))
    print("  clusters excluded        : %d  (%d genomes)" % (len(dropped), n_drop_g))
    if sizes:
        print("  kept sizes largest->     : %s%s"
              % (sizes[:10], " ..." if len(sizes) > 10 else ""))
    print("  singletons in full fit   : %d"
          % sum(1 for v in members.values() if len(v) == 1))
    if not kept:
        sys.exit("ERROR: no cluster met the minimum size of %d -- nothing to "
                 "analyse. Lower --min_cluster_size or check the fit."
                 % a.min_cluster_size)


if __name__ == "__main__":
    main()
