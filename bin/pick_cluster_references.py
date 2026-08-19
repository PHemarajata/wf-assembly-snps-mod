#!/usr/bin/env python3
"""
Choose a per-cluster mapping reference, preferring a COMPLETE member and
borrowing the nearest complete genome when a cluster has none.

Emits the `cluster_id<TAB>reference_path` contract that --cluster_references
accepts, so the clustering front end (Mash or PopPUNK) and the reference choice
are decoupled from everything downstream.

WHY A SEPARATE SELECTION FROM THE MEDOID REPRESENTATIVE. The representative used
for the BACKBONE and the reference used for MAPPING are different objects with
different requirements, and conflating them is what previously picked a
135-contig draft as a cluster reference. Two hard constraints:

  1. Gubbins cannot use a multi-contig reference at all, and --split_replicons
     needs one FASTA per replicon. On this collection 2,541 of 2,802 genomes
     (90.7%) are drafts, so a plain medoid will almost always be unusable.
  2. The reference must be inside, or very close to, its own cluster. The
     measured reference effect is +28% SNPs in a diffuse cluster and +630% in a
     tight one, with ~87% of calls false in the latter.

SELECTION LOGIC, and the distinction matters:
  * Completeness (<= --max-contigs) is a GATE, not a ranking. Pass/fail.
  * Among genomes that pass the gate, rank by CENTRALITY -- mean, then max,
    Mash distance to cluster members. Distance to members is what drives
    mismapping and callable fraction.
  * N50 and total length are NOT ranking criteria once the gate is passed. A
    longer reference merely adds positions other members cannot fill, raising
    missing data and risking --filter-percentage exclusions.
  * A cluster with no complete member BORROWS the complete genome (from
    anywhere in the collection) with the lowest mean Mash distance to its
    members. The manual analysis resolved 45 units this way: 12 internal,
    33 borrowed.

Ported from pick_cluster_references_bp.py. The distance matrix is read through
the pipeline's own bin/mash_matrix_io.py rather than a private parser, so this
accepts every format the rest of the pipeline does -- lower-triangular Phylip
(what MASH_TRIANGLE emits), labelled square TSV (what the archived
mash_matrix_2802.tsv is), and raw `mash dist` edges. A private phylip-only
reader silently rejected the square TSV.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mash_matrix_io import read_any  # noqa: E402


def load_matrix(path):
    """(matrix, index) where index maps every plausible spelling of a name.

    The matrix may carry filenames (GCA_x_1.fasta) while membership carries bare
    sample ids, so index under each variant or the join drops every taxon.
    """
    labels, matrix = read_any(path)
    idx = {}
    for i, nm in enumerate(labels):
        nm = str(nm)
        for v in {nm,
                  nm.rsplit(".", 1)[0] if "." in nm else nm,
                  nm.replace(".fasta", "").replace(".fa", "").replace(".fna", "")}:
            idx.setdefault(v, i)
    return matrix, idx


def contig_stats(path):
    """(n_contigs, n50, total_len) in one pass, sequences not retained."""
    lens, cur = [], 0
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    if cur:
                        lens.append(cur)
                    cur = 0
                else:
                    cur += len(line.strip())
        if cur:
            lens.append(cur)
    except OSError:
        return None
    if not lens:
        return None
    total = sum(lens)
    lens.sort(reverse=True)
    acc, n50 = 0, lens[-1]
    for L in lens:
        acc += L
        if acc >= total / 2:
            n50 = L
            break
    return (len(lens), n50, total)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--membership", required=True,
                    help="clusters.tsv: cluster_id<TAB>sample_id")
    ap.add_argument("--fasta-dir", required=True)
    ap.add_argument("--phylip", required=True,
                    help="Mash distance matrix (MASH_TRIANGLE output)")
    ap.add_argument("--id-col", default="sample_id")
    ap.add_argument("--cluster-col", default="cluster_id")
    ap.add_argument("--max-contigs", type=int, default=4)
    ap.add_argument("--min-size", type=int, default=2)
    ap.add_argument("--reference-pool",
                    help="Optional file of allowed borrow candidates, one name "
                         "or path per line. WITHOUT it the pool is every "
                         "complete genome in --fasta-dir, which finds closer "
                         "references but diverges from a curated shortlist: on "
                         "this collection the open pool (192 complete genomes) "
                         "picked a reference 35-45%% closer in mean Mash "
                         "distance than the manual analysis's 14-reference "
                         "shortlist for 16 of 36 units. Closer is better for "
                         "mismapping, but changing the reference changes r/m, "
                         "so pass the shortlist when results must stay "
                         "comparable to existing baselines.")
    ap.add_argument("--reference-blocklist",
                    help="File of accession prefixes never to use as a mapping "
                         "reference, one per line, '#' comments allowed. Applies "
                         "to internal picks and borrow candidates alike. NO "
                         "REFERENCE IS CURRENTLY KNOWN TO NEED THIS. Three were "
                         "once listed here after every cluster mapped against "
                         "them failed Gubbins with 'Unable to fit model to data' "
                         "while 23 others gave 28/28 successes -- but the cause "
                         "was raxmlHPC v8 segfaulting at a -n run id of 128+ "
                         "characters, which Gubbins builds from the unit name, "
                         "which comes from the reference's FASTA defline. Those "
                         "three just had long filenames; holding the alignment "
                         "bytes identical and shortening only the basename made "
                         "all 12 previously-failing replicon-units succeed. A "
                         "reference that 'fails' is a claim about the genome: "
                         "isolate the variable before believing it.")
    ap.add_argument("--out", default="cluster_references.tsv")
    ap.add_argument("--report", default="reference_selection.tsv")
    a = ap.parse_args()

    blocked = []
    if a.reference_blocklist:
        with open(a.reference_blocklist) as fh:
            for line in fh:
                line = line.split("#", 1)[0].strip()
                if line:
                    blocked.append(line)
        sys.stderr.write("reference blocklist: %d prefixes -- %s\n"
                         % (len(blocked), ", ".join(blocked)))

    def is_blocked(stem):
        return any(stem.startswith(p) for p in blocked)

    with open(a.membership, newline="") as fh:
        head = fh.readline()
        fh.seek(0)
        d = "\t" if "\t" in head else ","
        mem = list(csv.DictReader(fh, delimiter=d))
    clusters = defaultdict(list)
    for r in mem:
        cid = (r.get(a.cluster_col) or "").strip()
        sid = (r.get(a.id_col) or "").strip()
        if cid and sid:
            clusters[cid].append(sid)
    if not clusters:
        sys.exit("ERROR: no cluster membership parsed from %s" % a.membership)

    listing = {}
    for fn in os.listdir(a.fasta_dir):
        if fn.endswith((".fasta", ".fa", ".fna")):
            listing[fn.rsplit(".", 1)[0]] = os.path.abspath(
                os.path.join(a.fasta_dir, fn))

    _mat, _idx = load_matrix(a.phylip)

    def centrality(members, cand):
        """(mean, max) Mash distance from cand to all members."""
        if cand not in _idx:
            return None
        i = _idx[cand]
        ds = [float(_mat[i][_idx[m]]) for m in members if m in _idx]
        return (sum(ds) / len(ds), max(ds)) if ds else None

    cache = {}

    def stats(sid):
        if sid in cache:
            return cache[sid]
        p = listing.get(sid)
        if p is None:
            for stem, q in listing.items():
                if stem.startswith(sid):
                    p = q
                    break
        cache[sid] = (contig_stats(p) if p else None, p)
        return cache[sid]

    # Every complete genome available is a borrow candidate -- the whole input
    # collection, NOT just genomes that landed in an analysed cluster. This
    # distinction is large: on the 2,802-genome collection 261 genomes are
    # complete, but only 26 of them sit inside the 1,051 clustered ones. Drawing
    # the pool from cluster members alone would hand most clusters a needlessly
    # distant reference, which is exactly the effect (+630% SNPs in a tight
    # cluster) the whole selection exists to avoid.
    all_ids = sorted({s for ms in clusters.values() for s in ms})
    allowed = None
    if a.reference_pool:
        allowed = set()
        with open(a.reference_pool) as fh:
            for line in fh:
                nm = line.strip()
                if nm:
                    allowed.add(os.path.splitext(os.path.basename(nm))[0])
    complete_pool = []
    for stem in sorted(listing):
        if allowed is not None and stem not in allowed:
            continue
        if is_blocked(stem):
            continue
        v = contig_stats(listing[stem])
        if v and v[0] <= a.max_contigs:
            cache.setdefault(stem, (v, listing[stem]))
            complete_pool.append(stem)
    if allowed is not None and not complete_pool:
        sys.exit("ERROR: --reference-pool matched no complete genome in %s"
                 % a.fasta_dir)

    chosen, report = {}, []
    for cid, sids in sorted(clusters.items()):
        if len(sids) < a.min_size:
            report.append((cid, len(sids), "", "", "", "SKIPPED_below_min_size", "", ""))
            continue

        internal_complete = []
        for s in sids:
            if is_blocked(s):
                continue
            v, p = stats(s)
            if v and v[0] <= a.max_contigs and p:
                internal_complete.append((s, v[0], v[1], p))

        if internal_complete:
            ranked = []
            for sid_, nc_, n50_, path_ in internal_complete:
                c = centrality(sids, sid_)
                if c is not None:
                    ranked.append((c[0], c[1], nc_, sid_, path_))
            if ranked:
                ranked.sort()                      # mean distance, then max
                mean_d, max_d, nc, sid, path = ranked[0]
                chosen[cid] = path
                report.append((cid, len(sids), sid, nc, len(internal_complete),
                               "internal", "%.6f" % mean_d, "%.6f" % max_d))
                continue
            # No centrality available (taxon absent from the matrix): still
            # prefer a complete internal member over borrowing.
            internal_complete.sort(key=lambda t: (t[1], -t[2]))
            sid, nc, n50, path = internal_complete[0]
            chosen[cid] = path
            report.append((cid, len(sids), sid, nc, len(internal_complete),
                           "internal_no_centrality", "", ""))
            continue

        # ---- borrow: nearest complete genome anywhere in the collection ----
        best, bestd, bestmax = None, None, None
        for cand in complete_pool:
            if cand in sids:
                continue
            c = centrality(sids, cand)
            if c is None:
                continue
            if bestd is None or c[0] < bestd:
                best, bestd, bestmax = cand, c[0], c[1]
        if best:
            _v, p = stats(best)
            chosen[cid] = p
            report.append((cid, len(sids), best, _v[0], 0, "borrowed",
                           "%.6f" % bestd, "%.6f" % bestmax))
        else:
            report.append((cid, len(sids), "", "", 0, "NO_REFERENCE", "", ""))

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["cluster_id", "reference_path"])
        for cid in sorted(chosen):
            w.writerow([cid, chosen[cid]])

    with open(a.report, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["cluster_id", "n", "reference", "ref_contigs",
                    "n_complete_members", "source", "mean_mash", "max_mash"])
        for row in report:
            w.writerow(row)

    counts = defaultdict(int)
    for r in report:
        counts[r[5]] += 1
    print("=" * 70)
    print("PER-CLUSTER REFERENCE SELECTION  (<= %d contigs required)"
          % a.max_contigs)
    print("=" * 70)
    print("  clusters considered      : %d" % len(clusters))
    print("  complete genomes in pool : %d / %d" % (len(complete_pool), len(all_ids)))
    for k in sorted(counts):
        print("  %-24s : %d" % (k, counts[k]))
    print("  references written       : %d" % len(chosen))

    missing = [r[0] for r in report if r[5] == "NO_REFERENCE"]
    if missing:
        print("\nERROR: %d cluster(s) could not be given a reference: %s"
              % (len(missing), missing[:10]), file=sys.stderr)
        sys.exit("No usable reference for %d cluster(s). Every cluster needs one "
                 "or the per-cluster chain cannot run." % len(missing))


if __name__ == "__main__":
    main()
