#!/usr/bin/env python3
"""
Pool Gubbins per-branch SNP counts into a per-unit r/m, with the external
mapping reference's branches excluded.

WHY THE EXCLUSION IS NOT OPTIONAL
---------------------------------
The pipeline keeps the mapping reference as a taxon in the Gubbins input on
purpose: it keeps the alignment full-length and the invariant-site counts
honest. But Gubbins then reconstructs substitutions along the branch leading to
that reference, and because the reference sits OUTSIDE the population by
construction, that branch is enormous. Measured on strain_18_L1_1:

    (Reference:3859.7,( ...7 real genomes... )Node_6:3774.9)Node_7:0.0

every real genome sits on a branch of 4-52; the reference's is 3,860. Those
substitutions are divergence between the population and an outgroup, not
evolution within the population, and Gubbins scores nearly all of them as
"outside recombination" -- so they land in r/m's denominator:

    strain_18_L1_1  including reference branches  r/m = 0.42
                    excluding them                r/m = 8.73
                    manual analysis (7 taxa)      r/m = 9.14

Across the 2,070-genome run, 52% of ALL outside-recombination SNPs came from
reference branches, and the pooled median r/m moved 1.85 -> 6.30. The effect
scales with reference distance because the branch length IS the reference
distance, which is why it first masqueraded as a caller x distance interaction.

WHICH BRANCHES ARE DROPPED
--------------------------
Gubbins emits an unrooted tree written with an arbitrary root.

  * If the `Reference` leaf is a child of that root, its divergence is SPLIT
    between the leaf and the sibling clade's stem -- above, 3859.7 and 3774.9,
    two halves of one quantity. Dropping only the leaf would leave half the
    inflation behind, so BOTH children of the root go.
  * Otherwise the root sits inside the ingroup and the whole outgroup
    divergence is on the `Reference` leaf alone, so only that branch goes. The
    unit is still noted, because where the root landed is worth seeing.

REPORTING THE LONGEST REMAINING BRANCH
--------------------------------------
`max_kept_branch` / `max_kept_branch_len` name the longest branch that SURVIVES
the exclusion. A unit whose longest surviving branch is orders of magnitude
above its siblings has a divergent MEMBER, not a reference artefact, and its
r/m is depressed for a reason this correction does not address. That is the raw
per-item view; do not diagnose a low r/m without looking at it.

Ported from the validated `exclude_reference_branches_bp.py`, which produced
RM_RESULTS_L1_CORRECTED.tsv. The drop rule is byte-for-byte the same rule.
"""

import argparse
import collections
import csv
import glob
import os
import re
import sys


# ----------------------------------------------------------------------------
# Newick handling -- deliberately string-level, so this needs no dependencies
# beyond the stdlib and runs in the plain python container.
# ----------------------------------------------------------------------------

def _split_top_level(inner):
    """Split a Newick child list on commas at nesting depth zero."""
    parts, depth, cur = [], 0, []
    for ch in inner:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    parts.append("".join(cur))
    return parts


def root_children(newick):
    """
    Labels of the children of the root.

    Each child's label is its trailing token: a leaf name, or the Node_N label
    a clade carries after its closing bracket.
    """
    s = newick.strip().rstrip(";").strip()
    if not s.startswith("("):
        return []
    depth = 0
    inner = None
    for i, ch in enumerate(s):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                inner = s[1:i]
                break
    if inner is None:
        return []

    labels = []
    for p in _split_top_level(inner):
        p = p.strip()
        tail = p[p.rfind(")") + 1:] if ")" in p else p
        lab = tail.split(":")[0].strip()
        if lab:
            labels.append(lab)
    return labels


def branch_lengths(newick):
    """{label: branch_length} for every labelled branch in the tree."""
    out = {}
    for lab, blen in re.findall(r"([A-Za-z0-9_.\-]+):(-?[0-9.eE+]+)", newick):
        try:
            out[lab] = float(blen)
        except ValueError:
            continue
    return out


def branches_to_drop(tree_path, ref_name):
    """(set_of_labels_to_drop, note)"""
    try:
        with open(tree_path) as fh:
            newick = fh.read()
    except OSError:
        return set(), "no tree", {}
    lengths = branch_lengths(newick)
    if ref_name not in newick:
        return set(), "no %s taxon in tree" % ref_name, lengths
    kids = root_children(newick)
    if ref_name in kids:
        return set(kids), "reference at root; dropped both root children", lengths
    return {ref_name}, "REFERENCE NOT AT ROOT -- only its own branch dropped", lengths


# ----------------------------------------------------------------------------

def pooled(path, drop):
    """(inside, outside, dropped_inside, dropped_outside) from one stats file."""
    inside = outside = dropped_in = dropped_out = 0.0
    try:
        fh = open(path)
    except OSError:
        return inside, outside, dropped_in, dropped_out
    with fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                i = float(row["Number of SNPs Inside Recombinations"])
                o = float(row["Number of SNPs Outside Recombinations"])
            except (KeyError, TypeError, ValueError):
                continue
            if row.get("Node") in drop:
                dropped_in += i
                dropped_out += o
            else:
                inside += i
                outside += o
    return inside, outside, dropped_in, dropped_out


def unit_of(replicon_unit_id):
    """
    SPLIT_REFERENCE_REPLICONS names each replicon-unit `<cluster_id>__<replicon>`.
    r/m pools over a cluster's replicons, so the cluster id is everything before
    the first `__`. Without replicon splitting there is no `__` and the id is
    already the unit.
    """
    return replicon_unit_id.split("__", 1)[0]


def sort_key(unit):
    """Natural order on strain_<a>_L1_<b> ids, falling back to plain string."""
    m = re.match(r"^\D*(\d+).*?_L1_(\d+)$", unit)
    if m:
        return (0, int(m.group(1)), int(m.group(2)), unit)
    return (1, 0, 0, unit)


def median(vals):
    s = sorted(vals)
    n = len(s)
    if not n:
        return float("nan")
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stats-dir", required=True,
                    help="directory of <unit>.per_branch_statistics.csv")
    ap.add_argument("--trees-dir", required=True,
                    help="directory of <unit>.node_labelled.final_tree.tre")
    ap.add_argument("--clusters-tsv", default=None,
                    help="clusters.tsv (cluster_id, sample_id) for the n column")
    ap.add_argument("--reference-taxon", default="Reference",
                    help="name Gubbins gives the mapping reference (default: Reference)")
    ap.add_argument("--output", required=True)
    a = ap.parse_args()

    sizes = collections.Counter()
    if a.clusters_tsv and os.path.isfile(a.clusters_tsv):
        with open(a.clusters_tsv) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                cid = row.get("cluster_id")
                if cid:
                    sizes[cid] += 1

    per = collections.defaultdict(dict)
    notes = collections.defaultdict(list)
    longest = {}

    stats_files = sorted(glob.glob(os.path.join(a.stats_dir, "*.per_branch_statistics.csv")))
    if not stats_files:
        # The file is `<unit>.per_branch_statistics.csv`, never a bare
        # `per_branch_statistics.csv` -- a glob that misses looks exactly like
        # "Gubbins never produced it", which is a far scarier conclusion.
        sys.stderr.write("ERROR: no *.per_branch_statistics.csv under %s\n" % a.stats_dir)
        return 1

    for path in stats_files:
        ru = os.path.basename(path)[:-len(".per_branch_statistics.csv")]
        if os.path.getsize(path) == 0:
            # Gubbins was skipped for this replicon-unit (too few sequences, or
            # an empty alignment); it contributes no SNPs, not a zero r/m.
            notes[unit_of(ru)].append("skipped: empty per-branch statistics")
            continue

        tre = os.path.join(a.trees_dir, "%s.node_labelled.final_tree.tre" % ru)
        drop, note, lengths = branches_to_drop(tre, a.reference_taxon)
        notes[unit_of(ru)].append(note)
        per[unit_of(ru)][ru] = pooled(path, drop) + (sorted(drop),)

        kept = {k: v for k, v in lengths.items() if k not in drop}
        if kept:
            lab = max(kept, key=kept.get)
            u = unit_of(ru)
            if u not in longest or kept[lab] > longest[u][1]:
                longest[u] = (lab, kept[lab])

    rows, flagged = [], []
    for unit in sorted(per, key=sort_key):
        reps = per[unit]
        i = sum(v[0] for v in reps.values())
        o = sum(v[1] for v in reps.values())
        di = sum(v[2] for v in reps.values())
        do = sum(v[3] for v in reps.values())
        if any("NOT AT ROOT" in n for n in notes[unit]):
            flagged.append(unit)
        lab, blen = longest.get(unit, ("", ""))
        rows.append({
            "unit": unit,
            "n": sizes.get(unit, ""),
            "n_replicons": len(reps),
            "rm_corrected": "%.4f" % (i / o) if o else "NA",
            "rm_uncorrected": "%.4f" % ((i + di) / (o + do)) if (o + do) else "NA",
            "snps_in_recomb": int(i),
            "snps_outside": int(o),
            "ref_branch_snps_inside": int(di),
            "ref_branch_snps_outside": int(do),
            "dropped_branches": ";".join(sorted({b for v in reps.values() for b in v[4]})),
            "max_kept_branch": lab,
            "max_kept_branch_len": ("%.4f" % blen) if blen != "" else "",
            "note": "; ".join(sorted(set(notes[unit]))),
        })

    cols = ["unit", "n", "n_replicons", "rm_corrected", "rm_uncorrected",
            "snps_in_recomb", "snps_outside", "ref_branch_snps_inside",
            "ref_branch_snps_outside", "dropped_branches", "max_kept_branch",
            "max_kept_branch_len", "note"]
    with open(a.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow(row)

    vals = [float(r["rm_corrected"]) for r in rows if r["rm_corrected"] != "NA"]
    unc = [float(r["rm_uncorrected"]) for r in rows if r["rm_uncorrected"] != "NA"]
    print("units pooled: %d (from %d replicon-units)"
          % (len(rows), sum(r["n_replicons"] for r in rows)))
    if vals:
        print("  r/m corrected    median %.2f   range %.2f-%.2f"
              % (median(vals), min(vals), max(vals)))
    if unc:
        print("  r/m uncorrected  median %.2f   range %.2f-%.2f"
              % (median(unc), min(unc), max(unc)))
    tot_out = sum(r["snps_outside"] + r["ref_branch_snps_outside"] for r in rows)
    ref_out = sum(r["ref_branch_snps_outside"] for r in rows)
    if tot_out:
        print("  outside-recombination SNPs on reference branches: %d of %d (%.1f%%)"
              % (ref_out, tot_out, 100.0 * ref_out / tot_out))
    if flagged:
        print("  reference NOT at root in %d unit(s): %s%s"
              % (len(flagged), ", ".join(flagged[:5]),
                 " ..." if len(flagged) > 5 else ""))
    print("wrote %s" % a.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
