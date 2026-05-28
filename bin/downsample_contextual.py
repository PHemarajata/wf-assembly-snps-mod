#!/usr/bin/env python3
"""
Systematically downsample contextual assemblies for the BP SNP/clustering workflow.

Strategy (hybrid: dereplicate, then per-country floor/cap):
  1. Map each contextual FASTA to its Country_Final / Subregion / collection date
     from the metadata TSV (unmatched files -> country "Unknown").
  2. Compute pairwise Mash distances among the contextual genomes (cheap, one-time)
     and collapse near-identical genomes WITHIN each country using single-linkage at
     --derep-threshold, keeping one representative per redundancy group. This breaks
     the single-linkage "chaining" that otherwise melts distinct clades into one
     mega-cluster downstream, and auto-scales how many genomes are kept per region to
     the actual diversity present there.
  3. If the representatives still exceed the budget (--target-total minus the CDC
     isolates, which are always kept in full), allocate the budget across countries
     with a per-country floor and cap so rare regions stay represented and an
     over-sampled region (e.g. Thailand) cannot dominate.
  4. Emit a `sample,file` samplesheet that the pipeline's --input already accepts,
     plus human-readable reports of what was kept and why.

Depends only on the Python standard library and the `mash` binary (>=2.x). Run it
once before the pipeline; the expensive --mash_sketch_size 100000 run then happens
only on the reduced set inside the normal workflow.
"""

import argparse
import csv
import os
import random
import re
import subprocess
import sys
from collections import defaultdict

FASTA_EXTS = (".fasta", ".fas", ".fna", ".fsa", ".fa")
UNKNOWN_COUNTRY = "Unknown"


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #
def normalize_name(path):
    """Basename without extension, matching bin/cluster_mash.py conventions."""
    return os.path.splitext(os.path.basename(path))[0]


def sanitize_sample(name):
    """Pipeline converts spaces to underscores; do it up front for stable IDs."""
    return name.strip().replace(" ", "_")


def parse_year(value):
    """Extract a 4-digit year for ordering; missing/garbage sorts last."""
    if not value:
        return 10**9
    digits = ""
    for ch in str(value).strip():
        if ch.isdigit():
            digits += ch
            if len(digits) == 4:
                break
        elif digits:
            break
    return int(digits) if len(digits) == 4 else 10**9


def is_known_subregion(sub):
    if sub is None:
        return False
    s = sub.strip().lower()
    return s not in ("", "unknown", "na", "n/a", "none")


def list_fastas(directory):
    out = []
    for fn in sorted(os.listdir(directory)):
        if fn.lower().endswith(FASTA_EXTS):
            out.append(os.path.join(directory, fn))
    return out


def dedup_by_canonical(files):
    """Collapse files that share a canonical accession key (same genome stored
    under two filename styles, e.g. 'GCA_x.2.fasta' and 'GCA_x_2.fasta').
    Prefer the underscore (TSV-style) spelling, then the lexicographically
    smallest name. Returns (kept_files, n_removed)."""
    by_key = defaultdict(list)
    for f in files:
        by_key[canonical_key(f)].append(f)
    kept, removed = [], 0
    for key, group in by_key.items():
        if len(group) == 1:
            kept.append(group[0])
            continue
        group.sort(key=lambda f: (os.path.basename(f).count("."), os.path.basename(f)))
        kept.append(group[0])
        removed += len(group) - 1
    return sorted(kept), removed


def eprint(*args):
    print(*args, file=sys.stderr, flush=True)


# --------------------------------------------------------------------------- #
# Union-Find (single-linkage clustering without scipy)
# --------------------------------------------------------------------------- #
class DisjointSet:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x):
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:  # path compression
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


# --------------------------------------------------------------------------- #
# Metadata
# --------------------------------------------------------------------------- #
def canonical_key(fasta_or_path):
    """Key that ignores the version *separator style*, so disk names like
    'GCA_000182195.1.fasta' and TSV names like 'GCA_000182195_1.fasta' collapse
    to the same key ('GCA_000182195_1'). Still version-specific."""
    return normalize_name(fasta_or_path).replace(".", "_")


def accession_key(fasta_or_path):
    """Version-ignoring key: strips a trailing '.N' or '_N' (NCBI assembly
    version, 1-3 digits) so 'GCA_000182195.1' and 'GCA_000182195_2' both map to
    'GCA_000182195'. Country/region metadata is version-independent, so this
    recovers genomes whose on-disk version differs from the metadata's."""
    return re.sub(r"[._]\d{1,3}$", "", normalize_name(fasta_or_path))


def load_metadata(tsv_path, country_col, subregion_col, fasta_col,
                  accession_col, sample_col):
    """Return a multi-key index {key: {'country','subregion','date'}}. Each row
    is indexed under the version-ignoring accession key AND the canonical key for
    every available identifier column, maximising match recovery."""
    idx = {}
    with open(tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if country_col not in reader.fieldnames or subregion_col not in reader.fieldnames:
            sys.exit(f"ERROR: '{country_col}'/'{subregion_col}' not found in "
                     f"{tsv_path}. Available: {reader.fieldnames}")
        id_cols = [c for c in (accession_col, fasta_col, sample_col)
                   if c and c in reader.fieldnames]
        if not id_cols:
            sys.exit(f"ERROR: none of the identifier columns "
                     f"({accession_col}, {fasta_col}, {sample_col}) found in {tsv_path}.")
        date_col = "final_collection_dates" if \
            "final_collection_dates" in reader.fieldnames else None
        for row in reader:
            info = {
                "country": (row.get(country_col) or "").strip() or UNKNOWN_COUNTRY,
                "subregion": (row.get(subregion_col) or "").strip(),
                "date": (row.get(date_col) or "").strip() if date_col else "",
            }
            for col in id_cols:
                v = (row.get(col) or "").strip()
                if not v:
                    continue
                idx.setdefault(accession_key(v), info)
                idx.setdefault(canonical_key(v), info)
    return idx


def annotate_contextual(files, idx):
    """Attach metadata to each contextual file, trying the version-ignoring
    accession key first, then the canonical key. Returns (records, n_matched)."""
    records = []
    matched = 0
    for path in files:
        info = idx.get(accession_key(path)) or idx.get(canonical_key(path))
        is_matched = info is not None
        if is_matched:
            matched += 1
        else:
            info = {"country": UNKNOWN_COUNTRY, "subregion": "", "date": ""}
        records.append({
            "sample": sanitize_sample(normalize_name(path)),
            "path": os.path.abspath(path),
            "name": os.path.basename(path),
            "country": info["country"],
            "subregion": info["subregion"],
            "date": info["date"],
            "matched": is_matched,
        })
    return records, matched


# --------------------------------------------------------------------------- #
# Mash
# --------------------------------------------------------------------------- #
def run_mash_triangle(files, out_path, sketch_size, kmer, threads):
    """Run `mash triangle` over a file list -> lower-triangular distances."""
    listfile = out_path + ".filelist.txt"
    with open(listfile, "w") as fh:
        for p in files:
            fh.write(os.path.abspath(p) + "\n")
    cmd = ["mash", "triangle", "-p", str(threads), "-k", str(kmer),
           "-s", str(sketch_size), "-l", listfile]
    eprint("Running:", " ".join(cmd), "->", out_path)
    with open(out_path, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)
    return out_path


def iter_mash_edges(path):
    """Yield (name_a, name_b, distance) from a mash triangle OR `mash dist` file.

    Auto-detects format: triangle starts with a bare integer count and lists each
    new name with distances to all earlier names; `mash dist` is long-format with
    id1, id2, dist in the first three tab-separated columns.
    """
    with open(path) as fh:
        first = fh.readline()
        while first.strip() == "":
            first = fh.readline()
            if first == "":
                return
        if first.strip().lstrip("\t").isdigit():
            # ---- mash triangle ----
            names = []
            for line in fh:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                name = normalize_name(parts[0])
                for j, val in enumerate(parts[1:]):
                    if val == "":
                        continue
                    try:
                        d = float(val)
                    except ValueError:
                        continue
                    yield name, names[j], d
                names.append(name)
        else:
            # ---- long `mash dist` format ----
            def emit(line):
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    return None
                a, b = normalize_name(parts[0]), normalize_name(parts[1])
                if a == b:
                    return None
                try:
                    return a, b, float(parts[2])
                except ValueError:
                    return None
            res = emit(first)
            if res:
                yield res
            for line in fh:
                res = emit(line)
                if res:
                    yield res


# --------------------------------------------------------------------------- #
# Dereplication
# --------------------------------------------------------------------------- #
def dereplicate(records, mash_path, threshold):
    """Single-linkage grouping WITHIN country at `threshold`. Returns
    group_id list aligned to records (index into records)."""
    idx_of = {r["sample"]: i for i, r in enumerate(records)}
    ds = DisjointSet(len(records))
    edges_used = 0
    for a, b, d in iter_mash_edges(mash_path):
        if d > threshold:
            continue
        ia, ib = idx_of.get(a), idx_of.get(b)
        if ia is None or ib is None:
            continue
        if records[ia]["country"] != records[ib]["country"]:
            continue  # never merge across countries
        ds.union(ia, ib)
        edges_used += 1
    groups = [ds.find(i) for i in range(len(records))]
    eprint(f"  dereplication: applied {edges_used} within-country edges "
           f"at threshold {threshold}")
    return groups


def representative_key(rec, rand_val):
    """Lower sorts first = preferred representative."""
    return (
        0 if is_known_subregion(rec["subregion"]) else 1,
        0 if rec["name"].upper().startswith("GCF") else 1,
        parse_year(rec["date"]),
        rand_val,
    )


def choose_representatives(records, groups, rng):
    """Pick one representative per group; return (rep_indices, group_members)."""
    rand_vals = {i: rng.random() for i in range(len(records))}
    members = defaultdict(list)
    for i, g in enumerate(groups):
        members[g].append(i)
    reps = []
    for g, idxs in members.items():
        best = min(idxs, key=lambda i: representative_key(records[i], rand_vals[i]))
        reps.append(best)
    return reps, members, rand_vals


# --------------------------------------------------------------------------- #
# Floor/cap allocation
# --------------------------------------------------------------------------- #
def allocate(reps_by_country, budget, min_per, max_per):
    """Decide how many reps to keep per country to land near `budget`.
    Returns {country: keep_count}."""
    countries = list(reps_by_country.keys())
    avail = {c: len(reps_by_country[c]) for c in countries}
    total = sum(avail.values())
    if total <= budget:
        return dict(avail)  # keep everything; dereplication sufficed

    cap = {c: min(avail[c], max_per) for c in countries}
    keep = {c: min(avail[c], min_per) for c in countries}  # floor first

    used = sum(keep.values())
    if used > budget:
        # Floors alone exceed budget: trim fairly, largest allocations first,
        # never below 1 while budget allows.
        order = sorted(countries, key=lambda c: keep[c], reverse=True)
        i = 0
        while used > budget:
            c = order[i % len(order)]
            if keep[c] > 1 or (keep[c] == 1 and used > budget):
                if keep[c] > 0:
                    keep[c] -= 1
                    used -= 1
            i += 1
            if all(keep[c] == 0 for c in countries):
                break
        return keep

    # Distribute the remaining budget proportionally to leftover capacity,
    # using largest-remainder passes until exhausted or no capacity remains.
    remaining = budget - used
    while remaining > 0:
        leftover = {c: cap[c] - keep[c] for c in countries if cap[c] - keep[c] > 0}
        if not leftover:
            break
        pool = sum(leftover.values())
        shares = {c: remaining * leftover[c] / pool for c in leftover}
        # integer floor allocation
        added = 0
        for c in leftover:
            give = min(int(shares[c]), leftover[c])
            keep[c] += give
            added += give
        remaining -= added
        if added == 0:
            # hand out the last few by largest fractional remainder
            ranked = sorted(leftover, key=lambda c: shares[c] - int(shares[c]),
                            reverse=True)
            for c in ranked:
                if remaining == 0:
                    break
                if cap[c] - keep[c] > 0:
                    keep[c] += 1
                    remaining -= 1
            break
    return keep


def select_within_country(reps, records, keep_count, rand_vals):
    """Pick `keep_count` of a country's representatives, spread across subregions
    then years for geographic/temporal coverage."""
    if keep_count >= len(reps):
        return list(reps)
    buckets = defaultdict(list)
    for i in reps:
        key = records[i]["subregion"].strip().lower() \
            if is_known_subregion(records[i]["subregion"]) else "~unknown"
        buckets[key].append(i)
    for key in buckets:
        buckets[key].sort(key=lambda i: (parse_year(records[i]["date"]), rand_vals[i]))
    # round-robin across subregion buckets (largest buckets first)
    ordered_keys = sorted(buckets, key=lambda k: len(buckets[k]), reverse=True)
    chosen, pos = [], 0
    while len(chosen) < keep_count:
        progressed = False
        for key in ordered_keys:
            if pos < len(buckets[key]):
                chosen.append(buckets[key][pos])
                progressed = True
                if len(chosen) == keep_count:
                    break
        pos += 1
        if not progressed:
            break
    return chosen


# --------------------------------------------------------------------------- #
# Sweep mode
# --------------------------------------------------------------------------- #
def run_sweep(records, mash_path, thresholds, out_path):
    idx_of = {r["sample"]: i for i, r in enumerate(records)}
    country = [r["country"] for r in records]
    sets = {t: DisjointSet(len(records)) for t in thresholds}
    tmax = max(thresholds)
    for a, b, d in iter_mash_edges(mash_path):
        if d > tmax:
            continue
        ia, ib = idx_of.get(a), idx_of.get(b)
        if ia is None or ib is None or country[ia] != country[ib]:
            continue
        for t in thresholds:
            if d <= t:
                sets[t].union(ia, ib)
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["derep_threshold", "n_contextual", "n_groups_kept",
                    "fraction_removed"])
        n = len(records)
        for t in thresholds:
            roots = {sets[t].find(i) for i in range(n)}
            kept = len(roots)
            w.writerow([t, n, kept, f"{(n - kept) / n:.4f}" if n else "0"])
    eprint(f"Wrote sweep report to {out_path}")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    p = argparse.ArgumentParser(
        description="Downsample contextual assemblies (dereplicate + geo balance).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--contextual-dir", required=True)
    p.add_argument("--cdc-dir", required=True,
                   help="Study isolates; always kept in full, never dereplicated.")
    p.add_argument("--metadata", required=True, help="Metadata TSV.")
    p.add_argument("--country-col", default="Country_Final")
    p.add_argument("--subregion-col", default="Subregion")
    p.add_argument("--fasta-col", default="FASTA_name")
    p.add_argument("--accession-col", default="Assembly Accession",
                   help="Accession column; matched version-ignored (primary key).")
    p.add_argument("--sample-col", default="sample_id",
                   help="Extra identifier column used as a matching fallback.")
    p.add_argument("--target-total", type=int, default=1000,
                   help="Target total genomes (CDC + contextual).")
    p.add_argument("--derep-threshold", type=float, default=0.0005,
                   help="Mash distance at/below which genomes are 'redundant'.")
    p.add_argument("--derep-sketch-size", type=int, default=50000,
                   help="Sketch size for the selection mash run.")
    p.add_argument("--kmer", type=int, default=21)
    p.add_argument("--min-per-country", type=int, default=3)
    p.add_argument("--max-per-country", type=int, default=200)
    p.add_argument("--threads", type=int, default=max(1, os.cpu_count() or 1))
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--unmatched", choices=["drop", "unknown"], default="drop",
                   help="How to treat contextual FASTAs with no metadata row. "
                        "'drop' = exclude them (treat the TSV as the authoritative, "
                        "cleaned keep-list); 'unknown' = keep them in an 'Unknown' "
                        "country bucket.")
    p.add_argument("--mash-dist", default=None,
                   help="Reuse a precomputed mash triangle or `mash dist` file.")
    p.add_argument("--sweep", default=None,
                   help="Comma-separated thresholds; report group counts and exit.")
    p.add_argument("--no-derep", action="store_true",
                   help="Skip mash/dereplication (every genome is its own group). "
                        "Useful for testing the allocation logic.")
    p.add_argument("--outdir", default="downsample_out")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(args.seed)

    # --- inputs ---
    contextual_files = list_fastas(args.contextual_dir)
    cdc_files = list_fastas(args.cdc_dir)
    if not contextual_files:
        sys.exit(f"ERROR: no FASTA files found in {args.contextual_dir}")
    eprint(f"Found {len(contextual_files)} contextual and {len(cdc_files)} CDC FASTAs")
    contextual_files, n_dups = dedup_by_canonical(contextual_files)
    if n_dups:
        eprint(f"Removed {n_dups} duplicate-filename copies "
               f"({len(contextual_files)} unique contextual genomes remain)")

    meta = load_metadata(args.metadata, args.country_col, args.subregion_col,
                         args.fasta_col, args.accession_col, args.sample_col)
    records, matched = annotate_contextual(contextual_files, meta)
    unmatched_records = [r for r in records if not r["matched"]]
    n_unmatched = len(unmatched_records)
    eprint(f"Matched metadata for {matched}/{len(records)} contextual files "
           f"({n_unmatched} have no metadata row)")
    if n_unmatched:
        # Always document the metadata-less genomes for provenance.
        dropped_path = os.path.join(args.outdir, "unmatched_no_metadata.tsv")
        with open(dropped_path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["filename", "accession_key", "path", "disposition"])
            disp = "dropped" if args.unmatched == "drop" else "kept_as_Unknown"
            for r in sorted(unmatched_records, key=lambda r: r["name"]):
                w.writerow([r["name"], accession_key(r["name"]), r["path"], disp])
        eprint(f"  documented {n_unmatched} metadata-less genomes -> {dropped_path}")
        if args.unmatched == "drop":
            records = [r for r in records if r["matched"]]
            eprint(f"  --unmatched=drop: excluded {n_unmatched} files lacking "
                   f"metadata; {len(records)} contextual files remain")
        else:
            eprint(f"  --unmatched=unknown: keeping {n_unmatched} files in the "
                   f"'{UNKNOWN_COUNTRY}' bucket")
    if not records:
        sys.exit("ERROR: no contextual files left after metadata matching.")

    # --- mash distances ---
    mash_path = None
    if not args.no_derep:
        if args.mash_dist:
            mash_path = args.mash_dist
            eprint(f"Reusing distances from {mash_path}")
        else:
            mash_path = run_mash_triangle(
                [r["path"] for r in records],
                os.path.join(args.outdir, "mash_distances.tsv"),
                args.derep_sketch_size, args.kmer, args.threads)

    # --- sweep mode: report and exit ---
    if args.sweep:
        if mash_path is None:
            sys.exit("ERROR: --sweep needs mash distances (don't combine with --no-derep)")
        thresholds = sorted(float(t) for t in args.sweep.split(","))
        run_sweep(records, mash_path, thresholds,
                  os.path.join(args.outdir, "derep_sweep.tsv"))
        return

    # --- dereplicate ---
    if args.no_derep:
        groups = list(range(len(records)))
    else:
        groups = dereplicate(records, mash_path, args.derep_threshold)
    reps, members, rand_vals = choose_representatives(records, groups, rng)
    rep_set = set(reps)
    eprint(f"Dereplication: {len(records)} contextual -> {len(reps)} representatives")

    # --- allocate budget across countries ---
    reps_by_country = defaultdict(list)
    for i in reps:
        reps_by_country[records[i]["country"]].append(i)
    budget = args.target_total - len(cdc_files)
    if budget < 0:
        eprint(f"WARNING: target-total ({args.target_total}) < CDC count "
               f"({len(cdc_files)}); keeping only CDC isolates.")
        budget = 0

    keep_counts = allocate(reps_by_country, budget,
                           args.min_per_country, args.max_per_country)
    kept_idx = set()
    for c, reps_c in reps_by_country.items():
        chosen = select_within_country(reps_c, records, keep_counts.get(c, 0),
                                        rand_vals)
        kept_idx.update(chosen)
    eprint(f"Kept {len(kept_idx)} contextual representatives "
           f"(budget {budget}) + {len(cdc_files)} CDC")

    # --- write samplesheet ---
    samplesheet = os.path.join(args.outdir, "samplesheet.csv")
    with open(samplesheet, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["sample", "file"])
        for path in cdc_files:
            w.writerow([sanitize_sample(normalize_name(path)), os.path.abspath(path)])
        for i in sorted(kept_idx, key=lambda i: records[i]["sample"]):
            w.writerow([records[i]["sample"], records[i]["path"]])
    eprint(f"Wrote samplesheet ({len(cdc_files) + len(kept_idx)} rows) -> {samplesheet}")

    # --- per-genome selection report ---
    report = os.path.join(args.outdir, "selection_report.tsv")
    rep_of_group = {g: min(idxs, key=lambda i: representative_key(records[i], rand_vals[i]))
                    for g, idxs in members.items()}
    with open(report, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", "source", "country", "subregion", "date",
                    "group_id", "role", "kept", "reason"])
        for path in cdc_files:
            w.writerow([sanitize_sample(normalize_name(path)), "cdc", "", "", "",
                        "", "cdc", "yes", "study_isolate"])
        for i, rec in enumerate(records):
            g = groups[i]
            is_rep = i in rep_set
            kept = i in kept_idx
            if is_rep:
                role = "representative"
                reason = "kept" if kept else "dropped_by_country_cap"
            else:
                role = "redundant"
                rep_i = rep_of_group[g]
                reason = f"redundant_to:{records[rep_i]['sample']}"
            w.writerow([rec["sample"], "contextual", rec["country"], rec["subregion"],
                        rec["date"], g, role, "yes" if kept else "no", reason])
    eprint(f"Wrote per-genome report -> {report}")

    # --- per-country summary ---
    summary = os.path.join(args.outdir, "country_summary.tsv")
    n_input = defaultdict(int)
    for r in records:
        n_input[r["country"]] += 1
    n_final = defaultdict(int)
    for i in kept_idx:
        n_final[records[i]["country"]] += 1
    with open(summary, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["country", "n_input", "n_after_derep", "n_final"])
        for c in sorted(n_input, key=lambda c: n_input[c], reverse=True):
            w.writerow([c, n_input[c], len(reps_by_country.get(c, [])),
                        n_final.get(c, 0)])
    eprint(f"Wrote country summary -> {summary}")


if __name__ == "__main__":
    main()
