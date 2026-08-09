#!/usr/bin/env python3
"""
Roll up per-cluster Gubbins diagnostics + IQ-TREE outputs into a single
cluster_phylogeny_summary.csv.

Inputs (all positional or via flags):
  --clusters-tsv    clusters.tsv (cluster_id, sample_id) from CLUSTER_GENOMES
  --diagnostics-dir directory containing <cluster_id>.diagnostics.log files
  --trees-dir       directory containing <cluster_id>.final.treefile and .final.iqtree
  --filtered-dir    directory containing <cluster_id>.filtered_polymorphic_sites.fasta
  --recomb-dir      directory containing <cluster_id>.recombination_predictions.gff
  --output          path to cluster_phylogeny_summary.csv

Columns produced:
  cluster_id, n_isolates, alignment_total_chars, alignment_length_per_seq,
  seq_count_in_alignment, starting_tree_present, gubbins_status,
  gubbins_exit_code, n_recombination_blocks, n_filtered_polymorphic_sites,
  iqtree_status, treefile_size_bytes, iqtree_log_size_bytes,
  confidence_tier, notes

Note on the two alignment size columns:
  GUBBINS_CLUSTER's diagnostics log reports "Alignment length (non-header
  chars)", which is the TOTAL sequence character count summed over every
  record in the input alignment -- not the number of alignment columns. That
  raw total is reported as `alignment_total_chars`; for a 49-taxon cluster of
  a ~7 Mb organism it is ~340 Mb, which is not an alignment length in any
  useful sense. `alignment_length_per_seq` is the per-sequence alignment
  length (total / seq_count_in_alignment), which is the number readers
  usually want. If the total is not an exact multiple of the sequence count
  the alignment is ragged; the per-seq value is then left blank and a
  `ragged_input_alignment` note is emitted.
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path


def parse_diagnostics(path):
    """Parse a single .diagnostics.log file produced by GUBBINS_CLUSTER."""
    info = {
        "alignment_total_chars": "",
        "alignment_length_per_seq": "",
        "seq_count_in_alignment": "",
        "starting_tree_present": "false",
        "gubbins_status": "unknown",
        "gubbins_exit_code": "",
        "notes": [],
    }
    if not path.is_file():
        info["gubbins_status"] = "no_diagnostics"
        return info

    text = path.read_text(errors="replace").splitlines()
    skipped_too_few = False
    failed = False
    for line in text:
        m = re.match(r"Sequence count:\s*(\d+)", line)
        if m:
            info["seq_count_in_alignment"] = m.group(1)
        # This is the total character count over ALL records, not the number
        # of alignment columns -- see the module docstring.
        m = re.match(r"Alignment length \(non-header chars\):\s*(\d+)", line)
        if m:
            info["alignment_total_chars"] = m.group(1)
        m = re.match(r"Starting tree present:.*\((\d+) bytes\)", line)
        if m:
            info["starting_tree_present"] = "true"
        if line.startswith("Starting tree: (none provided)"):
            info["starting_tree_present"] = "false"
        m = re.match(r"Gubbins exit code:\s*(\S+)", line)
        if m:
            info["gubbins_exit_code"] = m.group(1)
        if "skipping Gubbins" in line:
            skipped_too_few = True
        if line.startswith("ERROR: Gubbins failed"):
            failed = True
        if line.startswith("ERROR: Alignment file missing or empty"):
            info["notes"].append("input_alignment_empty")

    # Derive the per-sequence alignment length from the raw total. A rectangular
    # alignment divides exactly; anything else is ragged and gets flagged rather
    # than reported as a rounded-off length.
    total = info["alignment_total_chars"]
    seqs = info["seq_count_in_alignment"]
    if total.isdigit() and seqs.isdigit() and int(seqs) > 0:
        total_i, seqs_i = int(total), int(seqs)
        if total_i % seqs_i == 0:
            info["alignment_length_per_seq"] = str(total_i // seqs_i)
        else:
            info["notes"].append("ragged_input_alignment")

    code = info["gubbins_exit_code"]
    if skipped_too_few:
        info["gubbins_status"] = "skipped_low_information"
    elif failed:
        info["gubbins_status"] = "failed_fallback_to_sanitized_alignment"
    elif code == "0":
        info["gubbins_status"] = "completed"
    elif code in ("", "NA"):
        info["gubbins_status"] = "unknown"
    else:
        info["gubbins_status"] = f"failed_exit_{code}"

    return info


def count_fasta_records(path):
    if not path.is_file() or path.stat().st_size == 0:
        return 0
    n = 0
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def alignment_n_columns(path):
    """Return the alignment column count by reading the first record's sequence length."""
    if not path.is_file() or path.stat().st_size == 0:
        return 0
    cols = 0
    in_first = False
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if in_first:
                    break
                in_first = True
                continue
            if in_first:
                cols += len(line.strip())
    return cols


def count_gff_recombination_blocks(path):
    if not path.is_file() or path.stat().st_size == 0:
        return 0
    n = 0
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            n += 1
    return n


def determine_iqtree_status(treefile, iqtree_log, seq_count):
    """Classify the final IQ-TREE outcome for a cluster."""
    if treefile is None or not treefile.is_file() or treefile.stat().st_size == 0:
        return "failed_no_tree"
    text = treefile.read_text(errors="replace").strip()
    if not text:
        return "failed_empty_tree"
    if seq_count is not None and seq_count < 3:
        return "skipped_low_information"
    if iqtree_log is None or not iqtree_log.is_file() or iqtree_log.stat().st_size == 0:
        return "minimal_tree_no_log"
    return "completed"


def confidence_tier(n_isolates, gubbins_status, iqtree_status, n_filtered_sites):
    """Map status flags to handoff §11 tiers."""
    if n_isolates <= 1:
        return "Tier4_singleton_or_uninterpretable"
    if iqtree_status not in ("completed", "minimal_tree_no_log"):
        return "Tier4_no_meaningful_tree"
    if gubbins_status == "completed" and iqtree_status == "completed" and n_isolates >= 6 and n_filtered_sites >= 10:
        return "Tier1_high_confidence"
    if gubbins_status in ("completed", "skipped_low_information", "failed_fallback_to_sanitized_alignment"):
        if n_isolates >= 4 and iqtree_status == "completed":
            return "Tier2_moderate_confidence"
    return "Tier3_low_confidence_exploratory"


def load_cluster_sizes(clusters_tsv):
    sizes = defaultdict(int)
    if not clusters_tsv.is_file():
        return sizes
    with clusters_tsv.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cid = row.get("cluster_id")
            if cid:
                sizes[cid] += 1
    return sizes


def discover_cluster_ids(diagnostics_dir, trees_dir):
    ids = set()
    if diagnostics_dir.is_dir():
        for p in diagnostics_dir.glob("*.diagnostics.log"):
            ids.add(p.name.replace(".diagnostics.log", ""))
    if trees_dir.is_dir():
        for p in trees_dir.glob("*.final.treefile"):
            ids.add(p.name.replace(".final.treefile", ""))
    return sorted(ids)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters-tsv", type=Path, required=True)
    ap.add_argument("--diagnostics-dir", type=Path, required=True)
    ap.add_argument("--trees-dir", type=Path, required=True)
    ap.add_argument("--filtered-dir", type=Path, required=True)
    ap.add_argument("--recomb-dir", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    cluster_sizes = load_cluster_sizes(args.clusters_tsv)
    cluster_ids = discover_cluster_ids(args.diagnostics_dir, args.trees_dir)

    # Fall back: include any cluster mentioned in clusters.tsv even if it has no diagnostic/tree
    for cid in cluster_sizes:
        if cid not in cluster_ids:
            cluster_ids.append(cid)
    cluster_ids = sorted(set(cluster_ids))

    fieldnames = [
        "cluster_id",
        "n_isolates",
        "alignment_total_chars",
        "alignment_length_per_seq",
        "seq_count_in_alignment",
        "starting_tree_present",
        "gubbins_status",
        "gubbins_exit_code",
        "n_recombination_blocks",
        "n_filtered_polymorphic_sites",
        "iqtree_status",
        "treefile_size_bytes",
        "iqtree_log_size_bytes",
        "confidence_tier",
        "notes",
    ]

    with args.output.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()

        for cid in cluster_ids:
            diag_path = args.diagnostics_dir / f"{cid}.diagnostics.log"
            tree_path = args.trees_dir / f"{cid}.final.treefile"
            log_path = args.trees_dir / f"{cid}.final.iqtree"
            filt_path = args.filtered_dir / f"{cid}.filtered_polymorphic_sites.fasta"
            gff_path = args.recomb_dir / f"{cid}.recombination_predictions.gff"

            diag = parse_diagnostics(diag_path)
            seq_count = int(diag["seq_count_in_alignment"]) if diag["seq_count_in_alignment"].isdigit() else None

            n_filtered = count_fasta_records(filt_path)
            filtered_cols = alignment_n_columns(filt_path)
            n_blocks = count_gff_recombination_blocks(gff_path)

            iqtree_status = determine_iqtree_status(tree_path, log_path, seq_count)
            treefile_size = tree_path.stat().st_size if tree_path.is_file() else 0
            log_size = log_path.stat().st_size if log_path.is_file() else 0

            n_iso = cluster_sizes.get(cid, n_filtered or 0)
            tier = confidence_tier(n_iso, diag["gubbins_status"], iqtree_status, filtered_cols)

            writer.writerow({
                "cluster_id": cid,
                "n_isolates": n_iso,
                "alignment_total_chars": diag["alignment_total_chars"],
                "alignment_length_per_seq": diag["alignment_length_per_seq"],
                "seq_count_in_alignment": diag["seq_count_in_alignment"],
                "starting_tree_present": diag["starting_tree_present"],
                "gubbins_status": diag["gubbins_status"],
                "gubbins_exit_code": diag["gubbins_exit_code"],
                "n_recombination_blocks": n_blocks,
                "n_filtered_polymorphic_sites": filtered_cols,
                "iqtree_status": iqtree_status,
                "treefile_size_bytes": treefile_size,
                "iqtree_log_size_bytes": log_size,
                "confidence_tier": tier,
                "notes": ";".join(diag["notes"]),
            })

    print(f"Wrote {args.output} with {len(cluster_ids)} cluster rows", file=sys.stderr)


if __name__ == "__main__":
    main()
