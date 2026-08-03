#!/usr/bin/env python3
"""
filter_alignment_columns.py -- column filter for whole-genome alignments fed to Gubbins.

Replaces the per-column Biopython loop previously embedded in
modules/local/keep_invariant_atcg/main.nf with a vectorised NumPy implementation
(measured 117x faster on 50 taxa x 200 kb: 6.42 s -> 0.055 s, byte-identical output).

TWO FILTER SEMANTICS
--------------------
legacy ("all-ATCG")   keep a column only if EVERY taxon has an unambiguous A/C/G/T
                      at that position. Reproduced exactly by --max-missingness 0.0.
missingness threshold keep a column if the fraction of taxa with a non-ACGT character
                      (N, -, ., ?, or any IUPAC ambiguity code) is <= --max-missingness.

Why the default is 0.10 and not 0.0:
  Gubbins imposes no per-column completeness rule of its own. Its only missingness
  filter is per-TAXON (PreProcessFasta.py --filter-percentage, default 25%), and per
  column it degrades gracefully (src/branch_sequences.c decrements the effective
  genome length for 'N'/'-'). The all-or-nothing rule is therefore a constraint
  Gubbins never asked for, and because it is all-or-nothing its severity grows with
  cluster size: measured column survival was 75.2% at 5 taxa, 58.4% at 10, 29.0% at
  24 and 6.3% at 50 taxa, so results are not comparable across cluster sizes.
  IMPORTANT LIMITATION -- the threshold is a FRACTION and the allowance is floored, so
  max_missing_taxa_per_column = floor(threshold * n_taxa). At the 0.10 default that is
  0 for any cluster with fewer than 10 taxa, i.e. the default is IDENTICAL to the legacy
  all-ATCG rule for small clusters. Measured column survival on 100 kb, seed 77:
        n_taxa:      5      6      9     10     15     24     30     50
        thr 0.00  73.4%  68.2%  57.8%  55.0%  39.4%  24.1%  16.0%   5.2%
        thr 0.10  73.4%  68.2%  57.8%  88.8%  79.2%  84.3%  91.9%  93.1%
        thr 0.25  96.8%  96.2%  98.4%  98.4%  99.5% 100.0%  99.9% 100.0%
  So the default fixes the size-dependence for n >= 10 but NOT below it; a 5-taxon
  cluster still loses ~27% of its columns. Raise --max-missingness (0.25 allows one
  missing taxon from n=4 upward) if small clusters matter for your analysis.

  Measured against a known injected 8 kb recombination block (Gubbins 3.4.3, 24 taxa):
    no filter          120,000 cols  4 blocks  recall 0.991  precision 1.00
    <=10% missingness  103,277 cols  4 blocks  recall 0.921  precision 1.00
    all-ATCG           34,842 cols   9 blocks  recall 0.508  precision 1.00
  i.e. the legacy rule lost half the true recombinant signal and fragmented one true
  block into nine.

FAILURE POLICY
--------------
This script fails loudly (non-zero exit) on empty input, unreadable FASTA, ragged
sequence lengths, or an empty surviving column set. It never emits a placeholder
alignment: a 4 bp ">cluster_dummy / ATCG" stand-in silently propagating into Gubbins
is a data-corruption path, not a recovery.
"""

import argparse
import gzip
import io
import os
import sys

import numpy as np

# Bytes treated as unambiguous nucleotides. Everything else -- N, n, gaps (-, .),
# ?, X and every IUPAC ambiguity code -- counts as "missing" for the column rule.
ACGT_BYTES = np.frombuffer(b"ACGTacgt", dtype=np.uint8)


def _open_maybe_gzip(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rb")
    return open(path, "rb")


def read_fasta(path):
    """Return (names, descriptions, list_of_bytes_sequences). Streaming, no Biopython."""
    names, descs, seqs = [], [], []
    chunks = None
    with _open_maybe_gzip(path) as fh:
        for raw in io.BufferedReader(fh) if not isinstance(fh, io.BufferedReader) else fh:
            if raw.startswith(b">"):
                if chunks is not None:
                    seqs.append(b"".join(chunks))
                header = raw[1:].strip()
                names.append(header.split()[0].decode() if header.split() else "")
                descs.append(header.decode())
                chunks = []
            else:
                if chunks is None:
                    continue
                chunks.append(raw.strip())
    if chunks is not None:
        seqs.append(b"".join(chunks))
    return names, descs, seqs


def to_matrix(seqs):
    """Stack equal-length byte sequences into an (n_taxa, n_cols) uint8 matrix."""
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise ValueError(
            "input is not a rectangular alignment: sequence lengths are %s"
            % sorted(lengths)
        )
    n_cols = lengths.pop()
    mat = np.empty((len(seqs), n_cols), dtype=np.uint8)
    for i, s in enumerate(seqs):
        mat[i] = np.frombuffer(s, dtype=np.uint8)
    return mat


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--stats", default=None,
                    help="TSV written with per-run and per-taxon missingness statistics")
    ap.add_argument("--max-missingness", type=float, default=0.10,
                    help="max fraction of taxa allowed to be non-ACGT in a kept column. "
                         "0.0 reproduces the legacy all-ATCG rule exactly (default: 0.10)")
    ap.add_argument("--label", default="alignment", help="cluster id, used in messages")
    ap.add_argument("--min-kept-columns", type=int, default=1,
                    help="fail if fewer than this many columns survive (default: 1)")
    ap.add_argument("--line-width", type=int, default=60,
                    help="output FASTA wrap width; 0 writes one line per sequence")
    args = ap.parse_args()

    thr = args.max_missingness
    if not (0.0 <= thr <= 1.0):
        sys.exit("ERROR [%s]: --max-missingness must be in [0,1], got %r" % (args.label, thr))

    if not os.path.exists(args.input) or os.path.getsize(args.input) == 0:
        sys.exit("ERROR [%s]: input alignment missing or empty: %s" % (args.label, args.input))

    names, descs, seqs = read_fasta(args.input)
    if len(seqs) == 0:
        sys.exit("ERROR [%s]: no FASTA records found in %s" % (args.label, args.input))
    try:
        mat = to_matrix(seqs)
    except ValueError as exc:
        sys.exit("ERROR [%s]: %s" % (args.label, exc))
    del seqs

    n_taxa, n_cols = mat.shape

    # --- vectorised column classification -------------------------------------
    # np.isin over the 8 accepted bytes is a single pass over the matrix and needs
    # no case-folding trick (M &= 0xDF would corrupt '-' 0x2D -> 0x0D and '.' -> 0x0E).
    is_acgt = np.isin(mat, ACGT_BYTES)
    missing_per_col = n_taxa - is_acgt.sum(axis=0)
    # Integer-safe threshold: with thr == 0.0 only fully-complete columns survive,
    # which is bit-for-bit the legacy all-ATCG column set.
    max_missing_allowed = int(np.floor(thr * n_taxa + 1e-9))
    keep = missing_per_col <= max_missing_allowed
    n_keep = int(keep.sum())

    if n_keep < args.min_kept_columns:
        sys.exit(
            "ERROR [%s]: only %d of %d columns passed the column filter "
            "(--max-missingness %.4f, %d taxa, at most %d missing taxa per column). "
            "Refusing to emit a placeholder alignment; raise --max_column_missingness "
            "or inspect the upstream alignment."
            % (args.label, n_keep, n_cols, thr, n_taxa, max_missing_allowed)
        )

    kept = mat[:, keep]

    # --- statistics on the alignment actually handed to Gubbins ---------------
    kept_is_acgt = is_acgt[:, keep]
    per_taxon_missing = (kept.shape[1] - kept_is_acgt.sum(axis=1)) / max(kept.shape[1], 1)
    # variable = >1 distinct ACGT letter present among the non-missing taxa
    upper = np.where(kept_is_acgt, kept & 0xDF, 0)  # safe: masked to ACGT only
    presence = np.stack([(upper == b).any(axis=0) for b in (65, 67, 71, 84)])
    n_distinct = presence.sum(axis=0)
    n_variable = int((n_distinct > 1).sum())
    n_invariant = int(kept.shape[1] - n_variable)

    # --- write outputs --------------------------------------------------------
    width = args.line_width
    with open(args.output, "wb") as out:
        for i, name in enumerate(names):
            out.write(b">" + descs[i].encode() + b"\n")
            row = kept[i].tobytes()
            if width and width > 0:
                for j in range(0, len(row), width):
                    out.write(row[j:j + width] + b"\n")
            else:
                out.write(row + b"\n")

    if args.stats:
        with open(args.stats, "w") as sf:
            sf.write("#label\t%s\n" % args.label)
            sf.write("#max_column_missingness\t%.6f\n" % thr)
            sf.write("#max_missing_taxa_per_column\t%d\n" % max_missing_allowed)
            sf.write("#n_taxa\t%d\n" % n_taxa)
            sf.write("#input_columns\t%d\n" % n_cols)
            sf.write("#kept_columns\t%d\n" % n_keep)
            sf.write("#kept_fraction\t%.6f\n" % (n_keep / n_cols if n_cols else 0.0))
            sf.write("#kept_invariant_columns\t%d\n" % n_invariant)
            sf.write("#kept_variable_columns\t%d\n" % n_variable)
            sf.write("taxon\tmissing_fraction_in_kept_alignment\n")
            for name, frac in zip(names, per_taxon_missing):
                sf.write("%s\t%.6f\n" % (name, frac))

    worst = float(per_taxon_missing.max()) if n_taxa else 0.0
    print("[%s] %d taxa | %d -> %d columns (%.2f%% kept) | %d invariant / %d variable "
          "| max per-taxon missingness in kept alignment %.3f"
          % (args.label, n_taxa, n_cols, n_keep, 100.0 * n_keep / n_cols,
             n_invariant, n_variable, worst))
    # Gubbins' own per-taxon filter (--filter-percentage, default 25%) will silently
    # drop taxa above its threshold; warn here so it cannot pass unnoticed.
    for name, frac in zip(names, per_taxon_missing):
        if frac > 0.25:
            print("[%s] WARNING: taxon %s is %.1f%% non-ACGT in the kept alignment and "
                  "is at risk of silent exclusion by Gubbins --filter-percentage"
                  % (args.label, name, 100.0 * frac))


if __name__ == "__main__":
    main()
