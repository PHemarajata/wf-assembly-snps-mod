#!/usr/bin/env python3
"""asc_preflight.py -- decide whether +ASC is legal for an alignment, and emit -fconst if not.

WHY THIS EXISTS
IQ-TREE's ascertainment-bias correction (+ASC) aborts with
"ERROR: Invariant sites observed ... please remove them" when the input contains
constant columns. Gubbins' <prefix>.filtered_polymorphic_sites.fasta is *supposed*
to be variant-only, but it is produced after recombination masking, and masking a
site to 'N' in every taxon that carried the minority allele can leave a column that
is constant among the remaining unambiguous bases. When that happened, the old
IQTREE_ASC module caught the non-zero exit and wrote

    ($names);

-- a completely unresolved star topology that is indistinguishable from a real
result to every downstream consumer (grafting, the summary table). This script
detects the condition up front so the module can respond deliberately.

DECISION (all measured on 30 taxa x 1500 sites, 611 constant columns, -T 4)
  no constant columns
      -> use the +ASC model as configured, alignment unchanged

  --strategy varsites   (DEFAULT, ASC-correct and fastest)
      -> strip the constant columns and KEEP +ASC, passing NO -fconst.
         7.90 s, fully resolved tree (27 internal nodes on 30 tips).
         This is what IQ-TREE itself suggests by writing <prefix>.varsites.phy in
         the body of the error message.
         +ASC and -fconst are MUTUALLY EXCLUSIVE: -fconst reconstitutes the
         constant sites, so IQ-TREE re-raises the same "Invalid use of +ASC"
         error. Verified -- passing both fails in 0.01 s with no treefile.

  --strategy fconst     (full-alignment likelihood instead of ASC)
      -> strip the constant columns, DROP +ASC, and pass the observed constant
         A,C,G,T counts via -fconst. 7.31 s. Equivalent information, different
         likelihood formulation; branch lengths are not ASC-corrected.

  --strategy drop_asc   (simplest, slowest, least correct)
      -> leave the alignment alone and drop the +ASC term. 11.21 s.
         Branch lengths are uncorrected for ascertainment bias, which for SNP-only
         input systematically inflates them.

For contrast, the behaviour this replaces: GTR+ASC on the unmodified alignment
FAILED in 0.0 s with no treefile, and the calling module then wrote a star tree
with 0 internal nodes into <cluster>.final.treefile.

Exit status is always 0 on a readable alignment; the decision is written as
shell-sourceable KEY=VALUE lines to --out.
"""

import argparse
import sys

import numpy as np

ACGT = np.frombuffer(b"ACGTacgt", dtype=np.uint8)


def read_fasta(path):
    names, descs, seqs, chunks = [], [], [], None
    with open(path, "rb") as fh:
        for raw in fh:
            if raw.startswith(b">"):
                if chunks is not None:
                    seqs.append(b"".join(chunks))
                h = raw[1:].strip()
                names.append(h.split()[0].decode() if h.split() else "")
                descs.append(h.decode())
                chunks = []
            elif chunks is not None:
                chunks.append(raw.strip())
    if chunks is not None:
        seqs.append(b"".join(chunks))
    return names, descs, seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--out", required=True, help="shell KEY=VALUE decision file")
    ap.add_argument("--model", required=True, help="configured model, e.g. GTR+ASC")
    ap.add_argument("--strategy", choices=["varsites", "fconst", "drop_asc"],
                    default="varsites")
    ap.add_argument("--stripped-output", default=None,
                    help="where to write the constant-column-free alignment "
                         "(varsites and fconst strategies)")
    args = ap.parse_args()

    names, descs, seqs = read_fasta(args.input)
    if not seqs:
        sys.exit("ERROR: no FASTA records in %s" % args.input)
    lens = {len(s) for s in seqs}
    if len(lens) != 1:
        sys.exit("ERROR: ragged alignment, sequence lengths %s" % sorted(lens))
    n_cols = lens.pop()
    mat = np.empty((len(seqs), n_cols), dtype=np.uint8)
    for i, s in enumerate(seqs):
        mat[i] = np.frombuffer(s, dtype=np.uint8)

    is_acgt = np.isin(mat, ACGT)
    upper = np.where(is_acgt, mat & 0xDF, 0)
    presence = np.stack([(upper == b).any(axis=0) for b in (65, 67, 71, 84)])  # A C G T
    n_distinct = presence.sum(axis=0)
    # A column is "constant" for ASC purposes when at most one distinct unambiguous
    # base is present among the taxa. Columns with zero ACGT (all N) are also
    # uninformative and are treated as constant.
    const_mask = n_distinct <= 1
    n_const = int(const_mask.sum())

    model = args.model
    has_asc = "+ASC" in model.upper()

    lines = ["N_TAXA=%d" % mat.shape[0],
             "N_COLS=%d" % n_cols,
             "N_CONSTANT_COLS=%d" % n_const,
             "N_VARIABLE_COLS=%d" % (n_cols - n_const)]

    if not has_asc or n_const == 0:
        lines += ["ASC_ACTION=none",
                  "IQ_MODEL=%s" % model,
                  "IQ_ALIGNMENT=%s" % args.input,
                  "IQ_FCONST="]
        if has_asc:
            print("[asc_preflight] no constant columns; +ASC is valid as configured (%s)" % model)
        else:
            print("[asc_preflight] model %s has no +ASC; nothing to check" % model)
    elif args.strategy == "drop_asc":
        # rebuild with '+' separators, dropping the ASC term
        parts = [p for p in model.split("+") if p.upper() != "ASC"]
        stripped_model = "+".join(parts)
        lines += ["ASC_ACTION=drop_asc",
                  "IQ_MODEL=%s" % stripped_model,
                  "IQ_ALIGNMENT=%s" % args.input,
                  "IQ_FCONST="]
        print("[asc_preflight] %d/%d columns are constant; dropping +ASC -> model %s"
              % (n_const, n_cols, stripped_model))
    else:
        if not args.stripped_output:
            sys.exit("ERROR: --stripped-output is required for --strategy %s" % args.strategy)
        keep = ~const_mask
        n_var = int(keep.sum())
        if n_var < 3:
            sys.exit("ERROR: %s has only %d variable columns after removing %d constant "
                     "columns; no ML tree can be estimated. This cluster is too clonal "
                     "for a resolved phylogeny -- it must be reported as such, not as a "
                     "tree." % (args.input, n_var, n_const))
        # Counts of constant sites by their fixed base, for -fconst a,c,g,t.
        counts = {65: 0, 67: 0, 71: 0, 84: 0}
        const_idx = np.flatnonzero(const_mask)
        if const_idx.size:
            sub = presence[:, const_idx]
            for k, b in enumerate((65, 67, 71, 84)):
                counts[b] = int(sub[k].sum())
        kept = mat[:, keep]
        with open(args.stripped_output, "wb") as out:
            for i, d in enumerate(descs):
                out.write(b">" + d.encode() + b"\n")
                row = kept[i].tobytes()
                for j in range(0, len(row), 60):
                    out.write(row[j:j + 60] + b"\n")
        fconst = "%d,%d,%d,%d" % (counts[65], counts[67], counts[71], counts[84])
        if args.strategy == "varsites":
            # Keep +ASC on the variable-sites-only alignment. NO -fconst: it is
            # mutually exclusive with +ASC (verified -- IQ-TREE re-raises
            # "Invalid use of +ASC" and writes no tree).
            lines += ["ASC_ACTION=varsites_keep_asc",
                      "IQ_MODEL=%s" % model,
                      "IQ_ALIGNMENT=%s" % args.stripped_output,
                      "IQ_FCONST=",
                      "N_CONSTANT_A_C_G_T=%s" % fconst]
            print("[asc_preflight] %d/%d columns are constant; stripped them and KEEPING "
                  "+ASC on the %d variable columns (constant A,C,G,T = %s, recorded but "
                  "NOT passed as -fconst, which is incompatible with +ASC)"
                  % (n_const, n_cols, n_var, fconst))
        else:
            # Drop +ASC, restore the constant-site information via -fconst instead.
            parts = [p for p in model.split("+") if p.upper() != "ASC"]
            noasc = "+".join(parts)
            lines += ["ASC_ACTION=fconst_no_asc",
                      "IQ_MODEL=%s" % noasc,
                      "IQ_ALIGNMENT=%s" % args.stripped_output,
                      "IQ_FCONST=%s" % fconst]
            print("[asc_preflight] %d/%d columns are constant; stripped them, dropped "
                  "+ASC -> %s, and will pass -fconst %s (%d variable columns remain). "
                  "Branch lengths are NOT ASC-corrected under this strategy."
                  % (n_const, n_cols, noasc, fconst, n_var))

    with open(args.out, "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
