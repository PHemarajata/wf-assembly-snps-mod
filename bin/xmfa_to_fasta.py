#!/usr/bin/env python3
"""Concatenate a parsnp XMFA core alignment into a single-block FASTA.

    xmfa_to_fasta.py backbone_alignment.xmfa backbone_core.fasta

BUILD_BACKBONE_TREE publishes `backbone_alignment.xmfa` -- parsnp's multi-block
format, one record per (block, genome). Most tools, Gubbins included, want a
single-block FASTA alignment with one record per genome. This concatenates the
blocks in order.

Verified on the 2,802-genome run: 4,168 blocks x 76 genomes -> 76 taxa x
5,993,126 columns, composition ACGTN- only.

Every LCB in this file carries all 76 sequences, so blocks concatenate directly
with no gap-filling. Records are keyed by the XMFA sequence INDEX (the integer
before the colon in '>1:594-2830 + cluster1 s1:p594'), and index -> name comes
from the '##SequenceIndex' / '##SequenceHeader' header comments.
"""
import sys, collections

xmfa, out = sys.argv[1], sys.argv[2]

names, idx = {}, None
parts = collections.defaultdict(list)
cur_idx, cur_seq = None, []
block_lens, nblocks = [], 0

def flush():
    global cur_idx, cur_seq
    if cur_idx is not None:
        parts[cur_idx].append("".join(cur_seq))
    cur_idx, cur_seq = None, []

with open(xmfa) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith("##SequenceIndex"):
            idx = int(line.split()[1])
        elif line.startswith("##SequenceHeader") and idx is not None:
            names[idx] = line.split(">", 1)[1].strip()
        elif line.startswith("#"):
            continue
        elif line.startswith("="):
            flush(); nblocks += 1
        elif line.startswith(">"):
            flush()
            cur_idx = int(line[1:].split(":", 1)[0])
        elif cur_idx is not None:
            cur_seq.append(line.strip())
flush()

lens = {i: sum(len(s) for s in v) for i, v in parts.items()}
L = set(lens.values())
print(f"  sequences: {len(parts)}   blocks: {nblocks}")
print(f"  concatenated length(s): {sorted(L)[:3]}{' ...' if len(L)>3 else ''}")
if len(L) != 1:
    print("  ERROR: ragged concatenation — blocks do not carry equal-length records")
    sys.exit(1)

with open(out, "w") as fh:
    for i in sorted(parts):
        fh.write(f">{names.get(i, f'seq{i}')}\n")
        s = "".join(parts[i])
        for j in range(0, len(s), 60):
            fh.write(s[j:j+60] + "\n")
print(f"  wrote {out}: {len(parts)} taxa x {sorted(L)[0]:,} columns")
