\
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Standalone tree grafting CLI (Bio.Phylo-based)
----------------------------------------------
Algorithm: leaf-expansion grafting
- Remove the representative tip from each *cluster* tree.
- Attach the remaining cluster subtree as a sister to the representative tip
  at the matching tip in the *backbone* tree.

This avoids Nextflow/gotree coupling and performs structural tree surgery
purely in Python using Biopython's Bio.Phylo.

Requirements
------------
pip install biopython

Usage (examples)
----------------
python graft_trees.py \
  --backbone backbone.treefile \
  --clusters 'Clusters/**/cluster_*.final.treefile' \
  --reps cluster_representatives.tsv \
  --out-tree Final_Results/global_grafted.treefile \
  --report   Final_Results/grafting_report.txt \
  --log      Final_Results/grafting_log.txt

Author: ChatGPT (generated)
License: MIT
"""
import argparse
import glob
import os
import sys
from typing import Dict, List, Optional, Set

# We depend on Biopython for this implementation.
try:
    from Bio import Phylo
    from Bio.Phylo.BaseTree import Tree, Clade
except Exception as e:
    sys.stderr.write(
        "ERROR: This script requires Biopython.\\n"
        "Install with: pip install biopython\\n"
    )
    raise

def normalize_label(s: str) -> str:
    """Normalize tip labels to improve matching:
    - strip surrounding quotes
    - remove common extensions
    - replace spaces with underscores
    """
    if s is None:
        return ""
    s = s.strip()
    if len(s) >= 2 and ((s[0] == '"' and s[-1] == '"') or (s[0] == "'" and s[-1] == "'")) :
        s = s[1:-1]
    for ext in (".fa", ".fasta", ".fa.gz", ".fna", ".treefile"):
        if s.endswith(ext):
            s = s[: -len(ext)]
    s = s.replace(" ", "_")
    return s

def load_tree(path: str) -> 'Tree':
    """Load a Newick tree using Bio.Phylo; preserve underscores."""
    with open(path, "rt") as fh:
        tree = Phylo.read(fh, "newick")
    return tree

def all_tip_labels(tree: 'Tree') -> List[str]:
    return [t.name for t in tree.get_terminals() if t.name is not None]

def build_parent_map(tree: 'Tree') -> Dict['Clade', Optional['Clade']]:
    """Return a mapping: child -> parent for all nodes (root's parent is None)."""
    parent = {tree.root: None}
    stack = [tree.root]
    while stack:
        node = stack.pop()
        for ch in getattr(node, "clades", []):
            parent[ch] = node
            stack.append(ch)
    return parent

def find_tip_by_label(tree: 'Tree', label: str) -> Optional['Clade']:
    """Find first terminal with exact label (case-sensitive)."""
    for t in tree.get_terminals():
        if t.name == label:
            return t
    return None

def infer_representative(cluster_tree: 'Tree',
                         backbone_map: Dict[str, str]) -> Optional[str]:
    """Pick a representative leaf from cluster that also exists (normalized) in backbone.
    If none, return the first terminal label.
    """
    cl_tips = all_tip_labels(cluster_tree)
    if not cl_tips:
        return None
    for lab in cl_tips:
        n = normalize_label(lab)
        if n in backbone_map:
            return lab
    # fallback: first tip
    return cl_tips[0]

def rename_conflicts(cluster_tree: 'Tree',
                     forbidden: Set[str],
                     rep_label: str,
                     suffix: str = "_dup") -> Dict[str, str]:
    """Rename any cluster tip labels (except rep_label) that collide with forbidden set.
    Returns a dict of {old -> new} names applied.
    """
    applied = {}
    for t in cluster_tree.get_terminals():
        if t.name is None:
            continue
        if t.name == rep_label:
            continue
        if t.name in forbidden:
            base = t.name
            k = 1
            new = f"{base}{suffix}{k}"
            while new in forbidden or find_tip_by_label(cluster_tree, new) is not None:
                k += 1
                new = f"{base}{suffix}{k}"
            applied[t.name] = new
            t.name = new
    return applied

def prune_representative(cluster_tree: 'Tree', rep_label: str) -> bool:
    """Remove the representative tip from the cluster tree; collapse unary nodes.
    Returns True if pruning succeeded or the tip wasn't present (treated as already pruned).
    """
    try:
        rep_tip = find_tip_by_label(cluster_tree, rep_label)
        if rep_tip is None:
            return True
        cluster_tree.prune(target=rep_tip)
        return True
    except Exception as e:
        return False

def choose_backbone_label(backbone_map: Dict[str, str],
                          rep_label: str) -> Optional[str]:
    """Return the *exact* backbone tip label to graft at, given a representative tip label
    from the cluster tree. Matching tries:
      1) normalized exact match
      2) fuzzy contains match (either direction)
    """
    rep_norm = normalize_label(rep_label)
    if rep_norm in backbone_map:
        return backbone_map[rep_norm]
    # fuzzy
    for k, v in backbone_map.items():
        if rep_norm and (k.find(rep_norm) != -1 or rep_norm.find(k) != -1):
            return v
    return None

def graft_cluster_into_backbone(
    backbone: 'Tree',
    cluster_subtree: 'Tree',
    at_tip_label: str,
    parent_edge_mode: str = "keep"
) -> bool:
    """Attach cluster_subtree.root as sister to the 'at_tip_label' leaf in backbone.

    parent_edge_mode:
      - "keep": preserve the original branch length of the representative leaf on the
                new internal node (join), set child (rep leaf) branch length to 0.0
      - "zero": set the new join branch length to 0.0 and keep the old rep leaf length
                (default for leaf stays untouched)
    """
    # Locate target leaf and its parent
    rep_leaf = find_tip_by_label(backbone, at_tip_label)
    if rep_leaf is None:
        return False

    parent_map = build_parent_map(backbone)
    parent = parent_map.get(rep_leaf)
    if parent is None:
        # leaf is at the root: make a new root
        new_root = Clade(name=None, branch_length=0.0)
        # old root becomes child
        new_root.clades.append(backbone.root)
        # cluster subtree root becomes sibling
        if cluster_subtree.root.branch_length is None:
            cluster_subtree.root.branch_length = 0.0
        new_root.clades.append(cluster_subtree.root)
        backbone.root = new_root
        return True

    # Create a new internal "join" clade
    old_len = getattr(rep_leaf, "branch_length", None)
    join = Clade(name=None)

    # Replace rep_leaf with join under parent
    try:
        idx = parent.clades.index(rep_leaf)
    except ValueError:
        return False
    parent.clades[idx] = join

    # Configure branch lengths per mode
    if parent_edge_mode == "keep":
        join.branch_length = old_len  # edge from parent -> join
        rep_leaf.branch_length = 0.0  # edge from join -> rep
    else:  # "zero"
        join.branch_length = 0.0
        # keep rep_leaf.branch_length as-is

    # Ensure cluster root has a defined length
    if cluster_subtree.root.branch_length is None:
        cluster_subtree.root.branch_length = 0.0

    # Attach children to join
    join.clades = [rep_leaf, cluster_subtree.root]
    return True

def read_reps_map(tsv: Optional[str]) -> Dict[str, str]:
    """Read cluster->representative mapping from TSV: cluster_id<tab>rep_label"""
    mapping = {}
    if not tsv:
        return mapping
    if not os.path.exists(tsv) or os.path.getsize(tsv) == 0:
        return mapping
    with open(tsv, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            cid, rep = parts[0], parts[1]
            mapping[cid] = rep
            mapping[f"{cid}.final"] = rep
    return mapping

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Graft cluster trees into a backbone using leaf-expansion algorithm (Bio.Phylo)."
    )
    ap.add_argument("--backbone", required=True, help="Backbone Newick treefile")
    ap.add_argument("--clusters", required=True, action="append",
                    help="Glob for cluster treefiles (repeatable). Example: 'Clusters/**/cluster_*.final.treefile'")
    ap.add_argument("--reps", default=None, help="TSV: cluster_id<TAB>representative_label (optional)")
    ap.add_argument("--out-tree", required=True, help="Output combined tree (Newick)")
    ap.add_argument("--report", required=True, help="Text report output path")
    ap.add_argument("--log", required=True, help="Detailed log output path")
    ap.add_argument("--parent-edge-mode", choices=["keep", "zero"], default="keep",
                    help="How to set branch lengths at the join node [default: keep].")
    ap.add_argument("--rename-conflicts", action="store_true",
                    help="Rename cluster labels that collide with backbone labels (excluding representative)")
    ap.add_argument("--dry-run", action="store_true", help="Plan only; do not write out-tree")
    args = ap.parse_args(argv)

    # Expand cluster globs
    cluster_files: List[str] = []
    for pat in args.clusters:
        cluster_files.extend(glob.glob(pat, recursive=True))
    # Deduplicate keep order
    seen = set()
    cluster_files = [f for f in cluster_files if not (f in seen or seen.add(f))]

    os.makedirs(os.path.dirname(args.out_tree) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(args.report) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(args.log) or ".", exist_ok=True)

    with open(args.log, "wt") as log, open(args.report, "wt") as rep:
        log.write("Starting tree grafting (leaf-expansion, Bio.Phylo)\n")
        log.write(f"Backbone: {args.backbone}\n")
        log.write(f"Clusters (expanded): {len(cluster_files)}\n")
        for f in cluster_files[:20]:
            log.write(f"  - {f}\n")
        if len(cluster_files) > 20:
            log.write(f"  ... and {len(cluster_files)-20} more\n")

        # Load backbone
        try:
            backbone = load_tree(args.backbone)
        except Exception as e:
            log.write(f"ERROR: Failed to read backbone: {e}\n")
            return 2

        bb_labels = all_tip_labels(backbone)
        log.write(f"Backbone tips: {len(bb_labels)}\n")
        backbone_map = {}
        for lab in bb_labels:
            backbone_map.setdefault(normalize_label(lab), lab)

        reps_map = read_reps_map(args.reps)
        if reps_map:
            log.write(f"Loaded representatives from {args.reps}: {len(reps_map)} entries\n")
        else:
            log.write("No representatives TSV provided or empty; will infer per cluster.\n")

        grafted_ok = 0
        grafted_fail = 0

        for cf in cluster_files:
            if not os.path.exists(cf) or os.path.getsize(cf) == 0:
                log.write(f"SKIP: Missing/empty cluster file: {cf}\n")
                grafted_fail += 1
                continue

            cluster_id = os.path.splitext(os.path.basename(cf))[0]
            log.write(f"\nProcessing cluster: {cluster_id}\n")

            try:
                ctree = load_tree(cf)
            except Exception as e:
                log.write(f"  ERROR: Failed to read cluster tree: {e}\n")
                grafted_fail += 1
                continue

            cl_labels = all_tip_labels(ctree)
            log.write(f"  Cluster tips: {len(cl_labels)}\n")

            rep_label = reps_map.get(cluster_id)
            if not rep_label:
                rep_label = infer_representative(ctree, backbone_map)
                if rep_label:
                    log.write(f"  Rep inferred: {rep_label}\n")
            else:
                log.write(f"  Rep from map: {rep_label}\n")

            if not rep_label:
                log.write("  ERROR: Could not determine representative.\n")
                grafted_fail += 1
                continue

            bb_tip_label = choose_backbone_label(backbone_map, rep_label)
            if not bb_tip_label:
                log.write(f"  ERROR: No matching tip in backbone for rep '{rep_label}'.\n")
                grafted_fail += 1
                continue
            log.write(f"  Backbone tip to graft at: {bb_tip_label}\n")

            rep_in_cluster = None
            if rep_label in cl_labels:
                rep_in_cluster = rep_label
            else:
                rn = normalize_label(rep_label)
                norm_map = {normalize_label(x): x for x in cl_labels}
                rep_in_cluster = norm_map.get(rn)
                if not rep_in_cluster:
                    for k, v in norm_map.items():
                        if rn and (k.find(rn) != -1 or rn.find(k) != -1):
                            rep_in_cluster = v
                            break
            if not rep_in_cluster:
                log.write(f"  ERROR: Representative '{rep_label}' not found in cluster tree.\n")
                grafted_fail += 1
                continue
            log.write(f"  Cluster representative label resolved: {rep_in_cluster}\n")

            if args.rename_conflicts:
                applied = rename_conflicts(ctree, set(bb_labels), rep_in_cluster, suffix="_dup")
                if applied:
                    log.write(f"  Renamed {len(applied)} conflicting labels in cluster:\n")
                    for old, new in list(applied.items())[:10]:
                        log.write(f"    {old} -> {new}\n")
                    if len(applied) > 10:
                        log.write(f"    ... and {len(applied)-10} more\n")

            if not prune_representative(ctree, rep_in_cluster):
                log.write("  ERROR: Failed to prune representative from cluster tree.\n")
                grafted_fail += 1
                continue

            if len(all_tip_labels(ctree)) == 0:
                log.write("  WARNING: Cluster subtree empty after pruning; skipping.\n")
                grafted_fail += 1
                continue

            ok = graft_cluster_into_backbone(
                backbone=backbone,
                cluster_subtree=ctree,
                at_tip_label=bb_tip_label,
                parent_edge_mode=args.parent_edge_mode,
            )
            if not ok:
                log.write("  ERROR: Grafting operation failed.\n")
                grafted_fail += 1
                continue

            bb_labels = all_tip_labels(backbone)
            backbone_map = {normalize_label(l): l for l in bb_labels}
            log.write("  SUCCESS: grafted cluster subtree.\n")
            grafted_ok += 1

        if not args.dry_run:
            with open(args.out_tree, "wt") as fh:
                Phylo.write(backbone, fh, "newick")

        rep.write("TREE GRAFTING REPORT\n")
        rep.write("===================\n")
        rep.write(f"Backbone tree: {args.backbone}\n")
        rep.write(f"Total cluster trees: {len(cluster_files)}\n")
        rep.write(f"Successfully grafted: {grafted_ok}\n")
        rep.write(f"Failed grafting: {grafted_fail}\n")
        if not args.dry_run:
            rep.write(f"Final tree: {args.out_tree}\n")
            rep.write(f"Final tree leaves: {len(all_tip_labels(backbone))}\n")

        log.write("\nDone.\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
