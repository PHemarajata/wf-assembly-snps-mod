# MIGRATION_NOTES.md — which defaults changed scientific output, and how to revert

Every behaviour change on this branch is gated by a parameter whose legacy value
restores the prior output. This file maps each one.

If you need to reproduce results from before this branch, use the
[legacy invocation](#reproducing-pre-branch-output) at the bottom.

---

## Parameters that alter scientific output

| parameter | new default | legacy value | what changes |
|---|---|---|---|
| `--max_column_missingness` | `0.10` | `0.0` | Column filter on the alignment fed to Gubbins. `0.0` is the old all-or-nothing A/T/C/G rule and reproduces prior output byte-identically. |
| `--cluster_split_method` | `similarity` | `order` | How components larger than `max_cluster_size` are split. `order` is the legacy contiguous-slice behaviour, arbitrary with respect to phylogeny. |
| `--iqtree_support` | `false` | `true` | Branch support on the final per-cluster tree. Previously always on and unconfigurable. **Turn it on for anything publishable** — support values are how a reader judges a clade. |
| `--iqtree_starting_model` | `GTR+G` | `MFP` | Model for the throwaway starting tree Gubbins discards. The old `MFP` ran ModelFinder over up to 968 models and then selected plain `JC`. |
| `--drop_reference_taxon` | `null` (drop when it is a duplicate) | `false` | Removes snippy-core's `Reference` taxon when it duplicates the cluster medoid. `false` restores the duplicated genome and the repeated `Reference` tips. |
| `--gubbins_filter_percentage` | `25` | `25` (but previously **unset**) | Same numeric value, but it was never passed before, so Gubbins' invisible default applied with no record. Setting it explicitly makes exclusions attributable. `100` disables taxon dropping entirely. |
| `--mash_threshold` | `0.03` (or `auto`) | any fixed number | `auto` derives the threshold from the collection. Pass an explicit number to pin it. |

### Not a parameter, but a behaviour change

| change | effect | how to detect the old behaviour |
|---|---|---|
| Gubbins now **fails** instead of writing empty placeholders | A crashed Gubbins stops the run rather than producing a tree with no recombination masking | old runs: `find <results> -name "*.recombination_predictions.gff" -empty` |
| IQ-TREE no longer writes star trees on `+ASC` failure | Degenerate `(a,b,c);` trees are gone except where <3 taxa genuinely leaves no topology | old runs: count internal nodes, see README |
| `INFILE_HANDLING_UNIX` makes contig IDs unique | Files whose IDs are already unique are reproduced **byte-for-byte**; only colliding IDs are rewritten to `<id>_<n>` | affects only inputs that would otherwise crash snippy |

---

## Parameters you should set explicitly

These are not migration issues, but they are easy to get wrong:

| parameter | why |
|---|---|
| `--gubbins_tree_builder rapidnj` | `params.config` defaults to `iqtree`. Only `low_spec` sets `rapidnj`, so on any other profile you must pass it. **`veryfasttree` is not a value Gubbins accepts and will fail every cluster.** |
| `--publish_dir_mode link` | Default is `copy`, which duplicates every published file on disk. `link` hard-links instead. Requires `outdir` and `work/` on one filesystem. |
| `--max_cluster_size 50` | `params.config` defaults to `100`; 50 is the tested production value. |
| `--merge_singletons false` | Default is already `false`. Do **not** turn it on: it pools all leftover genomes into a single bin with no similarity criterion. That bin is not a clade, and Gubbins' recombination calls inside it are not interpretable. |

---

## Reproducing pre-branch output

```bash
nextflow run . -profile low_spec,docker \
    --recombination_aware_mode true \
    --input  /path/to/assemblies \
    --outdir results_legacy \
    --mash_threshold 0.002 \
    --cluster_split_method order \
    --max_column_missingness 0.0 \
    --iqtree_support true \
    --drop_reference_taxon false \
    --starting_tree_builder iqtree \
    --iqtree_starting_model MFP
```

Cluster assignments and filtered alignments should match a pre-branch run.
**Trees will still differ where the old code emitted star trees** — it could not
produce a real tree there, so there is nothing to reproduce.

They will also differ if the old run hit a silent Gubbins failure, because the
old code produced no recombination masking in that case and this branch refuses
to continue.
