# AUDIT_REPORT.md — defects, evidence, and measurements

Every entry here was found by **running the workflow and reading its output**,
not by code inspection. Each records what was measured, on what data.

Reference datasets:

- **112 genomes** — *B. pseudomallei* subset, 6.96–7.39 Mb, median 91 contigs.
- **2,795 genomes** — full de-duplicated collection, 18.9 GB, total pairwise Mash
  span **0.01088**.

---

## A. Silent success — the defect class that matters most

The original code repeatedly caught a failure, wrote an empty or placeholder
output, and exited 0. Nextflow reported success and downstream processes
consumed meaningless data.

### A.1 `GUBBINS_CLUSTER` wrote placeholders on failure

```groovy
echo "ERROR: Gubbins failed; creating placeholder outputs." >> "${DIAG}"
: > "${cluster_id}.filtered_polymorphic_sites.fasta"
: > "${cluster_id}.recombination_predictions.gff"
: > "${cluster_id}.node_labelled.final_tree.tre"
# ...no non-zero exit
```

Nextflow showed `GUBBINS_CLUSTER 2 of 2 ✔` for a Gubbins that had failed on both
clusters. Downstream cannot distinguish "no recombination detected" from
"Gubbins crashed" — both are an empty GFF.

**Fixed:** both the crash path and the exited-0-but-empty path now fail and dump
the diagnostics log to stderr. The two *pre-flight* guards (missing alignment,
<3 sequences) still exit 0 — those are policy skips with a recorded reason, not
crashes.

**Consequence of the fix:** this is what exposed A.2, C.2 and C.3 below. Under
the old behaviour all three would have produced plausible-looking output.

### A.2 IQ-TREE wrote star trees

`GTR+ASC` aborts on any invariant column, which happens routinely in clonal
alignments. The module caught the failure and wrote `($names);` — zero internal
nodes — into `*.final.treefile`, which grafting consumed as a success.

**Fixed** with a constant-site preflight (`bin/asc_preflight.py`) that keeps
`+ASC` on variable sites only. The star-tree fallback is gone.

---

## B. Gubbins never ran under `low_spec`

`conf/profiles/low_spec.config` set:

```groovy
gubbins_first_tree_builder = 'veryfasttree'
gubbins_tree_builder       = 'veryfasttree'
```

Gubbins 3.3.5 rejects it:

```
run_gubbins.py: error: argument --tree-builder/-t: invalid choice: 'veryfasttree'
(choose from 'raxml', 'raxmlng', 'iqtree', 'iqtree-fast', 'fasttree', 'hybrid', 'rapidnj')
```

Gubbins exited 2 on **every** cluster. Three things hid it:

1. A.1 converted the exit-2 into a green checkmark.
2. `gubbins_tree_builder` is listed in `schema_ignore_params`, so schema
   validation skips it — **even a real `nextflow run` would not have caught it**.
3. Only `run_gubbins.py` knows the valid set.

The comment justifying the setting claimed veryfasttree was *"~14% faster than
rapidnj as the first-tree builder (371.1 s vs 431.0 s at n=50)"*. That cannot
have been measured through Gubbins. Treat it as a standalone VeryFastTree timing
written into config as if it were a Gubbins setting.

**Fixed:** `rapidnj`, the fastest builder Gubbins accepts.

### B.1 Measured Gubbins cost with `rapidnj`

112-genome set, 8 clusters, 22 cores:

| taxa | wall | peak RSS |
|---|---|---|
| 27 | 2m 35s | 916 MB |
| 24 | 3m 06s | 428 MB |
| 21 | 2m 09s | 416 MB |
| 12 | 35 s | 650 MB |
| 3 | 19–25 s | 431–863 MB |

**Cost is not monotonic in taxon count** — 24 taxa cost more than 27 — so the
`0.110·n^1.97` pure-*n* fit does not describe this workload. SNP density matters
at least as much as *n*.

### B.2 Gubbins does not scale with threads

27-taxon alignment, rapidnj, identical input:

| threads | 1 | 2 | 4 | 8 |
|---|---|---|---|---|
| wall | 164 s | 144 s | 140 s | 144 s |

1→2 buys 12%; past 2 there is nothing. Allocating 10 cpus per Gubbins task
wastes cores that should run more clusters concurrently.

---

## C. Container and environment defects

### C.1 Containers ran as root

`docker.runOptions` was written as a bare key inside `profiles { docker { … } }`.
That block is a **profile named `docker`**, not the docker config **scope**, so
dotted keys landed and bare ones did not:

```
before:  docker { enabled, cacheDir }
after:   docker { enabled, runOptions, cacheDir }
```

Containers therefore ran as root, task outputs were root-owned, and
`publish_dir_mode='link'` failed:

```
java.nio.file.FileSystemException: ... Operation not permitted
```

**The cause is ownership, not filesystems.** Everything was on one ext4 volume;
`fs.protected_hardlinks=1` forbids hardlinking a file you do not own. Reproduced
directly with `ln` outside Nextflow.

Present in all five docker-enabling profiles, so DGX runs were also running
containers as root.

### C.2 `/dev/shm` is 64 MB in Docker

Gubbins' `pyjar` opens a `multiprocessing.Pool`, whose semaphores and shared
arrays live in `/dev/shm`:

```
OSError: [Errno 28] No space left on device
  File "multiprocessing/synchronize.py", line 57, in __init__
    sl = self._semlock = _multiprocessing.SemLock(
```

**This reads like a full disk and is not one** — 243 GB was free, and 2,794
alignments plus the backbone tree had already completed on that filesystem.

Hit on the 2,795-genome run at cluster_26: 46 taxa, **63 MB** snp-sites
alignment against the **64 MB** default, after 21 smaller clusters had passed.
34 of the 61 clusters carry 50 taxa, so nearly every remaining Gubbins task
would have failed the same way.

**Fixed** with `--shm-size=2g`, a cap on a tmpfs that reserves nothing.

### C.3 The IQ-TREE container has no python

`IQTREE_ASC` declared `conda "… conda-forge::python=3.10 …"` and
`container "quay.io/biocontainers/iqtree:2.2.6"`, then shelled out to
`asc_preflight.py`. The container has neither `python` nor `python3`:

```
env: can't execute 'python3': No such file or directory
```

So the module worked under `-profile conda` and could never work in a container.
The two execution modes declared different toolsets for the same process.

**Fixed** by moving the preflight into its own process (`ASC_PREFLIGHT`) using
the pyseer image already pinned elsewhere in the repo. IQ-TREE stays pinned at
2.2.6 — the obvious alternative, the Gubbins container, ships IQ-TREE 2.3.0 and
would silently change the version building the final publishable trees.

### C.4 Duplicate contig IDs

17 of 112 assemblies are headed:

```
>SAMPLE.fasta 1
>SAMPLE.fasta 2
```

Unique as *lines*, but a FASTA record's ID is the first whitespace token, so all
contigs share one ID. snippy refuses such a reference outright. Affected files
include 134- and 227-contig drafts where every contig collapses to a single ID.

**Fixed** in `INFILE_HANDLING_UNIX` via `bin/dedupe_contig_ids.awk`. Only
colliding IDs are rewritten; verified that all 95 conforming assemblies are
reproduced **byte-for-byte**.

---

## D. Naming and identity defects

### D.1 Representative filename collision

`SELECT_CLUSTER_REPRESENTATIVE` emitted the exact name `representative.fa`,
which collides the moment two clusters' outputs are staged together. The
collision was the symptom; the real defect was that `COLLECT_REPRESENTATIVES`
derived each sequence header from the **filename**:

```bash
rep_id=$(basename $rep_file .fa)          # now always "representative"
awk -v rep_id="$rep_id" '/^>/ { print ">" rep_id; next }'
```

`representatives.fa` feeds `BUILD_BACKBONE_TREE`, whose headers become backbone
**tip labels**, and `GRAFT_TREES` joins those against a map built from
`representative_id.txt`. **No tip could ever have matched its cluster.**
Nextflow's collision check is the only reason this surfaced as an error rather
than a quietly wrong backbone tree.

**Fixed** with cluster-scoped filenames and by reading the real label from
`representative_id.txt` — the same source grafting uses, so the two agree by
construction.

### D.2 snippy-core's duplicate `Reference` taxon

snippy-core emits the reference as an extra taxon named `Reference`. With the
default medoid reference that taxon is **by construction a duplicate of a sample
already in the cluster**:

```
cluster_1 samples : 2002721684, 2010007509, 2011756296
medoid            : 2010007509
core.full.aln     : 2002721684, 2010007509, 2011756296, Reference
                                └────── same genome ──────┘
```

Four consequences: Gubbins receives a zero-divergence pair; IQ-TREE builds an
artificial zero-length cherry and the duplicate shifts other branch lengths
(1701: 0.0073281644 → 0.0074756798 once removed); the taxon count feeding ASC
correction is inflated; and because grafting prunes the *representative* tip and
not `Reference`, the published global tree carried **one `Reference` tip per
cluster, each a different organism**.

**Fixed**, gated by `--drop_reference_taxon`. Kept when
`use_global_reference=true`, where the reference is a genuine extra taxon.

---

## E. Resource and profile defects

### E.1 Profiles targeted a pipeline that no longer exists

Auditing `withName:` selectors against processes actually on the
recombination-aware path:

| profile | processes with a matching selector |
|---|---|
| `dgx_station_a100_updated` | **4 of 21** |
| `local_workstation_rtx4070` | 12 of 21 |
| `low_spec` | 15 of 21 |

Dead selectors: `SNIPPY_ALIGN`, `MASH_SKETCH`, `MASH_DIST`, `SKA_BUILD`,
`GPU_ACCELERATED`. Most importantly **no block for `SNIPPY_SCATTER`**, which is
89% of total CPU, so it fell through to a generic label.

### E.2 Over-allocated memory caps concurrency

Nextflow's local executor packs tasks by **declared** memory, so an
over-allocation directly limits how many tasks run:

| process | declared | measured peak | over |
|---|---|---|---|
| `GUBBINS_CLUSTER` | 22 GB | 0.89 GB | 25× |
| `IQTREE_ASC` | 18 GB | 0.05 GB | 360× |
| `IQTREE_FAST` | 12 GB | 0.87 GB | 14× |

At 22 GB each, only two Gubbins fit in 58 GB — which is exactly why
`maxForks = 2` was set. **Right-sizing memory is what makes concurrency
possible.**

### E.3 `BUILD_BACKBONE_TREE` was under-allocated

parsnp over one representative per cluster. Measured:

| representatives | peak RSS | wall |
|---|---|---|
| 8 | 5.30 GB | — |
| 138 | **12.80 GB** | 929 s |

Scaling is sub-linear. `low_spec` allocated **12 GB**, below the measured peak.

Its memory arithmetic also excluded this process "because nothing else is
running", which is false: `BUILD_BACKBONE_TREE` consumes
`COLLECT_REPRESENTATIVES`, which needs only the representatives — ready early —
so the backbone **overlaps** the per-cluster work. Confirmed by observation
during the 2,795-genome run.

### E.4 Profile precedence follows definition order, not command-line order

```
-profile bp,low_spec   → max_cluster_size = 25
-profile low_spec,bp   → max_cluster_size = 25
```

Earlier documentation claimed the reverse order flips the winner. It does not,
because these profiles merge via `includeConfig` and `low_spec` is defined later
in `nextflow.config`. **Composing `low_spec` with a larger machine profile
silently clamps you to 16 cores / 30 GB.**

---

## F. `--alignment_method ska` undercalls recombination

SKA is 2.7× faster end to end (9m11s vs 24m53s on 112 genomes; 1.0 vs 3.0
CPU-hours) and **wrong for this workflow**.

A split k-mer matches only when both flanking half-k-mers match exactly, so a
second SNP inside the flank destroys the match. The loss is a function of **SNP
spacing**. Ratio of SKA variable sites to snippy variable sites on the same
cluster, same samples, same reference:

| gap to next SNP | cluster_12 (21 taxa) | cluster_9 (12 taxa) |
|---|---|---|
| 1–10 bp | **0.11×** | **0.05×** |
| 11–20 bp | 0.41× | 0.34× |
| 21–31 bp | 0.75× | 0.72× |
| 51–100 bp | 0.86× | 0.83× |
| 101–500 bp | 0.96× | 0.93× |
| 501–5000 bp | 1.09× | 1.06× |

Monotonic in spacing, at parity beyond ~100 bp, reproduced on two independent
clusters. That shape is the k-mer flank. **Ruled out as alternatives:** both
alignments are essentially complete (missing 0.163% snippy vs 0.132% SKA, zero
all-missing columns) and Gubbins retained every taxon in both runs.

Gubbins detects recombination as regions of **elevated SNP density**, so
deleting ~90% of the tightest SNP clusters deletes the signal. The error
amplifies: 23% fewer variable sites became **54% fewer recombination blocks**
(2,388 vs 5,148), and topologies then disagreed in 3 of the 4 clusters large
enough to have a non-trivial topology.

SKA also emits IUPAC ambiguity codes that Gubbins rejects outright
(7,027–58,370 per cluster); that is fixed, but it is why the path had never run
to completion before.

**Conclusion: use snippy.** SKA remains reasonable for recombination-free uses.

---

## G. Threshold selection

### G.1 A fixed threshold does not transfer

`0.003` was calibrated on 112 genomes. On the 2,795-genome collection it
silently drops **203 genomes** into components of <3 taxa, which cannot produce
a tree.

### G.2 Coherence measured across thresholds

Collection span (max pairwise Mash) = **0.01088**.

| threshold | clusters ≥3 | dropped | within-cluster max | % of span |
|---|---|---|---|---|
| 0.003 | 109 | 203 | 0.00698 | 64% |
| 0.005 | 67 | 58 | 0.00946 | 87% |
| **0.006425** | **61** | **1** | **0.00946** | **87%** |
| 0.006549 | 61 | 0 | 0.01061 | **98%** |

At 0.006549 a single cluster spans 98% of the entire collection. That is not a
clade, and Gubbins' assumption of limited diversity sharing a recent common
ancestor breaks — its recombination calls there are not interpretable.

### G.3 Acceptance is not monotonic

```
0.005802  dropped 6  within_max 0.010318  rejected
0.005927  dropped 5  within_max 0.009463  ok
0.006176  dropped 3  within_max 0.010458  rejected
0.006300  dropped 2  within_max 0.010458  rejected
0.006425  dropped 1  within_max 0.009463  ok        <- selected
0.006549  dropped 0  within_max 0.010458  rejected
```

Raising the threshold changes which component a borderline genome joins, and
average-linkage re-splitting can then *lower* a cluster's span. Accepted and
rejected states interleave.

**Consequence for implementation:** a bisection assumes a single boundary and
walks past these windows — an early implementation converged on 0.005162 and
dropped 38 genomes. `--mash_threshold auto` therefore **scans uniformly** and
picks the best over all candidates.

**Consequence for the user:** the selected threshold can sit in a narrow
accepted window, meaning the clustering is genuinely sensitive there. Read
`Summaries/threshold_sweep.tsv`.

### G.4 Cluster structure is sensitive to which genomes are in the collection

Removing 8 genomes from the 2,795-genome set (0.3% of the collection) changed
the auto-selected threshold and reorganised most clusters:

| | 2,795 genomes | 2,787 genomes (8 removed) |
|---|---|---|
| selected threshold | 0.006425 | **0.005194** |
| clusters ≥3 | 61 | 93 |
| genomes dropped | 1 | **34** |

Of 93 groupings in the smaller set, only 42 survive unchanged and only 25 keep
the same `cluster_id`.

**This is a real data effect, not an artefact of the selection algorithm.** The
collection span is identical in both cases (0.010879) and so is the coherence
limit (0.009791); the same grid region is evaluated. What differs is the
clustering itself — near 0.006 the 2,795-genome set yields a maximum
within-cluster span of 0.009463 (accepted) while the 2,787-genome set yields
0.010382 (rejected). Those 8 genomes act as BRIDGES: single-linkage components
chain through them, so their presence changes which genomes share a component
and therefore how far that component spans.

Two consequences:

1. **Incremental addition is not viable on this path.** Adding a handful of
   genomes can reorganise the clustering wholesale, so cached per-cluster work
   cannot be reused and a clean re-run is required. There is also no supported
   placement mode here — UShER exists only in `assembly_snps_scalable.nf`, not in
   the recombination-aware workflow.
2. **Cluster membership is a property of the sample set, not just of the
   population.** Two analyses of overlapping collections can partition the same
   genomes differently. Report the threshold and the collection alongside any
   cluster-based result, and do not treat cluster identity as stable across
   analyses.

### G.5 `merge_singletons` produces an incoherent bin

`merge_singletons=true` pools all leftover genomes into one cluster **with no
similarity criterion**. On the 112-genome set at threshold 0.002 that bin had
internal distances up to 0.006985 against a dataset maximum of 0.006995 — 4×
more divergent than a genuine cluster. It is not a clade. Leave it `false`.

---

## H. Performance measurements

### H.1 Where the time goes (112 genomes, 8 clusters)

| process | n | Σ CPU | share |
|---|---|---|---|
| `SNIPPY_SCATTER` | 96 | **7,333 s** | **89%** |
| `GUBBINS_CLUSTER` | 8 | 593 s | 7% |
| `BUILD_BACKBONE_TREE` | 1 | 78 s | 1% |

Within one snippy task: `bwa mem` + sort 63%, freebayes 24%, **`bwa index` only
6%**. Hoisting the reference index out of the per-sample loop — long assumed to
be the win — is worth ~5%.

### H.2 Snippy concurrency saturates early

Three full 112-genome runs, identical inputs, varying only `maxForks`:

| forks | mean/task | Σ CPU | alignment-phase wall | run wall | CPU-h |
|---|---|---|---|---|---|
| 8 | 76 s | 7,333 s | 917 s | 24m53s | 3.0 |
| 12 | 91 s | 8,718 s | **727 s** | 22m13s | 4.5 |
| 18 | 136 s | 13,027 s | 724 s | 21m24s | 5.8 |

Past 12 forks the alignment phase stops improving: 18 buys 3 s over 12 for 50%
more CPU. **12 is the knee.**

Cause is CPU topology, not memory (peak RSS 707 MB at 18 concurrent, swap
untouched). `nproc` reports 22 on a Core Ultra 9 185H, but that is **16 physical
cores** in a hybrid arrangement (6 P-cores with SMT, 8 E-cores, 2 LP-E). Past
~12 tasks work lands on E-cores that run `bwa mem` far slower. **Do not size
`maxForks` from `nproc` on a hybrid CPU.**

### H.3 Clustering was O(k³)

`consolidate_groups` recomputed every inter-group mean on every merge pass.
cProfile on the real 2,795-genome input:

```
ncalls        tottime   cumtime  function
        6      450.6s   3126.0s  consolidate_groups
289,988,081   476.1s   1213.4s  numpy ix_
289,988,081   382.6s   1130.1s  numpy _mean
```

290 million submatrix means. Block means compose exactly, so carrying **sums**
instead of means makes a merge an O(k) row update:

```
sum(A∪B, C) = sum(A,C) + sum(B,C)
within(A∪B) = within_sum(A) + within_sum(B) + sum(A,B)
```

**3,126 s → 10 s**, output byte-for-byte identical (282 clusters, 2,795 genomes).
Tie-breaking order is preserved deliberately, because Mash distances contain
exact ties (identical genomes at distance 0) and a different tie-break would
silently regroup genomes.

---

## I. Full-scale validation

2,795 genomes, 61 clusters, threshold 0.006425 auto-selected, 22-core laptop:

| metric | value |
|---|---|
| wall time | 8h 20m |
| CPU-hours | 108 |
| `work/` | 264 GB |
| peak RSS | 12.8 GB |
| clusters fully resolved | **61 / 61** |
| grafted tree tips | **2,794 unique, 0 duplicates** |
| `Reference` tips | **0** |
| Gubbins taxon retention | **2,794 in → 2,794 out** |
| `EXCLUDED_TAXON` lines | **0** |
| recombination blocks | **109,588** |

### Still unvalidated

- The A100 profile has never been run on an A100; its
  `SNIPPY_SCATTER maxForks=48` is reasoning, not measurement.
- Gubbins timings at 50 taxa are extrapolated from 27.
- Collections beyond 2,795 genomes.
- `-profile singularity` on this path.
- Bitwise reproducibility: topologies are stable, branch lengths differ in the
  7th decimal (IQ-TREE thread-count changes floating-point reduction order; no
  `-seed` is set).
