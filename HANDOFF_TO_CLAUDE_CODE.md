# Handoff: `wf-assembly-snps-mod` optimization → Claude Code

**Branch:** `perf/low-spec-optimization` (9 commits on `caa3013`) · **PR:** #2, open against `main`
**Repo:** https://github.com/PHemarajata/wf-assembly-snps-mod
**Everything below is pushed.** Nothing lives only in a sandbox.

---

## 1. Why this work exists

The user runs recombination-aware SNP analysis on *Burkholderia pseudomallei* assemblies (~7.2 Mb,
two chromosomes, highly clonal, recombinogenic). Processing 2000+ genomes required a DGX Station
A100. The goal was to make it run on a modest workstation, and the user's opening question was
whether a GPU could help.

**It cannot, and this is settled — do not re-litigate it.** IQ-TREE has never had CUDA support in any
version: the entire `iqtree3` source tree contains zero GPU code (the sole grep hit for `cuda` is a
contributor's name in a header comment), `-h` exposes no GPU option, and `ldd` on the binary links no
CUDA libraries. Same for RAxML-NG, FastTree, and Mash. VeryFastTree's `-ext CUDA` prints a
`with CUDA` banner and *then* exits 1 on the bioconda build — that is almost certainly what the user
tried. GPU-capable phylogenetics (BEAGLE) serves Bayesian MCMC, a different method that would be far
slower here. Moreover BEAGLE's own authors note GPU likelihood parallelism gains little when there
are few unique site patterns, which is exactly the clonal regime. The DGX was helping through CPU
cores and RAM, never its A100s.

## 2. What was wrong, and what the fixes rest on

Two categories. Keep them separate when talking to the user — they care about the distinction.

### Performance (does not change results)

Measured on the user's 112 real assemblies:

| Stage | Before | After | Speedup |
|---|---|---|---|
| Mash front end (112 genomes) | 14.3 s | 1.45 s | 9.9× |
| Column filter (50 taxa × 7.16 Mb) | ~5.0 min | 5.51 s | ~55× |
| `IQTREE_FAST` (50 × 419,329 SNPs) | >40 min, killed unfinished | 13.55 s | >175× |

Root causes: two `O(n²)` Python loops (`df.iterrows()` + `.at[]`, and `.iloc[i,j]` per cell);
`mash dist` never passed `-p` despite `cpus=4`; `IQTREE_FAST` computed a model variable, discarded
it, and ran full ModelFinder (**968 models tested, plain `JC` selected**) to build a starting tree
Gubbins throws away; `-nt AUTO` burned ~57 s/task probing *and* sized to all 22 host cores regardless
of `task.cpus`; Snippy ran every sample serially inside one task.

### Correctness (changes results — six defects)

1. **Star trees written as results.** `GTR+ASC` aborts on any invariant column, which happens
   routinely in clonal alignments (a 20-taxon test with no deliberate constant columns still had
   308). The module caught the failure and wrote `($names);` — **zero internal nodes** — into
   `*.final.treefile`, which grafting consumed as a success. Fixed with a constant-site preflight
   (`bin/asc_preflight.py`) that keeps `+ASC` on variable sites only; the star-tree fallback is gone
   and failures now fail. **The user should check archived results** — `RUNNING.md` has the loop.
2. **Cluster splitting by list order.** Oversized components were split by contiguous Python list
   slices, arbitrary w.r.t. phylogeny. 0/40 pure clusters on synthetic data, 0/4 on real. Now
   average-linkage on the Mash submatrix, gated by `--cluster_split_method` (`order` = legacy).
3. **Column filter discarded the informative sites.** See §4 — the reported magnitude was wrong but
   the defect is real.
4. **Sample→assembly binding by substring match** mis-bound 3 of 6 prefix-colliding IDs.
5. **`--filter-percentage` never set**, so Gubbins' invisible 25% default silently dropped taxa.
6. **Crash on modern pandas.** The original `mash_tab_to_matrix.py` dies with
   `ValueError: underlying array is read-only` under pandas 3.0 — latent because the old modules ran
   unpinned `pip install pandas numpy` at runtime.

## 3. Corrections I had to make to my own earlier claims

Carry these forward; they are the parts most likely to mislead.

- **The 6.3% figure was wrong.** I originally reported the all-ATCG filter left Gubbins seeing 6.3%
  of the genome at 50 taxa. On real data it keeps **89.03%**. My synthetic model assumed
  *independent* per-genome ambiguity; real missingness is strongly **spatially correlated** (same
  accessory regions absent across many genomes), so losses don't compound — independence predicts
  18% survival, observed 89%. **The defect is still real but for a different reason:** the discarded
  columns are disproportionately *variable*. In a 400k-column sample, all-ATCG retained 3,006
  variable sites vs 23,229 under ≤10% (**7.7×**); genome-wide 55,188 → 62,552 (+13.3%).
  The `0.508 → 0.921` Gubbins recall figure is simulation-only — treat as an upper bound.
- **My "0.005 threshold" suggestion was superseded.** It assumed the profile default of 0.028. The
  user actually runs **0.002**. Correct advice is **0.003** (§5).
- **I claimed SNIPPY_SCATTER executed successfully. It did not** — all 4 tasks FAILED with exit 127,
  `snippy: command not found`. Corrected in commit `05aa267`; `RUNNING.md` now carries a per-process
  validation table. Be precise about this: reaching the snippy call site proves the reference-FASTA
  wiring constructs arguments correctly, and nothing more.

## 4. The user's data — measured, not assumed

112 assemblies in `~/Downloads/subset_100/`, 6.96–7.39 Mb, median 91 contigs, 16 complete / 96 draft,
median N content 437 bases (0.006%). No Parsnp length-drops expected.

**Pairwise Mash distances: min 0.00000, median 0.00448, max 0.00700 across 6,216 pairs.** This single
fact drives most threshold reasoning — the population is extremely clonal.

## 5. Threshold guidance (`--mash_threshold`)

| threshold | clusters ≥3 taxa | genomes properly clustered | in merged/dropped bin | modelled Gubbins |
|---|---|---|---|---|
| 0.002 (user's current) | 8 | 75 | 37 | 4.2 min |
| **0.003 (recommended)** | 8 | **96** | **16** | 3.6 min |
| 0.005 | 1 | 110 | 2 | 19.3 min |
| 0.028 / 0.03 (profile defaults) | 1 | 112 | 0 | 20.0 min |

At 0.002, **30 of 38 clusters have <3 taxa** and cannot produce a tree. `merge_singletons=true`
(set in `bp.config`) then pools all 37 genomes into one `merged_small_clusters` bin **with no
similarity criterion** — internal distances up to 0.006985 against a dataset maximum of 0.006995,
4.0× more divergent than a genuine cluster. That bin is not a clade, and Gubbins assumes samples of
limited diversity sharing a recent common ancestor, so its recombination calls there are not
interpretable. I consider this a **seventh defect, not yet fixed** (§8).

0.003 is stable, not a lucky cut point: coverage rises smoothly to 96 genomes across 0.0029–0.0033
before components fuse past ~0.0035. **Recalibrate for any new collection** — `threshold_scan_real.csv`
reproduces the scan. The Gubbins cost column is a model (`0.110·n^1.97`), not a measurement.

## 6. Validation status — be rigorous about this

The development sandbox blocks AF_UNIX socket creation outright, which disables both the Docker
client (`/var/run/docker.sock` does not exist inside it) and Gubbins' `pyjar`. Rootless apptainer
also fails: its re-exec'd child hardcodes `/etc/apptainer/apptainer.conf`, which needs root.
**This is why you have containers and I did not.**

From `pipeline_info/execution_trace_*.txt` on 10 real assemblies:

| Process | Status |
|---|---|
| `INFILE_HANDLING_UNIX` (10), `MASH_SKETCH_BATCH`, `MASH_PASTE`, `MASH_TRIANGLE`, `CLUSTER_GENOMES`, `SELECT_CLUSTER_REPRESENTATIVE` (2) | **COMPLETED** |
| `SNIPPY_SCATTER` (4–5) | **FAILED, exit 127** — binary absent |
| everything downstream | **NEVER RAN** |

**Never validated, in priority order — this is your job:**

1. **Snippy end to end.** The scatter/gather restructuring was validated with a mock; `snippy-core`
   has never run. Also unquantified: scattering makes snippy re-index the reference once per
   *sample* rather than once per cluster.
2. **Gubbins.** Recall/precision figures are simulation-only. Check `EXCLUDED_TAXON:` lines.
3. **IQ-TREE ASC preflight on real Gubbins output.** Logic is tested; the real input shape is not.
4. **`publish_dir_mode = 'link'`** — needs `outdir` and `work/` on one filesystem.
5. **Scale.** Nothing above 112 genomes. `low_spec` memory figures are budgeted, not observed.
6. **The branch-support multiplier is unresolved** — I measured 1.66× and 8.0× under different
   models and could not determine which predicts production cost. `iqtree_support` defaults **off**;
   **turn it on for anything publishable**, since support values are how a reader judges a clade.

## 7. Traps that cost me time

- **`-profile bp,low_spec` silently changes science.** Nextflow applies profiles left to right;
  `low_spec` sets `max_cluster_size`, `iqtree_support`, `max_column_missingness`,
  `gubbins_filter_percentage`. With `-profile low_spec` alone, `mash_threshold` falls back to
  `params.config`'s 0.03 — above the user's maximum observed distance. **Always pass
  `--mash_threshold` explicitly.**
- **snippy 4.6.0 parses `samtools --version` as a float** and compares to `1.3`, so samtools **1.23
  evaluates as < 1.3** and snippy aborts. `samtools=1.9` is pinned in `environment.yml` and the
  module's conda directive. Do not bump it casually.
- **`-profile conda` was advertised but resolved nothing** — 0 of 41 modules declared a `conda`
  directive. All 11 on the active path now do.
- **`nextflow config` does not run the schema validator.** Four bugs reached the user because I
  validated that way. Only a real `nextflow run` catches them.
- **The user's `java` is Java 11** from a miniforge3 base that auto-activates. Nextflow needs 17–24;
  `JAVA_HOME` overrides `PATH` (verified). Launch from a clean shell.
- **`.nextflow.log` was tracked** via a `!.nextflow.log` negation in `.gitignore`, breaking
  `git pull` for anyone who had run the pipeline. Fixed in `125a117`.

## 8. Open work, in the order I would do it

1. **Run the ~10-genome smoke test with Docker** (`RUNNING.md` Step 3) and fix what surfaces. This
   is the highest-value action available — items 1–4 of §6 all close here.
2. **Similarity-aware `merge_singletons`** (the seventh defect, §5). Currently pools unrelated
   genomes into one incoherent bin. Should group leftovers into coherent bins or mark them
   unanalysable rather than fabricating a pseudo-cluster. **Not started.**
3. **GitHub Actions CI.** A smoke test on 3–5 tiny genomes per push would have caught four of my
   bugs. The user's goal is independent use by others; this matters more than it looks.
4. **Resolve the branch-support multiplier** (§6.6) with a real post-Gubbins alignment.
5. **Scale test** at 500 → 2000 genomes, watching RSS against the `low_spec` budget.

## 9. Working with this user

An experienced bioinformatician who runs what you give them immediately and reports exact errors —
four of my bugs came back within minutes. They ask good scope questions before acting
("should I deactivate conda first?"), so **flag environment assumptions explicitly**. They care that
changes are scientifically defensible and reversible: every behaviour change is gated by a param
whose legacy value restores prior output, and `MIGRATION_NOTES.md` documents the mapping. They are
moving to Claude Code specifically for container access. Their end goal is a repo **others can use
independently** — weigh portability and documentation accordingly.

## 10. Reference documents on the branch

- `RUNNING.md` — three execution modes, per-process validation table, star-tree check, legacy-flag
  invocation that reproduces prior results
- `AUDIT_REPORT.md` — per-defect evidence; §G covers the threshold analysis and Docker diagnosis;
  the addendum carries the 6.3%→89% correction
- `MIGRATION_NOTES.md` — which defaults alter output and how to revert each
- `environment.yml` — native install, no container engine required

**One thing I would tell the user again:** any `.final.treefile` from the old code with zero internal
nodes carries no phylogenetic information but looks like a normal result. If they have published or
archived output from before this branch, it is worth checking.
