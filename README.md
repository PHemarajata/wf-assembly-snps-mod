# wf-assembly-snps-mod

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.04.3-23aa62.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/docker-blue.svg)](https://www.docker.com/)
[![Singularity](https://img.shields.io/badge/singularity-blue.svg)](https://sylabs.io/docs/)

Recombination-aware SNP phylogenetics for large bacterial assembly collections,
built for *Burkholderia pseudomallei* (~7.2 Mb, two chromosomes, highly clonal
and recombinogenic).

The workflow clusters genomes by Mash distance, builds a whole-genome alignment
per cluster, masks recombination with Gubbins, infers a per-cluster ML tree with
IQ-TREE, then grafts the cluster trees onto a backbone tree of cluster
representatives.

> Modified fork of [bacterial-genomics/wf-assembly-snps](https://github.com/bacterial-genomics/wf-assembly-snps).

**Status:** validated end to end on **2,795 real assemblies** on a 22-core
laptop — 61 clusters, 8h20m, 108 CPU-hours, every genome retained through
Gubbins, all output checks passed. See
[Validation](#validation-what-has-actually-been-run).

---

## Contents

- [Quick start](#quick-start)
- [Read this first: silent-success failures](#read-this-first-silent-success-failures)
- [Defects found and fixed](#defects-found-and-fixed)
- [Improvements](#improvements)
- [How the workflow works](#how-the-workflow-works)
- [Choosing a profile](#choosing-a-profile)
- [Choosing `--mash_threshold`](#choosing---mash_threshold)
- [Resource budget](#resource-budget)
- [Validating your output](#validating-your-output)
- [Troubleshooting](#troubleshooting)
- [Parameter reference](#parameter-reference)
- [Output structure](#output-structure)
- [Validation: what has actually been run](#validation-what-has-actually-been-run)
- [Standalone tools](#standalone-tools)

---

## Quick start

### Prerequisites

```bash
nextflow -version    # 23.04+; tested on 25.04.6
java -version        # MUST be 17-24
docker info          # or singularity
```

> **Java is the most common first failure.** An active conda/miniforge base
> usually puts Java 11 on `PATH`, and `JAVA_HOME` wins over `PATH`. Set it:
>
> ```bash
> export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64   # adjust to your JDK
> export PATH="$JAVA_HOME/bin:$PATH"
> ```

### On a DGX Station A100 (128 cores, 512 GB)

```bash
git clone https://github.com/PHemarajata/wf-assembly-snps-mod.git
cd wf-assembly-snps-mod

export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
export PATH="$JAVA_HOME/bin:$PATH"

nextflow run . \
    -profile dgx_station_a100_updated,docker \
    --recombination_aware_mode true \
    --input  /path/to/assemblies \
    --outdir results_run \
    --mash_threshold auto \
    --max_cluster_size 50 \
    --gubbins_tree_builder rapidnj \
    --publish_dir_mode link
```

**Run a 10-genome smoke test first.** Copy ten assemblies into a directory and
run the same command against it — a few minutes, and it exercises every process
in the pipeline.

> **The A100's GPUs cannot accelerate this workflow.** No tree builder on this
> path has CUDA support: IQ-TREE contains no GPU code, and RAxML-NG, FastTree
> and Mash have none either. VeryFastTree's `-ext CUDA` prints a "with CUDA"
> banner and then exits 1 on the bioconda build. What makes that machine fast
> here is 64 cores and 512 GB of RAM. There is no GPU switch to find — this was
> investigated thoroughly and is settled.

### On a 22-core / 64 GB laptop

```bash
nextflow run . -profile local_workstation_rtx4070,docker \
    --recombination_aware_mode true \
    --input /path/to/assemblies --outdir results_run \
    --mash_threshold auto --max_cluster_size 50 \
    --gubbins_tree_builder rapidnj --publish_dir_mode link
```

Full walkthrough in [`RUNNING.md`](RUNNING.md).

---

## Read this first: silent-success failures

This branch began as a performance exercise — make a workload that needed a DGX
run on a workstation — and became a correctness audit. **Several defects produced
scientifically wrong output while reporting success.** Those matter far more than
the speedups.

The original code repeatedly caught a failure, wrote an empty or placeholder
output, and exited 0. Nextflow showed a green checkmark and the pipeline
continued on meaningless data.

Before this branch:

- **Gubbins** could fail on *every* cluster and the run would still "succeed",
  producing a global tree with **no recombination masking at all** — and nothing
  in the output indicating it. On a recombinogenic organism, that is the one
  error this workflow exists to prevent.
- **IQ-TREE** wrote star trees (`(a,b,c,d);` — zero internal nodes) when `+ASC`
  aborted on constant columns. These look like ordinary Newick files and carry no
  phylogenetic information.

Both now fail loudly with diagnostics. **This is why the measurements in this
README are trustworthy and older results may not be.**

> ### ⚠️ If you have results from before this branch, check them
>
> ```bash
> # 1. Star trees: 0 internal nodes means no phylogenetic information
> for f in <old_results>/**/*.final.treefile; do
>   echo "$(( $(tr -cd '(' < "$f" | wc -c) - 1 )) internal nodes  $f"
> done | sort -n | head
>
> # 2. Empty recombination output means Gubbins failed silently
> find <old_results> -name "*.recombination_predictions.gff" -empty
> ```

---

## Defects found and fixed

Each was found by running the workflow and reading its output, not by code
inspection. Evidence and measurements are in the commit messages and
[`AUDIT_REPORT.md`](AUDIT_REPORT.md).

### Scientific correctness

| # | Defect | Consequence |
|---|---|---|
| 1 | `GUBBINS_CLUSTER` wrote empty placeholders and exited 0 on failure | **Runs with zero recombination masking reported as successful** |
| 2 | `gubbins_tree_builder = 'veryfasttree'` is not a value Gubbins accepts | Gubbins exited 2 on **every** cluster under `low_spec`; with #1 hiding it, that profile had never produced a real recombination analysis |
| 3 | `+ASC` aborts on constant columns; the module caught it and wrote a **star tree** | Trees with zero internal nodes published as results |
| 4 | snippy-core emits the reference as an extra taxon `Reference`, duplicating the cluster medoid | Every cluster carried a duplicated genome; the grafted tree held one `Reference` tip per cluster, each a *different* organism |
| 5 | `COLLECT_REPRESENTATIVES` derived sequence IDs from a filename that had become a constant | Every backbone tip would have been labelled `representative`, so no tip could match its cluster during grafting |
| 6 | Oversized clusters split by arbitrary list-order slices | Clusters that ignore phylogeny (0/40 pure on synthetic data) |
| 7 | Gubbins' `--filter-percentage` was never set, so its invisible 25% default applied | Taxa silently dropped from alignments |

### Environment and portability

| # | Defect | Consequence |
|---|---|---|
| 8 | `docker.runOptions` written as a bare key inside `profiles { docker { … } }` — which is a *profile* named `docker`, not the docker config **scope** | Containers ran as **root**; outputs were root-owned and `publish_dir_mode='link'` failed with `EPERM` (`fs.protected_hardlinks=1` forbids hardlinking a file you don't own) |
| 9 | Docker defaults `/dev/shm` to **64 MB** | Gubbins' `pyjar` `multiprocessing.Pool` overran it on large clusters: `OSError: [Errno 28] No space left on device` — which reads like a full disk and is not one. Fixed with `--shm-size=2g` |
| 10 | Assemblies headed `>SAMPLE.fasta 1`, `>SAMPLE.fasta 2` share one sequence **ID** (the first whitespace token) | snippy refuses the reference: `Duplicate sequence … in <ref>`. 17 of 112 test assemblies affected |
| 11 | The IQ-TREE container ships **no python**, but the module shelled out to `asc_preflight.py` | The ASC preflight worked under `-profile conda` and could never work in a container |
| 12 | Profiles referenced processes that no longer exist (`SNIPPY_ALIGN`, `MASH_SKETCH`, `MASH_DIST`, `SKA_BUILD`, `GPU_ACCELERATED`) | The A100 profile matched **4 of 21** processes; `SNIPPY_SCATTER` — 89% of runtime — had no resource block and fell through to a generic label |
| 13 | Allocations were budgeted, not measured, and 12–360× over | Nextflow packs tasks by *declared* memory, so over-allocation **directly caps concurrency**: 22 GB declared for a process using 0.89 GB meant only two could run at once |

### Documented behaviour that turned out to be wrong

- **Profile ordering.** Earlier notes claimed `-profile bp,low_spec` and
  `-profile low_spec,bp` behave differently. They do not — precedence follows the
  order profiles are **defined** in `nextflow.config`, not the order you list
  them. Both give `max_cluster_size=25`. **Composing `low_spec` with a larger
  machine profile silently clamps you to 16 cores / 30 GB.**
- **The `publish_dir_mode='link'` failure** was documented as a same-filesystem
  problem. It was file **ownership** (defect #8).
- **Reference re-indexing** was assumed to be the snippy bottleneck. Measured,
  `bwa index` is 6% of a task; `bwa mem` is 63%.

---

## Improvements

### `--mash_threshold auto` (new)

A fixed threshold does not transfer between collections. `0.003`, calibrated on
112 genomes, silently discards **203 of 2,795** genomes on a wider set because
they land in components of <3 taxa and cannot produce a tree.

`--mash_threshold auto` sweeps candidates and selects the one analysing the most
genomes **without producing a cluster spanning more than 90% of the collection's
own largest pairwise distance**. The bound is *relative*, which is what lets it
transfer between datasets.

Writes `Summaries/threshold_sweep.tsv` and `Summaries/chosen_threshold.txt`.

### Performance

| change | effect |
|---|---|
| Batched, vectorised Mash front end | 2,795 genomes: sketch 43 s, triangle 41 s |
| `consolidate_groups` O(k³) → O(k²) | **3,126 s → 10 s**, byte-identical output |
| Snippy scatter/gather | one task per cluster-sample instead of serial per cluster |
| Resource blocks sized from measurement | ~4× Gubbins throughput on the same hardware |

### Correct-by-construction resource profiles

`local_workstation_rtx4070.config` and `dgx_station_a100_updated.config` carry
measured peak RSS and measured concurrency knees with the reasoning inline.
Clone and run; no tuning required.

---

## How the workflow works

```
INFILE_HANDLING_UNIX      normalise input, make contig IDs unique
   ↓
MASH_SKETCH_BATCH → MASH_PASTE → MASH_TRIANGLE      pairwise distances
   ↓
CLUSTER_GENOMES           single-linkage components, split to max_cluster_size
   ↓                      by average linkage on the Mash submatrix
SELECT_CLUSTER_REPRESENTATIVE     medoid per cluster
   ↓
SNIPPY_SCATTER (per cluster-sample) → SNIPPY_CORE_GATHER
   ↓                                  whole-genome alignment per cluster
KEEP_INVARIANT_ATCG       column filter (retains invariant sites for Gubbins)
   ↓
IQTREE_FAST               cheap starting tree (Gubbins discards it)
   ↓
GUBBINS_CLUSTER           recombination detection and masking
   ↓
ASC_PREFLIGHT → IQTREE_ASC        final per-cluster ML tree
   ↓
COLLECT_REPRESENTATIVES → BUILD_BACKBONE_TREE       parsnp over representatives
   ↓
GRAFT_TREES               cluster trees grafted onto the backbone
```

**Why cluster at all:** Gubbins assumes samples of limited diversity sharing a
recent common ancestor. Run across a globally diverse collection that assumption
breaks and its recombination calls stop being interpretable. Clustering keeps
each Gubbins run inside a coherent group — which is also why the threshold
matters so much.

---

## Choosing a profile

| profile | target | notes |
|---|---|---|
| `dgx_station_a100_updated` | 128 cores / 512 GB | blocks updated for the current pipeline; see caveat in [Validation](#validation-what-has-actually-been-run) |
| `local_workstation_rtx4070` | 22 cores / 64 GB | fully measured on that hardware |
| `low_spec` | 16 cores / 32 GB | conservative; **do not compose with a bigger profile** |

Always add a container profile: `,docker` or `,singularity`.

> **Do not list two machine profiles together.** Precedence follows definition
> order in `nextflow.config`, not command-line order, so the result is usually
> not what you intend.

**`low_spec` also sets science parameters** (`max_cluster_size`,
`gubbins_tree_builder`, `max_column_missingness`, `gubbins_filter_percentage`,
`iqtree_support`, `publish_dir_mode`). Other profiles do not — so pass
`--gubbins_tree_builder rapidnj` and `--publish_dir_mode link` explicitly when
using a different profile.

---

## Choosing `--mash_threshold`

Use `auto`, then **read `Summaries/threshold_sweep.tsv`**.

Measured on 2,795 *B. pseudomallei* genomes whose total pairwise span is
**0.01088**:

| threshold | clusters ≥3 | dropped | within-cluster max | verdict |
|---|---|---|---|---|
| 0.003 | 109 | 203 | 0.00698 | coherent, but discards 7% of the data |
| 0.005 | 67 | 58 | 0.00946 (87%) | accepted |
| **0.006425** | **61** | **1** | **0.00946 (87%)** | **selected** |
| 0.006549 | 61 | 0 | 0.01061 (98%) | rejected — one cluster spans 98% of the entire collection |

**Acceptance is not monotonic.** Raising the threshold changes which component a
borderline genome joins, and average-linkage re-splitting can then *lower* a
cluster's span, so accepted and rejected values interleave. The chosen value can
sit in a narrow accepted window — 0.006425 lies between rejected neighbours at
0.006300 and 0.006549 — which means the clustering is genuinely sensitive there.
Worth knowing before building a phylogeny on it.

Override with an explicit value when the sweep shows something the criterion
cannot see.

---

## Resource budget

Measured on a completed 2,795-genome run (61 clusters) on a 22-core laptop:

| resource | measured |
|---|---|
| wall time | **8h 20m** |
| CPU-hours | 108 |
| `work/` for one clean run | **264 GB** |
| published results | 20 GB (hard links — no extra disk) |
| peak RSS | **12.8 GB** (`BUILD_BACKBONE_TREE`) |

**Keep 300 GB free.** On a DGX with more cores, expect proportionally less wall
time; disk is unchanged.

Per-task: `SNIPPY_SCATTER` ~62 MB, `SNIPPY_CORE_GATHER` ~630 MB,
`BUILD_BACKBONE_TREE` ~1.2 GB.

> **`-resume` is not free after a config change.** Editing `docker.runOptions` —
> or anything else feeding the task hash — invalidates cached tasks. A mid-run
> `--shm-size` fix re-ran all 2,794 alignments and left `work/` holding two
> copies of the alignment stage (510 GB, filling the disk to 99%). Batch
> container-option changes **before** a long run. `rm -rf work` does not delete
> published results when `publish_dir_mode='link'`.

---

## Validating your output

Exit code 0 is not sufficient — that is the central lesson of this fork. Run
these against every result set:

```bash
R=results_run

# 1. No star trees. A fully resolved tree has (parens - 1) == (tips - 3).
for f in $R/Clusters/*/*.final.treefile; do
  p=$(tr -cd '(' < "$f" | wc -c); t=$(( $(tr -cd ',' < "$f" | wc -c) + 1 ))
  echo "$((p-1)) internal / $t tips   $(basename $(dirname $f))"
done

# 2. Grafted tree: no duplicate tips, no leftover "Reference"
G=$R/Final_Results/global_grafted.treefile
tr ',()' '\n\n\n' < $G | sed "s/:.*//; s/'//g" | grep -vE '^$|^[0-9.]+$' \
  | sort | uniq -d          # must print nothing
grep -c Reference $G        # must be 0

# 3. Gubbins kept every taxon
grep -rh "Taxa in alignment" $R/Clusters/*/Gubbins/*.diagnostics.log
grep -rh "EXCLUDED_TAXON"    $R/Clusters/*/Gubbins/*.diagnostics.log

# 4. Recombination was actually detected
cat $R/Clusters/*/Gubbins/*.recombination_predictions.gff | grep -vc '^#'
```

Reference 2,795-genome run: **61/61 fully resolved, 2,794 unique tips, zero
duplicates, zero `Reference` tips, 2,794 taxa in → 2,794 out, 109,588
recombination blocks.**

---

## Troubleshooting

| symptom | cause | fix |
|---|---|---|
| `Cannot find Java or it's a wrong version` | conda base puts Java 11 on `PATH` | `export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64` |
| `OSError: [Errno 28] No space left on device` in Gubbins | `/dev/shm` is 64 MB — **not** the disk | already fixed via `--shm-size=2g`; verify with `nextflow config -profile <p>,docker \| grep runOptions` |
| `Duplicate sequence … in <ref>` | input contig IDs not unique | handled in `INFILE_HANDLING_UNIX`; if it recurs, inspect your headers |
| `Failed to publish file … [link]` `Operation not permitted` | containers running as root | check `runOptions` contains `-u $(id -u):$(id -g)` |
| `no process matching config selector: X` | stale selector in a profile | harmless, but process `X` has no resource block |
| Machine mostly idle mid-run | over-declared memory capping concurrency | compare declared memory against `peak_rss` in `pipeline_info/execution_trace_*.txt` |
| Gubbins reports success with an empty GFF | **should be impossible now** | if you see it, you are not on this branch |

**When a task fails**, the error names its work directory. Inside it:
`.command.err`, `.command.out`, `.command.sh`, and for Gubbins the published
`*.diagnostics.log`.

---

## Parameter reference

### Core

| parameter | default | meaning |
|---|---|---|
| `--input` | — | directory of assemblies, or a `sample,file` samplesheet |
| `--outdir` | — | results directory |
| `--recombination_aware_mode` | `false` | **set `true`** for this workflow |
| `--mash_threshold` | `0.03` | number, or `auto` (recommended) |
| `--max_cluster_size` | `100` | hard cap; 50 is the tested production value |
| `--alignment_method` | `snippy` | `snippy` \| `ska` \| `parsnp` — see warning below |
| `--publish_dir_mode` | `copy` | `link` avoids duplicating results on disk |
| `--min_input_filesize` | `45k` | minimum input file size |

### Clustering

| parameter | default | meaning |
|---|---|---|
| `--auto_threshold_coherence` | `0.90` | max fraction of collection span a cluster may span |
| `--auto_threshold_grid` | `null` | comma-separated candidates; default is 15 data quantiles |
| `--auto_refine` (script-level) | `12` | uniform refinement steps above the best accepted threshold |
| `--cluster_split_method` | `similarity` | `order` reproduces legacy output |
| `--merge_singletons` | `false` | pools small clusters into one bin — **not a clade**; leave off |
| `--mash_sketch_size` | `10000` | Mash sketch size |
| `--mash_kmer_size` | `21` | Mash k-mer size |
| `--mash_batch_size` | `200` | genomes per sketch task |

### Recombination and trees

| parameter | default | meaning |
|---|---|---|
| `--gubbins_tree_builder` | `iqtree` | **use `rapidnj`**; `veryfasttree` is invalid and fails |
| `--gubbins_first_tree_builder` | `rapidnj` | initial tree builder |
| `--gubbins_iterations` | `3` | Gubbins iterations |
| `--gubbins_filter_percentage` | `25` | per-taxon missingness; `100` disables taxon dropping |
| `--max_column_missingness` | `0.10` | column filter; `0.0` is the legacy all-ATCG rule |
| `--iqtree_support` | `false` | **turn on for anything publishable** (UFBoot/SH-aLRT) |
| `--iqtree_asc_model` | `GTR+ASC` | final-tree model |
| `--drop_reference_taxon` | `null` | drops snippy's duplicate `Reference`; `false` restores legacy |

Every parameter that alters scientific output is annotated in
`conf/params.config` with the measurement behind its default, and
[`MIGRATION_NOTES.md`](MIGRATION_NOTES.md) maps each to the value restoring prior
behaviour.

### ⚠️ `--alignment_method ska` undercalls recombination — do not use it here

SKA is 2.7× faster and **wrong for this workflow**. A split k-mer cannot call a
SNP whose flank contains another SNP, so SKA recovers only ~11% of SNPs within
10 bp of a neighbour while matching snippy beyond ~100 bp. Dense SNP clusters are
exactly the signal Gubbins keys on.

Measured on 112 genomes: **54% fewer recombination blocks** (2,388 vs 5,148) and
different topologies in 3 of the 4 informative clusters. The workflow warns when
you select it with `run_gubbins` on. SKA remains fine for recombination-free uses
such as quick clustering or distance estimates.

---

## Output structure

```
results_run/
├── Summaries/
│   ├── clusters.tsv                  cluster_id → sample_id
│   ├── chosen_threshold.txt          the threshold auto-selection picked
│   ├── threshold_sweep.tsv           every candidate, with the reason for rejection
│   ├── cluster_phylogeny_summary.csv per-cluster Gubbins/IQ-TREE status + tier
│   └── cluster_summary.txt
├── Clusters/
│   └── cluster_<id>/
│       ├── <id>.core.full.aln        whole-genome alignment
│       ├── <id>.final.treefile       per-cluster ML tree
│       ├── <id>.asc_preflight.txt    ASC decision for this cluster
│       └── Gubbins/
│           ├── <id>.recombination_predictions.gff
│           ├── <id>.filtered_polymorphic_sites.fasta
│           └── <id>.diagnostics.log  taxa in/out, exclusions, exit code
├── backbone.treefile                 tree over cluster representatives
├── representatives.fa
├── cluster_representatives.tsv
├── Final_Results/
│   └── global_grafted.treefile       ← the final tree
└── pipeline_info/
    ├── execution_trace_*.txt         per-task runtime and peak_rss
    ├── execution_report_*.html
    └── execution_timeline_*.html
```

`pipeline_info/execution_trace_*.txt` is the file to read when tuning resources:
it carries `realtime` and `peak_rss` per task.

---

## Validation: what has actually been run

| scale | status |
|---|---|
| 10 genomes | ✅ full pipeline, repeatedly |
| 112 genomes | ✅ full pipeline; snippy vs SKA method comparison |
| **2,795 genomes** | ✅ **full pipeline, 61 clusters, all output checks passed** |

Also exercised: cluster splitting (fires on ≥50-taxa components), auto-threshold
selection, contig-ID normalisation (17/112 files), Gubbins on 50-taxa clusters,
`BUILD_BACKBONE_TREE` over 138 representatives.

### Not yet validated — please help close these

- **The A100 profile has not been run on an A100.** Resource blocks were
  corrected and sized from laptop measurements, but `SNIPPY_SCATTER maxForks=48`
  is *reasoning, not measurement* — the laptop's knee of 12 was a hybrid-CPU
  artifact (P-cores vs E-cores) and a homogeneous EPYC should scale better.
  **Please re-measure it**, using the maxForks sweep recorded in
  `conf/profiles/local_workstation_rtx4070.config` as the method: run the same
  input at two or three `maxForks` values and compare `realtime` per task in the
  execution trace.
- Gubbins timings at 50 taxa are extrapolated from measurements at 27.
- Collections beyond 2,795 genomes.
- `-profile singularity` on the recombination-aware path.
- Bitwise reproducibility: topologies are stable across runs, but branch lengths
  differ in the 7th decimal because IQ-TREE's thread count changes
  floating-point reduction order and no `-seed` is set.

---

## Standalone tools

### Downsampling contextual genomes

When you have far more contextual (public) assemblies than you can analyse,
`bin/downsample_contextual.py` builds a balanced, de-redundant input **before**
running the pipeline. It collapses near-identical genomes within each country by
Mash distance, then enforces a per-country floor and cap. Study isolates are
always kept in full. Output is a `sample,file` samplesheet that `--input`
accepts.

```bash
./bin/downsample_contextual.py \
  --contextual-dir contextual_fasta --cdc-dir cdc_fasta \
  --metadata metadata.tsv \
  --derep-threshold 0.001 --target-total 1000 \
  --min-per-country 3 --max-per-country 200 --seed 42 \
  --outdir downsample_out
```

Key options: `--sweep` (report group counts across thresholds and exit),
`--mash-dist` (reuse cached distances), `--unmatched drop|unknown`.
Outputs `samplesheet.csv`, `selection_report.tsv`, `country_summary.tsv`.

### Standalone tree grafting

If the in-pipeline grafting step fails, `graft_trees.py` performs the same
leaf-expansion grafting outside Nextflow with detailed logging.

```bash
pip install biopython

./graft_trees.py \
  --backbone results_run/backbone.treefile \
  --clusters 'results_run/Clusters/**/cluster_*.final.treefile' \
  --reps     results_run/cluster_representatives.tsv \
  --out-tree results_run/final_grafted.treefile \
  --report   results_run/grafting_report.txt \
  --log      results_run/grafting_log.txt
```

Add `--rename-conflicts` for duplicate tip labels, `--dry-run` to plan only.

---

## Further reading

| document | contents |
|---|---|
| [`RUNNING.md`](RUNNING.md) | step-by-step execution, per-machine walkthroughs |
| [`AUDIT_REPORT.md`](AUDIT_REPORT.md) | per-defect evidence and measurements |
| [`MIGRATION_NOTES.md`](MIGRATION_NOTES.md) | which defaults changed output, how to revert |
| `conf/params.config` | every parameter with the measurement behind its default |
| `conf/profiles/*.config` | resource blocks with measured RSS and concurrency knees |

---

## License

Apache License 2.0 — see [LICENSE](LICENSE).

## Acknowledgments

- Original [bacterial-genomics/wf-assembly-snps](https://github.com/bacterial-genomics/wf-assembly-snps)
- [nf-core](https://nf-co.re/) community for workflow conventions
- Mash, snippy, Gubbins, IQ-TREE, Parsnp, SKA2 developers

## Support

Open an [issue](https://github.com/PHemarajata/wf-assembly-snps-mod/issues).

**Citation:** cite the upstream workflow and the underlying tools.
