# RUNNING.md — running `wf-assembly-snps-mod` on your own machine

Written for the `perf/low-spec-optimization` branch. Goal: anyone can clone the repo and run it
without needing anything explained to them in person.

---

## Which execution mode to pick

The workflow needs ~9 external tools. You have three ways to supply them; pick one.

| Mode | Command suffix | When to use |
|---|---|---|
| **Docker** (recommended) | `-profile low_spec,docker` | You have a working Docker daemon. Most reproducible. |
| **Singularity/Apptainer** | `-profile low_spec,singularity` | HPC without root. |
| **Conda, no containers** | `-profile low_spec` + `conda env create -f environment.yml` | No container engine available. |

All three now work. Before this branch, `-profile conda` was advertised in `nextflow.config` but
**no module declared a `conda` directive**, so it silently resolved nothing. All 11 modules on the
recombination-aware path now declare both `conda` and `container`.

---

## Step 0 — prerequisites

```bash
nextflow -version     # needs >= 23.04 (tested on 25.04.6)
java -version         # 17 recommended
docker info           # only if using -profile docker
```

If Nextflow is missing: `curl -s https://get.nextflow.io | bash && sudo mv nextflow /usr/local/bin/`

---

## Step 1 — get the branch

```bash
git clone https://github.com/PHemarajata/wf-assembly-snps-mod.git
cd wf-assembly-snps-mod

# from the bundle:
git fetch /path/to/wf-low-spec-optimization.bundle \
    perf/low-spec-optimization:perf/low-spec-optimization
git checkout perf/low-spec-optimization

chmod +x bin/*.py          # Nextflow puts bin/ on PATH; scripts must be executable
```

---

## Step 2a — the run, on a 22-core / 64 GB laptop (Core Ultra 9 185H + RTX 4070)

If you have this exact machine, use its profile and let the workflow pick the
threshold. Nothing needs tuning; every resource number in
`conf/profiles/local_workstation_rtx4070.config` is measured on this hardware.

```bash
export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
export PATH="$JAVA_HOME/bin:$PATH"

nextflow run . \
    -profile local_workstation_rtx4070,docker \
    --recombination_aware_mode true \
    --input  /path/to/assemblies \
    --outdir results_run \
    --mash_threshold auto \
    --max_cluster_size 50 \
    --gubbins_tree_builder rapidnj \
    --publish_dir_mode link \
    -resume
```

- **`JAVA_HOME` is required.** A miniforge base environment puts Java 11 on
  `PATH` and Nextflow needs 17–24; `JAVA_HOME` wins over `PATH`. Launch from a
  clean shell or set it as above.
- **`--mash_threshold auto`** sweeps candidate thresholds and picks one from the
  collection's own distance distribution, writing `Summaries/threshold_sweep.tsv`
  and `Summaries/chosen_threshold.txt`. A fixed threshold does not transfer
  between collections: 0.003 was derived on 112 genomes and drops 203 of 2795 on
  a wider set.

  **Read the sweep; do not take the number on faith.** Acceptance is not
  monotonic in the threshold — raising it changes which component a borderline
  genome lands in, and average-linkage re-splitting can then *lower* a cluster's
  span, so accepted and rejected values interleave. On the 2795-genome set the
  chosen 0.006425 (1 genome dropped) sits between rejected neighbours at 0.006300
  and 0.006549. That means the clustering is genuinely sensitive near the chosen
  value, which is worth knowing before you build a phylogeny on it. The sweep
  records every candidate with its coverage, its largest within-cluster distance,
  and why it was rejected.

  On the same collection this beat a careful manual choice (0.006, 5 dropped) at
  identical coherence — but it is a starting point for judgement, not a
  substitute for it. Override with an explicit `--mash_threshold <value>` when the
  sweep tells you something the criterion cannot see.
- **`--gubbins_tree_builder rapidnj`** is set explicitly because this profile does
  not own that parameter; without it you inherit `iqtree` from `params.config`.
- **`--publish_dir_mode link`** hard-links results instead of copying, so a run
  costs no extra disk. It needs `outdir` and `work/` on one filesystem.

Do **not** compose this with `low_spec`. Profile precedence follows the order the
profiles are DEFINED in `nextflow.config`, not the order you list them, so
`low_spec` always wins and would clamp this machine to 16 cores / 30 GB. (The
claim in `HANDOFF_TO_CLAUDE_CODE.md` §7 that reordering flips which one wins is
wrong: `-profile bp,low_spec` and `-profile low_spec,bp` both give
`max_cluster_size=25`.)

### Budget at ~2800 genomes on this laptop

MEASURED on a completed 2795-genome run (61 clusters, threshold 0.006425 chosen
by `auto`), not estimated:

| resource | measured |
|---|---|
| `work/` for ONE clean run | **264 GB** |
| published `results/` (hardlinked, no extra disk) | 20 GB |
| peak RSS | **12.8 GB** (`BUILD_BACKBONE_TREE`) |
| wall time | **~8–14 h** |
| CPU-hours | ~108 |

**Keep 300 GB free**, not 250. Per-task: `SNIPPY_SCATTER` ~62 MB × one per
genome, `SNIPPY_CORE_GATHER` ~630 MB × one per cluster, `BUILD_BACKBONE_TREE`
~1.2 GB.

> **`-resume` is not free after a config change.** Editing `docker.runOptions`
> (or anything else that feeds the task hash) invalidates cached tasks and
> re-runs them. On this run a mid-flight `--shm-size` fix re-ran all 2794
> alignments and left `work/` holding TWO copies of the alignment stage —
> 510 GB, which filled the disk to 99%. Batch container-option changes BEFORE a
> long run, and delete `work/` between full runs. Published results are hard
> links, so `rm -rf work` does **not** delete them.

---

## Step 2 — the run

For 112 genomes in `~/Downloads/subset_100/`, on a 16-core / 32 GB workstation:

```bash
nextflow run . \
    -profile low_spec,docker \
    --recombination_aware_mode true \
    --input  ~/Downloads/subset_100 \
    --outdir results_subset100 \
    --mash_threshold 0.003 \
    --max_cluster_size 50 \
    --snp_package parsnp \
    --recombination gubbins \
    -resume
```

**Why `--mash_threshold 0.003` and not the 0.002 you were using:** measured on your 112 genomes,
0.002 leaves 37 of them in clusters too small to build a tree, which `merge_singletons` then dumps into
one incoherent bin. 0.003 gives the same number of analysable clusters (8) but covers 96 genomes
instead of 75. See `AUDIT_REPORT.md` §G. Re-derive it for a new collection with:

```bash
mash sketch -p 8 -s 10000 -k 21 -o all $(ls /path/to/assemblies/*.fasta)
mash triangle -p 8 all.msh > dist.phylip
python3 bin/cluster_mash.py dist.phylip clusters.tsv --threshold 0.003 --max-cluster-size 50
```

### No container engine

```bash
conda env create -f environment.yml
conda activate wf-assembly-snps
nextflow run . -profile low_spec --recombination_aware_mode true \
    --input ~/Downloads/subset_100 --outdir results --mash_threshold 0.003
```

---

## Step 3 — validate before committing to the full set

Run on ~10 genomes first. Three things were never validated in development, and this is where they
surface — all three need a container engine or a full conda env, which the development sandbox lacked
(it blocks Unix-domain sockets outright, which disables both the Docker client and Gubbins' `pyjar`).

1. **Real snippy end to end.** The scatter/gather restructuring was validated with a mock snippy;
   the channel topology is proven, real `snippy-core` behaviour is not.
2. **Gubbins.** Confirm it completes and check `EXCLUDED_TAXON:` lines in the per-cluster
   diagnostics — Gubbins' `--filter-percentage` silently drops taxa above 25% missingness.
3. **`publish_dir_mode = 'hardlink'`.** Requires `outdir` and `work/` on the same filesystem;
   Nextflow falls back to copying otherwise.

```bash
mkdir -p /tmp/ten && ls ~/Downloads/subset_100/*.fasta | head -10 | xargs -I{} cp {} /tmp/ten/
nextflow run . -profile low_spec,docker --recombination_aware_mode true \
    --input /tmp/ten --outdir results_smoke --mash_threshold 0.003 -resume
```

Check afterwards:

```bash
# no cluster should produce a star tree (0 internal nodes) -- the old code emitted these silently
for f in results_smoke/**/*.final.treefile; do
  echo "$(( $(tr -cd '(' < "$f" | wc -c) - 1 )) internal nodes  $f"
done | sort -n | head

grep -r "EXCLUDED_TAXON" results_smoke/ | head        # taxa Gubbins dropped
```

---

## Validation status — what has and has not actually run

Measured on 10 real *B. pseudomallei* assemblies in the development sandbox, which has **no container
engine** and therefore no snippy/gubbins/iqtree on PATH. From `pipeline_info/execution_trace_*.txt`:

| Process | Status | Notes |
|---|---|---|
| `INFILE_HANDLING_UNIX` | COMPLETED (10/10) | |
| `MASH_SKETCH_BATCH` | COMPLETED | |
| `MASH_PASTE` | COMPLETED | |
| `MASH_TRIANGLE` | COMPLETED | publishes `Clustering/mash_distances.phylip` |
| `CLUSTER_GENOMES` | COMPLETED | publishes `Summaries/clusters.tsv` |
| `SELECT_CLUSTER_REPRESENTATIVE` | COMPLETED (2/2) | emits `representative.fa` |
| `SNIPPY_SCATTER` | **FAILED, exit 127** | `.command.sh: line 20: snippy: command not found` |
| everything downstream | **NEVER RAN** | blocked by the SNIPPY_SCATTER failure |

**What the SNIPPY_SCATTER failure does and does not tell us.** The task got past the module's own
reference guard (which previously aborted with `reference for cluster ... is missing or empty`) and
reached the `snippy` invocation on line 20, so the representative-FASTA wiring is fixed as far as
argument construction goes. It is **not** evidence that snippy runs correctly: exit 127 is the shell
failing to find the binary. `snippy-core`, Gubbins, IQ-TREE, grafting, and `publish_dir_mode='link'`
remain **entirely unvalidated**.

Step 3 above is what closes these gaps, and it needs a machine with Docker.

---

## Step 4 — reproduce your previous results

Every behaviour change is reversible. To confirm the optimizations alone did not alter your output:

```bash
nextflow run . -profile low_spec,docker --recombination_aware_mode true \
    --input ~/Downloads/subset_100 --outdir results_legacy \
    --mash_threshold 0.002 \
    --cluster_split_method order \
    --max_column_missingness 0.0 \
    --iqtree_support true \
    --starting_tree_builder iqtree --iqtree_starting_model MFP
```

Cluster assignments and filtered alignments should match your old run byte-for-byte. Trees will differ
only where the old code emitted star trees (it could not produce a real tree there).

---

## Making the repo usable by others

1. **Push the branch and open a PR** — the commits are self-documenting, each explaining what changed
   and what was measured.
2. **Add a CI smoke test.** GitHub Actions can run the workflow on 3–5 tiny genomes with
   `-profile test,docker` on every push. This catches breakage like the pandas 3.0 crash
   (`AUDIT_REPORT.md` §B) before a user hits it.
3. **State the tested versions in the README**: Nextflow 25.04.6, Java 17, and the pinned tool
   versions in `environment.yml`.
4. **Ship `AUDIT_REPORT.md` and `MIGRATION_NOTES.md` in `docs/`** so users understand which defaults
   changed scientific output and how to restore old behaviour.

A caution worth putting in your README: **existing results may contain star trees.** Any
`.final.treefile` with zero internal nodes carries no phylogenetic information but looks like a normal
result. Anyone re-using published output from the old code should check with the loop in Step 3.
