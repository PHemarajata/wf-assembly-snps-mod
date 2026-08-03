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
