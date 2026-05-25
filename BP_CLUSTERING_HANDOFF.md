# BP Clustering Project Handoff Note

## Executive summary

This project is a **BP population genomics / clustering workflow** focused on integrating approximately **3,000 contextual public assemblies** with approximately **326 study isolates** (260 patient and 66 environmental) from one province in Thailand. identifying related groups, building cluster-level phylogenies, handling genetic changes, and preparing outputs that can support epidemiologic interpretation.

The main analytical challenge has been that **BP population structure is complex**, and clustering behavior has been sensitive to method choice, sample naming, single-linkage chaining, singleton burden, and downstream phylogeny constraints. Several strategies were explored:

1. PopPUNK / PopPIPE-based population clustering.
2. Mash-distance threshold clustering in a custom Nextflow workflow.
3. Threshold selection using empirical Mash distance distributions and component sweeps.
4. Per-cluster phylogenetic reconstruction using SKA, Gubbins, and IQ-TREE.
5. Recombination-aware workflow modifications.
6. Robust error handling for singleton-heavy, duplicate-heavy, or low-information clusters.
7. Future tree-grafting / integration of cluster-level trees onto a broader backbone tree.

The project is currently best understood as being in a **workflow stabilization and interpretation-preparation stage**. Core issues around contig-name hygiene, PopPUNK/PopPIPE input consistency, Mash threshold tuning, Gubbins failure modes, IQ-TREE rerun behavior, and Nextflow module robustness have been substantially worked through. The remaining handoff need is to consolidate final outputs, determine which clustering framework is preferred for reporting, and document any CDC HQ feedback.

---

## Source material reviewed

This handoff was drafted from the uploaded transcript PDFs:

- `BP clustering - PopPIPE outputs overview.pdf`
- `BP clustering - Mash threshold determination.pdf`
- `BP clustering - Adjust clustering parameters.pdf`
- `BP clustering - Sanitize contig names.pdf`
- `BP clustering - Nextflow module error fix.pdf`
- `RE: Genomic analysis consultation for BP.pdf`
- `Re: CDC:APHL BP Genomics Zoom Meeting.pdf`

---

## 1. Project background

The project concerns **population genomic clustering of Burkholderia pseudomallei**. The working dataset described in the transcripts included:

- Approximately **3,000 contextual B. pseudomallei assemblies** from public repositories.
- Approximately **326 study isolates** from one province in Thailand (**260 patient** and **66 environmental**).
- Study isolates included both **human and environmental** samples.
- Assemblies were organized under a working directory such as `bp-megamix/assemblies_sanitized/`.
- PopPUNK / PopPIPE outputs were already available at one stage, including:
  - `combined_rfiles.txt` mapping `Taxon` to FASTA path.
  - `combined_clusters.clean.csv` mapping `Taxon` to PopPUNK cluster.
  - PopPUNK HDF5 database such as `bp_all_assignments.h5`.
  - A full/global tree file such as `full_tree.nwk`.
  - PopPIPE strain folders under paths resembling `output/strains/<cluster>/...`.

The primary goal was to build a defensible downstream analysis workflow from population clustering outputs toward recombination-aware cluster trees and ultimately interpretable population-genomic structure.

---

## 2. Core analytical question

The central question was not simply "which samples cluster together?" but rather:

> How should BP assemblies be grouped into analytically useful clusters so that each cluster is large enough and coherent enough for downstream phylogeny, but not so broad that single-linkage chaining, background population structure, or distant outliers collapse unrelated genomes into a misleading mega-cluster?

Operationally, the project needed to answer:

1. Which clustering strategy should be used: PopPUNK/PopPIPE, Mash-threshold clustering, or both?
2. How should thresholds be selected empirically rather than guessed?
3. How should singletons be handled?
4. How should clusters that are too small, too large, duplicate-heavy, or low-information be treated?
5. How should recombination-aware trees be produced per cluster?
6. How should downstream outputs be summarized so another analyst can review and continue the work?

---

## 3. Chronology of approaches and decisions

### 3.1 Initial PopPUNK / PopPIPE direction

The project initially proceeded from **PopPUNK / PopPIPE outputs**. The working assumption was that PopPUNK clusters could provide a biologically meaningful population-structure scaffold, and PopPIPE output folders could be mined for downstream phylogenetic products.

A proposed downstream plan was:

1. Merge `combined_rfiles.txt` with `combined_clusters.clean.csv` on `Taxon`.
2. Create `metadata/cluster_membership.tsv` with one row per isolate.
3. Create `metadata/cluster_sizes.tsv` with cluster sizes and pass/fail flags.
4. Select clusters above a minimum size threshold.
5. For each selected cluster:
   - Generate a manifest: `results/{cluster}/cluster_isolates.tsv`.
   - Build a SKA core SNP alignment.
   - Run Gubbins for recombination detection/masking.
   - Run IQ-TREE on the recombination-filtered alignment.
6. Produce a global `cluster_qc_summary.tsv` or `cluster_phylogeny_summary.csv`.
7. Use those outputs for a later grafting step onto `full_tree.nwk`.

A suggested initial `cluster_min_size` was **6**. The reasoning was that 4 is a bare minimum for a tree, but 6 gives more support and makes Gubbins more worthwhile.

### 3.2 PopPUNK/PopPIPE input consistency issues

A major lesson was that PopPIPE is strict about consistency among:

- PopPUNK rfile: one line per `Taxon` with FASTA path.
- PopPUNK cluster CSV: cluster assignment for the same exact `Taxon` IDs.
- PopPUNK HDF5 database: the database corresponding to the same exact genome set.

A specific failure mode involved Thai reference isolates appearing twice under inconsistent IDs, for example:

- `GCA_963564585`
- `GCA_963564585_1`

This produced mismatches where `combined_clusters.csv` contained one naming universe while `combined_rfiles.txt` contained another. PopPIPE then failed with `KeyError`-type errors because names existed in one input but not the others.

The resulting rule was:

> For PopPIPE, the rfile, clusters CSV, and HDF5 database must describe the same set of Taxon IDs, with exactly matching names and valid FASTA paths.

Best-practice recommendations developed during the work:

- Do not manually merge or edit `.h5` files.
- Use `poppunk_assign --update-db` if new genomes need to become part of the PopPUNK DB.
- Use `--write-references` to include reference isolates in cluster output without duplicating them as queries.
- If Thai references are already inside the PopPUNK reference database, do not also list them as queries under altered IDs.
- Build a single deduplicated combined rfile for PopPIPE.
- Rebuild `combined_clusters.clean.csv` so every `Taxon` appears in the rfile exactly once.

### 3.3 Conservative contig-name sanitization

Because tools including Snippy, BUSCO, Gubbins, QUAST, alignment tools, and phylogeny tools can be sensitive to sequence headers, a conservative contig-name sanitizer was developed.

The explicit goal was to avoid "oversanitizing," because a prior attempt had changed names too aggressively and caused downstream problems.

The conservative sanitization strategy was:

1. Strip leading/trailing whitespace.
2. Replace internal whitespace with a single underscore.
3. Remove only unsafe characters such as parentheses, commas, colons, semicolons, slashes, backslashes, pipes, and quotes.
4. Keep alphanumerics, underscores, dashes, and dots.
5. Preserve contig order and sequences.
6. Ensure uniqueness with minimal suffixing only when needed.
7. Write a mapping file from original to sanitized names.

This became an important quality-control foundation for downstream reproducibility.

### 3.4 Shift toward empirical Mash-threshold clustering

At one point the project moved away from relying only on PopPUNK and explored **simple Mash-threshold clustering**. The motivation was to use a direct, transparent distance threshold to define components for downstream analysis.

The general threshold-selection framework discussed was:

1. Generate stable Mash sketches, ideally using bacterial-genome-appropriate sketch parameters such as `-k 21` and larger sketch sizes.
2. Compute all-vs-all Mash distances.
3. Inspect the distance distribution using histograms or ECDFs.
4. Look for natural valleys between within-population and outlier/distant modes.
5. Sweep thresholds and measure:
   - number of components,
   - singleton fraction,
   - largest component size,
   - component-size distribution,
   - downstream tractability for Gubbins/IQ-TREE.
6. Choose the smallest threshold that captures meaningful close relationships without collapsing everything through single-linkage chaining.

### 3.5 First Mash-distance result and recommended threshold

One uploaded/attached `mash_dist.txt` was analyzed in the prior conversation. The observed distribution was summarized as:

- Strong main bulk around **0.004-0.006**.
- Median approximately **0.00445**.
- 98th percentile approximately **0.00599**.
- 99th percentile approximately **0.00646**.
- A large gap after approximately **0.00705** until approximately **0.0624-0.064**, interpreted as likely outgroup/junk/distant genomes.

The initial recommended Mash threshold was therefore:

- Start with `--mash_threshold 0.0068`.
- If needed, increase slightly to `0.0070`.
- Avoid values above approximately `0.008` because they would not add meaningful close edges until the much more distant ~0.06 group.

Recommended accompanying settings included:

- `--max_cluster_size 300-400` for large components before Gubbins, given a 12-core / 64 GB machine.
- `--merge_singletons false` while tuning.
- `--merge_singletons true` only as a pragmatic fallback when downstream trees were needed despite remaining singletons.
- `--mash_sketch_size 50000-100000` to reduce noise near tight thresholds.

### 3.6 All-singletons and one-mega-component behavior

Two opposite failure modes were encountered:

#### All singletons

When clustering returned all singletons, the interpretation was that the Mash cutoff was too strict. Suggested fixes were:

- Increase `--mash_threshold`.
- Use a threshold ladder rather than guessing.
- Increase `--max_cluster_size` if large clusters were being split too aggressively.
- Use `--merge_singletons true` only if the goal was to push samples through downstream analysis quickly.

#### One giant component

When everything collapsed into one group, the interpretation was that the threshold was too loose and single-linkage chaining was connecting the dataset through weak edges.

Recommended response:

- Tighten threshold in small steps:
  - `0.0060`
  - `0.0056`
  - `0.0052`
- Use `--max_cluster_size 150-200` as a practical guardrail.
- Keep `--merge_singletons false` during tuning so that singleton behavior remains visible.
- Target a largest cluster size around **150-250** for practical Gubbins performance on a 12-core / 64 GB workstation.

This was a key analytical insight: **for this dataset, Mash threshold selection is not just about close relatedness; it is also about controlling graph percolation.**

### 3.7 Later threshold sweep on a new dataset

A later sweep used a valley value:

- `VALLEY=0.005055`

The user ran a command resembling:

```
python sweep_components.py \
  --valley "$VALLEY" \
  --min-threshold 0.002 \
  --anchors 0.0005,0.001,0.005,0.01 \
  --recommend > mash_threshold_sweep.tsv
```

Visible sweep output included threshold/component/singleton-fraction rows such as:

- `t=0.003033`, `k=138`, `singleton_frac=0.047265`
- `t=0.003538`, `k=113`, `singleton_frac=0.041864`
- `t=0.004044`, `k=97`, `singleton_frac=0.035787`
- `t=0.004549`, `k=81`, `singleton_frac=0.02971`

This indicates that the project had progressed beyond ad hoc thresholding and was using a formal threshold sweep to study how components behaved across candidate cutoffs.

### 3.8 Cluster-size distribution after one clustering run

At one point, a clustering summary showed many small clusters:

- Median cluster size: **1**.
- 75th percentile: **2**.
- 90th percentile: **6**.
- 95th percentile: **11**.
- Clusters of size 2-5: **314 clusters**, representing approximately **21.0%** of samples.
- Clusters of size 6-20: **90 clusters**, representing approximately **23.6%** of samples.
- Clusters >50: **10 clusters**, representing approximately **24.9%** of samples.

Interpretation: the clustering result was not simply "wrong," but it had a high singleton/small-cluster burden. This created a downstream tension: tighter thresholds prevent false aggregation but leave many clusters too small for robust recombination-aware phylogeny.

### 3.9 Meaning of `--merge_singletons`

The `--merge_singletons` option was clarified as follows:

- After initial Mash clustering, clusters of size 1 are examined.
- Each singleton is assigned to the nearest eligible non-singleton cluster by Mash distance.
- The merge only occurs if the nearest cluster is within the Mash threshold and does not violate `--max_cluster_size`.
- If no eligible cluster is close enough, the sample remains a singleton.
- The tradeoff is fewer trivial clusters and better downstream throughput, at the cost of potentially broadening within-cluster diversity.

Practical recommendation:

- Use `--merge_singletons false` while determining thresholds.
- Consider `--merge_singletons true` only after threshold behavior is understood and if downstream completeness is more important than preserving singleton status.

---

## 4. Workflow implementation and troubleshooting history

### 4.1 PopPIPE already produced IQ-TREE outputs but not Gubbins in some cases

The transcripts note that PopPIPE had already been run and included IQ-TREE outputs, but did not necessarily run Gubbins where there was no transmission file. This prompted the question of what to do with existing strain folders and whether to retrofit recombination-aware per-cluster trees.

The approach became:

- Reuse existing PopPIPE strain folders where possible.
- Add a Snakemake layer to run or rerun Gubbins and IQ-TREE per qualifying cluster.
- Produce recombination-aware or sanitized-alignment fallback trees.
- Summarize outputs in `cluster_phylogeny_summary.csv`.
- Later graft cluster trees onto `full_tree.nwk`.

### 4.2 Snakemake / Gubbins / IQ-TREE layer

A Snakemake workflow was designed to process qualifying clusters. Major intended outputs included:

- `output/strains/<cluster>/gubbins/gubbins_q.tsv`
- `output/strains/<cluster>/gubbins/gubbins_<cluster>.filtered_polymorphic_sites.fasta`
- `output/strains/<cluster>/iqtree_recomb/<cluster>.recomb.treefile`
- `output/strains/<cluster>/iqtree_recomb/<cluster>.recomb.iqtree`
- `cluster_phylogeny_summary.csv`

Several issues were resolved or partially resolved:

#### Snakemake wildcard formatting issue

A rule contained literal `{cluster}` inside a shell block or comment where Snakemake attempted to expand it incorrectly. The fix was to consistently use:

- `{wildcards.cluster}` for wildcard references.
- `{params.*}` for parameters.

#### IQ-TREE checkpoint guard

IQ-TREE failed on clusters such as 19 and 26 because a previous checkpoint indicated successful completion and IQ-TREE refused to overwrite without `-redo`. The fix was to make the rule idempotent:

- If both the treefile and `.iqtree` log already exist, skip IQ-TREE.
- Otherwise, run IQ-TREE with `-redo`.

This prevents Snakemake from failing on harmless previous-completion state.

#### Gubbins duplicate/low-information failures

Gubbins failed for some clusters, for example cluster `31_72`, with errors consistent with too few sequences remaining after duplicate removal. The solution was a two-layer fallback:

1. Pre-run QC could mark `run_gubbins=0` and skip Gubbins for low-information clusters.
2. If Gubbins was attempted but failed at runtime, the workflow would:
   - log a warning,
   - copy `align_variants_sanitized.aln` to the expected filtered output name,
   - exit successfully so Snakemake could continue.

Interpretation rule:

> IQ-TREE can still build a tree from the sanitized alignment, but downstream summaries must flag whether recombination masking was actually performed or whether the cluster used the sanitized-alignment fallback.

This distinction matters for final interpretation.

---

## 5. Nextflow recombination-aware workflow work

The project also involved a separate or modified Nextflow workflow repository:

- `PHemarajata/wf-assembly-snps-mod`
- Related prior repo: `wf-assembly-snps-final`

A run command included parameters resembling:

```
nextflow run main.nf \
  -profile local_workstation_rtx4070,docker \
  --input /home/phemarajata/Downloads/bp_finalset \
  --outdir /home/phemarajata/Downloads/bp_finalset/results_recombination_aware \
  --recombination_aware_mode \
  --integrate_results \
  --mash_threshold 0.028 \
  --max_cluster_size 50 \
  --merge_singletons \
  --mash_sketch_size 50000 \
  --recombination gubbins \
  --snp_package parsnp \
  --run_gubbins \
  -resume
```

### 5.1 Module compilation error

A `BUILD_BACKBONE_TREE` module failed at parse/compile time. The suspected causes were:

- UTF-8 BOM or hidden characters at the top of the module.
- CRLF line endings.
- DSL2 not being explicitly recognized.
- Malformed first line or paste artifacts.
- Incorrect `include` path pulling the wrong file.

Recommended checks included:

- Inspecting first lines with `nl -ba`.
- Searching for the include line.
- Adding `nextflow.enable.dsl=2` to the top of the module.
- Removing BOM and CRLF line endings.
- Clearing `.nextflow/cache` before rerun.
- Checking Nextflow version.

### 5.2 Rewriting backbone tree module

The project then moved toward completely rewriting the failing `BUILD_BACKBONE_TREE` module using a comparable approach. The goal was to produce a robust backbone tree component that fits into the recombination-aware workflow and does not derail the pipeline due to brittle module code.

### 5.3 IQTREE_ASC process fixes

An `IQTREE_ASC` process produced multiple errors. Key lessons and fixes:

#### Invalid `--asc-corr` flag

IQ-TREE 2 expects ascertainment correction through the model string, for example:

- `GTR+ASC`

not through an invalid `--asc-corr` flag.

#### `args: unbound variable`

The process-generated bash script had an `args` variable that was referenced unsafely. The replacement process handled optional flags via `task.ext.args` with a safe default.

#### `bad substitution`

A later error arose from shell incompatibilities around command substitution inside a heredoc or use of bash-specific constructs in a container shell. The fix was to make the process more POSIX-friendly:

- Avoid bash arrays and `[[ ... ]]`.
- Detect `iqtree2` vs `iqtree` safely.
- Write `versions.yml` using `printf` after computing the version string.
- Create minimal expected outputs if IQ-TREE cannot run.

The desired behavior was to preserve workflow continuity while clearly flagging failures.

---

## 6. Main hypotheses and how they evolved

### Hypothesis 1: PopPUNK clusters are the best population-genomic units

This was plausible because PopPUNK has been used in BP studies and provides model-based population clustering. A paper using PopPUNK for BP motivated revisiting this path and attempting to reproduce comparable analytical steps.

Status: still viable, especially for population-structure framing. However, PopPUNK/PopPIPE requires strict ID consistency and careful handling of reference/query roles.

### Hypothesis 2: Simple Mash-threshold components are more transparent and operationally controllable

This approach became attractive because it allows direct threshold tuning and component-size control. It is easier to reason about when trying to manage all-singleton vs mega-component behavior.

Status: viable for operational clustering and batching into downstream phylogeny, but sensitive to threshold and single-linkage chaining.

### Hypothesis 3: There is a natural Mash threshold around the empirical distance valley

Early Mash-distance analysis suggested a threshold around `0.0068-0.0070`, based on the p99 value and a large gap before distant genomes. Later runs showed that this could still collapse to one large component depending on the dataset and graph connectivity.

Status: partially supported, but threshold must be dataset-specific and validated with component sweeps. For some data states, lower thresholds around `0.0060`, `0.0056`, `0.0052`, or around a valley near `0.005055` may be more appropriate.

### Hypothesis 4: Many clusters are too small or low-information for full recombination-aware phylogeny

Cluster-size summaries and Gubbins runtime failures support this. Many clusters are singletons or very small, and some clusters collapse after duplicate removal.

Status: strongly supported. The workflow should explicitly classify clusters by whether Gubbins ran, was skipped, or failed and fell back to sanitized alignment.

### Hypothesis 5: The best final workflow may combine PopPUNK population context with cluster-specific SNP/recombination phylogenies

This became the most sophisticated working model: use PopPUNK/PopPIPE or Mash components to define population units, generate cluster-level recombination-aware trees where possible, and then integrate or graft those trees onto a broader population scaffold.

Status: promising but not yet documented as complete in the supplied transcripts.

---

## 7. CDC HQ / APHL technical consultation

The uploaded consultation records document both the request for CDC technical input and the follow-up action plan after the CDC/APHL Zoom meeting.

### 7.1 Consultation participants and timing

The primary CDC technical contact was **Christopher Gulvik, CDC/NCEZID/DHCPP/BSPB**. The consultation request was sent by Peera Hemarajata/APHL on **January 12, 2026**, with CDC response on **January 14, 2026**. The Zoom meeting was scheduled for **Wednesday, January 21, 2026, 7:00-8:00 PM Thailand time**, corresponding to morning Atlanta time.

People copied or involved across the email chain included APHL Thailand and CDC/GHC/DGHP colleagues, including Kornthara Kawang, Saithip Bhengsri, Rachel Suzanne Beard, Pongpun Sawatwong, and Famui Mueanpai. Taylor Paisie was copied on the initial consultation request.

### 7.2 Questions posed to CDC

The consultation request framed the project as a technical review of the BP sequencing dataset and asked whether the team's approaches for clustering and relatedness were appropriate for BP and would avoid misleading interpretation as more analysis was done locally.

The dataset was described at that time as approximately **326 isolates**, including **260 patient isolates** and **66 environmental isolates**. The analysis used:

- `wf-paired-end-illumina-assembly`
- TheiaProk on Terra.bio
- ARDaP for AMR variant interpretation
- Mash for initial clustering
- PopPUNK for clustering
- PopPIPE for generating clusters/subclusters
- Gubbins for per-cluster recombination analysis
- ML trees for clusters
- an attempted subtree-grafting strategy

The team specifically noted that Mash-based groups appeared to reflect **background population structure rather than micro-transmission clusters**, with patient isolates distributed across broader diversity and environmental isolates grouping within different major clades. No direct patient-environment isolate pairings had been observed at that time.

### 7.3 CDC written feedback before the Zoom call

Chris Gulvik's written feedback provided several important interpretive guardrails:

1. **The dataset size was strong.** He noted that 300+ BP genomes was very impressive.
2. **Failure to find a direct patient-environment match was not unexpected.** He emphasized that even with a large sample size, finding a direct "match" between patient and environment is very rare for BP. This should not be interpreted as failure of effort; it reflects the biology and abundance of the organism in endemic environments.
3. **No strong environmental resistance signal was reassuring.** He described the absence of strong resistance signals in environmental BP genomes as positive for safety/public health interpretation.
4. **Theiagen and ARDaP were viewed as appropriate solutions.** He affirmed the tools the team had found and used.
5. **PopPUNK was endorsed as a good first-pass comparison tool.** He specifically described PopPUNK as useful and scalable because it is k-mer based, and he noted that its clustering feature is helpful.
6. **PopPUNK was not considered ideal for fine-scale BP resolution.** Because PopPUNK is k-mer based, Chris cautioned that it lacks the fine-scale resolution needed for outbreak investigations or detailed epidemiologic studies in BP. When genomes are very similar, k-mer differences may be dominated by sequencing errors or assembly artifacts rather than true evolutionary signal.
7. **SNP-based approaches were recommended for finer resolution.** He named Snippy, Parsnp, and `wf-assembly-snps` as more appropriate for fine-scale work and offered an attached `wf-assembly-snps.bash` script as a possible simple solution.
8. **Small test runs were recommended before scaling.** He suggested first verifying the workflow with a built-in test, then trying approximately 10 samples, and only then running all 300+ samples.
9. **He offered to run analysis if needed.** He indicated he would be willing to analyze the FASTA files if the team needed a fast answer.

### 7.3a CDC-provided `wf-assembly-snps.bash` script

The uploaded `wf-assembly-snps.bash` script is a standalone Bash workflow intended to run SNP-based comparative analysis on assembled genomes. It reinforces CDC's recommendation to use a SNP-based method for finer relatedness assessment rather than relying only on PopPUNK or Mash clustering.

Key features of the script:

- Accepts an input directory of assembly FASTA files ending in `fa`, `fas`, `fsa`, `fna`, or `fasta`, including gzipped files.
- Requires at least **3 genomes** for batch analysis.
- Includes two test modes:
  - `test`: four small PhiX genomes.
  - `test.full`: four BP genomes, including K96243, 1026b, ATS2021, and another contextual genome.
- Requires Conda and Bash 4+.
- Creates or reuses Conda environments for specific tools.
- Copies input files into an output working directory, decompresses gzipped inputs, strips FASTA extensions, performs conservative filename sanitization, and removes very small assemblies below a minimum size threshold of `5k`.
- Selects the largest input genome as the Parsnp reference.
- Runs Parsnp 2.1.3 to create a core genome alignment with `--skip-phylogeny`.
- Falls back to a modified Parsnp command with `--no-partition --curated` if the first Parsnp attempt fails.
- Uses HarvestTools to extract `parsnp.aln` from the Parsnp GGR output.
- Strips a `.ref` suffix from the alignment and Newick tree labels only if it occurs exactly once.
- Uses SNP-sites 2.5.1 to extract core SNPs and generate a VCF.
- Uses SNP-dists 0.8.2 to produce both SNP distance matrix and pairwise SNP distance outputs.
- Builds a tree using either:
  - IQ-TREE2 2.4.0 with `GTR+G` by default, or
  - FastTree 2.1.11 if `TREE_METHOD=fasttree` is set.
- Includes optional recombination modules controlled by a `RECOMBINATION` variable:
  - `skip` by default,
  - `gubbins`,
  - `clonalframeml`, or
  - `both`.

Important implementation notes:

- The default script setting is `RECOMBINATION="skip"`, so recombination analysis is not performed unless the script is edited or externally modified.
- The Gubbins branch uses Gubbins 3.1.4 and starts from the tree generated by the SNP workflow.
- The ClonalFrameML branch uses ClonalFrameML 1.12.
- The script is intentionally simple and useful for testing or demonstration, but it is not yet a fully production-hardened workflow for hundreds or thousands of genomes.
- It chooses the largest genome as reference, which is pragmatic but may not always be biologically optimal for final reporting.
- It creates Conda environments dynamically, which is convenient but may reduce reproducibility unless environment versions and logs are archived.
- It removes temporary copied input and reference folders at the end, so production runs should preserve external input metadata and logs separately.

Interpretation for this project:

> Chris's script supports the direction of using Parsnp/SNP-sites/SNP-dists/IQ-TREE, optionally with Gubbins or ClonalFrameML, as a more appropriate fine-scale analysis path for BP than PopPUNK/Mash alone. However, the script should be treated as a technical template or validation path, not as the sole final production pipeline without additional QC, metadata tracking, cluster batching, recombination status reporting, and reproducibility controls.

### 7.4 Post-call action plan from January 22, 2026

After the January 21 CDC/APHL Zoom call, the internal follow-up summary stated that the discussion went well and that the team's approach did not appear fundamentally off track. The action plan after discussion with Chris was:

1. **Confirm discrepant phenotypic AST vs WGS results using BD Phoenix.**
2. **Determine doxycycline MICs** for isolates where ARDaP predicted potentially reduced susceptibility. This was noted as an item of interest from David Su during a prior HQ call.
3. **Continue generating a tree with contextual data** to assess whether geographic clustering can be identified. This was considered particularly relevant given recent North American cases without travel history or known exposure to products that may harbor BP, such as the aromatherapy-oil-associated outbreak.
4. **Assess how different the Thailand isolates are from isolates in the continental United States and elsewhere.**
5. **Refine the clustering and SNP-tree algorithm** based on technical tips from Chris, with the goal of improving tree quality.
6. **Do not pursue further MLST analysis** for isolates that did not type in the initial analysis, provided the team confirms they are truly BP and not another organism.

### 7.5 Interpretation of CDC consultation for this project

The consultation supports the project's current hybrid direction:

- PopPUNK is useful as a broad, scalable, first-pass population-structure tool.
- PopPUNK/Mash groupings should not be overinterpreted as fine-scale transmission clusters.
- Fine-scale relatedness questions should rely on SNP-based workflows such as Snippy, Parsnp, or `wf-assembly-snps`, with recombination-aware analysis where appropriate.
- Direct patient-environment genome matches should not be expected, even with hundreds of genomes.
- Contextual phylogeny is important, especially for evaluating geographic structure and relationships to non-Thailand or North American isolates.
- MLST is not essential for isolates that fail initial typing, as long as species identity is confirmed.

### 7.6 Remaining CDC-related gaps

The `wf-assembly-snps.bash` script has now been reviewed and summarized in this handoff. Remaining CDC-related gaps are limited to any full verbal details from the January 21, 2026 Zoom call that were not captured in the emails or script, especially specific comments on parameter choices, reference selection, recombination handling, or how CDC would prioritize the analysis for reporting.

---

## 8. Current project status

Based on the transcripts, the project currently stands here:

### Completed or substantially addressed

- Established project context: BP population genomics with ~3,000 contextual and ~326 study isolates (260 patient, 66 environmental).
- Identified key PopPUNK/PopPIPE inputs and the need for strict `Taxon` ID consistency.
- Developed conservative contig-name sanitization logic.
- Designed a PopPUNK/PopPIPE downstream workflow for cluster membership, cluster sizes, per-cluster phylogeny, and QC summaries.
- Evaluated Mash threshold selection using distance distributions and component sweeps.
- Identified failure modes: all singletons, one mega-component, high singleton burden, Gubbins duplicate failures, IQ-TREE checkpoint failures.
- Developed robust Snakemake fixes for Gubbins and IQ-TREE reruns.
- Developed Nextflow troubleshooting and process rewrites for recombination-aware pipeline components.

### Partially complete / uncertain

- Final chosen clustering method: PopPUNK/PopPIPE, Mash threshold components, or hybrid.
- Final Mash threshold for the production dataset.
- Final handling of singletons.
- Final list of qualifying clusters.
- Final `cluster_phylogeny_summary.csv` contents.
- Whether tree grafting onto `full_tree.nwk` was implemented.
- Whether the Nextflow recombination-aware workflow now completes end-to-end.
- CDC HQ consultation findings.

### Not documented in supplied transcripts

- Final result tables.
- Final figures.
- Final cluster confidence interpretation.
- Full verbal details of the January 21, 2026 CDC/APHL Zoom consultation.
- Epidemiologic interpretation of specific human/environmental clusters.
- Any outbreak/actionability conclusion.

---

## 9. Practical recommendations for whoever takes over

### 9.1 First, freeze the sample universe

Before doing more analysis, define one canonical sample table with:

- `Taxon`
- FASTA path
- study/contextual flag
- human/environmental flag
- province/location metadata
- collection date or year
- source repository/accession
- QC metrics
- PopPUNK cluster, if available
- Mash component, if available

Do not allow the same biological isolate to appear under multiple Taxon IDs.

### 9.2 Decide whether PopPUNK, Mash, or both are the reporting scaffold

A reasonable final structure would be:

- Use **PopPUNK** for broad population context and comparability to published BP approaches.
- Use **Mash threshold components** for operational batching and sensitivity analysis.
- Use **cluster-level SNP/recombination phylogenies** for within-cluster resolution.

Do not present Mash threshold components as epidemiologic clusters without additional support.

### 9.3 Use threshold sweeps, not a single guessed cutoff

For Mash-based clustering, rerun or preserve a threshold sweep table containing:

- threshold,
- number of components,
- singleton fraction,
- largest component size,
- number of clusters above sizes 3, 6, 20, 50,
- percent of samples in singletons,
- percent of samples in clusters suitable for phylogeny,
- number of clusters exceeding practical Gubbins size.

Choose the threshold that balances biological resolution and tractability.

### 9.4 Preserve singleton status during threshold selection

Keep `--merge_singletons false` while evaluating candidate thresholds. Only enable singleton merging once the threshold is selected and the consequences are reviewed.

### 9.5 Track recombination handling explicitly

For every cluster tree, record one of the following statuses:

- `gubbins_completed`
- `gubbins_skipped_low_information`
- `gubbins_failed_fallback_to_sanitized_alignment`
- `iqtree_completed`
- `iqtree_skipped_existing_outputs`
- `iqtree_failed_minimal_tree_output`

This is essential for transparent downstream interpretation.

### 9.6 Do not overinterpret cluster membership

Given BP population structure and the observed threshold sensitivity, cluster membership should be described as **genomic grouping for analysis**, not proof of recent transmission or direct epidemiologic linkage.

Recommended wording:

> These clusters represent genomic groupings under the specified PopPUNK or Mash-threshold approach. They should be interpreted with the accompanying temporal, geographic, source, and phylogenetic evidence, and not as direct evidence of transmission by themselves.

---

## 10. Suggested final deliverables

A complete handoff package should include the following files or equivalents:

### Metadata and clustering

- `canonical_sample_metadata.tsv`
- `combined_rfiles.txt`
- `combined_clusters.clean.csv`
- `cluster_membership.tsv`
- `cluster_sizes.tsv`
- `mash_threshold_sweep.tsv`
- `clusters.tsv`
- `cluster_summary.txt`
- `singletons.txt`

### Phylogeny and recombination

- per-cluster FASTA manifests
- per-cluster sanitized alignments
- per-cluster Gubbins outputs
- per-cluster IQ-TREE outputs
- `cluster_phylogeny_summary.csv`
- `full_tree.nwk`
- grafted or integrated tree, if completed

### Workflow/code

- contig-name sanitization script
- Snakemake workflow for per-cluster Gubbins/IQ-TREE
- Nextflow workflow or modified modules
- config files used for final runs
- exact commands used for final production run

### Interpretation

- cluster-confidence table
- methods text
- limitations text
- CDC HQ consultation notes
- final summary memo or slides

---

## 11. Proposed cluster confidence framework

For downstream communication, classify clusters using evidence tiers.

### Tier 1: High-confidence genomic analysis cluster

Criteria:

- Stable under reasonable thresholds.
- Not created only by singleton merging.
- Adequate cluster size.
- Gubbins completed or recombination status is clear.
- IQ-TREE completed successfully.
- Cluster is coherent on phylogeny.
- Supported by temporal/geographic/source metadata, if available.

### Tier 2: Moderate-confidence genomic analysis cluster

Criteria:

- Reasonable genomic grouping, but some limitations.
- May be sensitive to nearby threshold choices.
- Gubbins may have been skipped or replaced by sanitized-alignment fallback.
- Epidemiologic metadata incomplete or mixed.

### Tier 3: Low-confidence / exploratory grouping

Criteria:

- Singleton-merged or threshold-sensitive grouping.
- Very small cluster.
- Duplicate-heavy or low-information alignment.
- No clear metadata support.
- Tree is trivial, minimal, or based on fallback behavior.

### Tier 4: Not interpretable as a cluster

Criteria:

- Singleton.
- Outlier/distant sample.
- Failed QC.
- No stable placement.
- No meaningful downstream phylogeny.

---

## 12. Recommended language for reporting

### General methods language

> We evaluated BP population structure using both model-based clustering outputs and Mash-distance component clustering. Because cluster assignments can be sensitive to threshold choice and sample naming consistency, clustering outputs were treated as analytical groupings and interpreted alongside phylogenetic, recombination, temporal, geographic, and source metadata where available.

### Mash threshold language

> Mash threshold candidates were evaluated using empirical distance distributions and threshold sweeps. We assessed the number of components, singleton fraction, largest component size, and downstream tractability for recombination-aware phylogenetic reconstruction. Thresholds that produced either all singletons or a single giant component were considered unsuitable for final interpretation.

### Gubbins/IQ-TREE language

> For qualifying clusters, recombination-aware phylogenetic analysis was attempted using Gubbins followed by IQ-TREE. Clusters with insufficient information, excessive duplicate sequences, or Gubbins runtime failure were retained using sanitized-alignment fallback outputs and flagged accordingly in downstream summaries.

### Cautionary interpretation language

> Genomic clustering alone does not establish recent transmission or direct epidemiologic linkage. Clusters should be interpreted as genomic groupings under the specified analytical approach, with confidence depending on threshold stability, phylogenetic coherence, recombination handling, and supporting epidemiologic metadata.

---

## 13. Immediate next steps

1. Confirm the canonical sample universe and remove duplicate biological isolates with inconsistent IDs.
2. Rebuild or verify the matched PopPUNK/PopPIPE trio: rfile, clusters CSV, and HDF5 database.
3. Preserve the final Mash threshold sweep output and choose a production threshold based on component behavior.
4. Decide whether singleton merging will be used in the final production run.
5. Generate a final `cluster_membership.tsv` and `cluster_sizes.tsv`.
6. Run or finalize per-cluster Gubbins/IQ-TREE workflow.
7. Generate `cluster_phylogeny_summary.csv` with explicit recombination/fallback status.
8. Decide whether to implement tree grafting onto `full_tree.nwk`.
9. Incorporate any remaining details from the January 21, 2026 CDC/APHL Zoom call, especially specific SNP-tree tips from Chris Gulvik that were not captured in email or the `wf-assembly-snps.bash` script.
10. Prepare a final technical memo with methods, limitations, and cluster-confidence interpretation.
