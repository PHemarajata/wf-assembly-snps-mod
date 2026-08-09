# Handoff: tip annotation, phylogeography and clock prep

**Branch:** `claude/workflow-metadata-mapping-2c2f98` · **PR:** [#4](https://github.com/PHemarajata/wf-assembly-snps-mod/pull/4), open against `main`
**Repo:** https://github.com/PHemarajata/wf-assembly-snps-mod
**Everything below is pushed.** The branch already carries a merge of `main` (PR #3, the
`alignment_length` rename), so it is current.

```bash
git fetch origin && git checkout claude/workflow-metadata-mapping-2c2f98
```

---

## 1. Why this work exists

The user has two things in flight:

1. **A completed analysis by a colleague** — a 2,773-tip grafted *B. pseudomallei* phylogeny
   produced by this workflow on 2026-08-08, plus an interpretation report, figures and cluster
   summaries. Delivered as a folder, not a run directory.
2. **Their own run of the same workflow**, in progress on a different collection.

The goal is **phylogeographic clustering** and, later, **molecular clock / dated phylogenetics**.
The immediate ask was to get corrected metadata mapped onto every tip so those analyses have
something trustworthy to consume, and to write the procedure down so it can be repeated on the
user's own results.

The user edited the metadata sheet themselves before this work started: the `Country_Final` column
for travel-associated isolates (`USA ex ...`) was rewritten to the **country of acquisition**, so US
isolates are attributed to where the infection was contracted rather than where it was diagnosed.
That edit is the reason the remap was needed.

---

## 2. What the PR adds

Three files, 1,109 lines, **no changes to the pipeline**. Both scripts are post-processing that read
the grafted tree plus `Summaries/`.

| file | what it does |
|---|---|
| `bin/build_tip_annotation.py` | joins metadata to tip labels; writes `tip_annotation.tsv`, TreeTime + LSD2 date files, a mugration/DTA traits file, a country change log, QC |
| `bin/clock_signal_check.py` | root-to-tip regression vs collection date, whole-tree and per-cluster |
| `docs/TIP_ANNOTATION.md` | method, column semantics, and the caveats that travel with the table — **read this first** |

Both use only `pandas` + `numpy`. Newick parsing is hand-written and dependency-free (no biopython),
unit-tested against quoted labels, internal-node labels, support values, `[&...]` comments,
scientific-notation branch lengths and missing branch lengths.

---

## 3. Where the data is

The colleague's results were on **macOS OneDrive**; on the Linux laptop these paths will differ.
Substitute a local `$BP` and the commands in `docs/TIP_ANNOTATION.md` work unchanged.

```
<BP>/Yuyi_final_results/
├── megamix_bestshot_cleaned_dropGCF_on_Fdups_APPENDED_acqcountry.tsv   # 2,804 rows, corrected countries
├── Summaries/                        # clusters.tsv, cluster_phylogeny_summary.csv,
│                                     # threshold_sweep.tsv, chosen_threshold.txt (0.005455)
├── BP_Global_interpret/              # PHYLOGENOMIC_REPORT.md, annotated_tips.csv (2,773 rows), figures
├── BP_Global_interpret_260808/       # same figures, fuller set
└── Tip_annotation_20260809/          # OUTPUT of this work — see §4
```

**The grafted tree itself is not in that folder.** `global_grafted.treefile`, `grafting_log.txt`,
`grafting_report.txt` and the colleague's `01_`–`07_` interpretation scripts were never copied out of
the run directory and are not on the user's machine. That is why `build_tip_annotation.py` has a
`--tips-from-csv FILE:COLUMN` mode — the tip labels were taken from `annotated_tips.csv`. **For the
user's own run, pass `--tree`; it is the correct input.**

---

## 4. What was produced

`Yuyi_final_results/Tip_annotation_20260809/` — regenerate with the two commands in
`docs/TIP_ANNOTATION.md` §1.

```
tip_annotation.tsv            2,773 rows x 43 cols — the master table
tip_dates_treetime.csv        2,416 dated tips, midpoint + bounds
tip_dates_lsd2.txt            same, b(lo,hi) for imprecise dates
traits_geography.tsv          geo_state / geo_state_pooled / geo_region / thai_province
country_change_log.tsv        17 travel reattributions + 2 ambiguous countries
clock_signal_by_cluster.tsv   per-cluster root-to-tip regression
tip_annotation_QC.txt         all counts
README.md                     folder-local summary
```

**Result: 2,773 / 2,773 tips matched to metadata. 0 unmatched, 0 duplicate labels.**

---

## 5. Settled findings — do not re-derive these

Each was measured, not assumed. Re-running the two scripts reproduces all of them.

**The join works and is not fragile.** Tip labels are not `sample_id` — they carry assembly versions,
assembler suffixes, FASTA extensions and appended free-text geography, and `GCF_`/`GCA_` are used
interchangeably between the tree and the metadata sheet. The normaliser in `build_tip_annotation.py`
reproduces the colleague's hand-reconciled key for **all 2,773 tips** without the four manual
identifier fixes that pass needed. Unmatched tips abort by default, deliberately.

**Do not date the grafted tree.** Root-to-tip vs collection date over all 2,416 dated tips:
**r = −0.005, R² = 0.0000**. That is the expected consequence of grafting per-cluster trees each
inferred from its own alignment, not evidence about the clock.

**Per-cluster regression is valid** — one alignment per cluster, and the shared backbone path cancels
out of both correlation and slope. 13 of 56 eligible clusters have r ≥ 0.3.

**The rate, not the correlation, is the discriminator.** `IQTREE_ASC` infers from Gubbins'
`filtered_polymorphic_sites.fasta` under `GTR+ASC` with no genome-wide `-fconst`, so a branch length
is substitutions per *variable* site. Converting with `n_filtered_polymorphic_sites` is
**equivalent to what `-fconst` would have done at inference time** — so the reported rates are already
on the corrected scale, and re-inferring with `-fconst` will not change them.

| cluster | r | variable sites | subs/genome/year |
|---|---|---|---|
| cluster_25 | 0.32 | 2,631 | **1.6** |
| cluster_11 | 0.50 | 2,026 | **4.7** |
| cluster_19 | 0.52 | 29,887 | 341 |
| cluster_46 | **0.64** | 106,459 | 1,364 |

Expectation is ~1–10 SNPs/genome/year. `cluster_46` has the best correlation in the set and the least
believable rate: a cluster holding 10⁵ core variable sites across 22 genomes is not one population,
and its apparent temporal signal is lineage membership rather than time. **Choose the dating unit by
depth, not by `cluster_id`.**

**`--iqtree_asc_fallback fconst` does not fix the scale.** `asc_preflight.py` counts only the constant
columns left *inside* the SNP alignment, not the ~6.8 Mb of invariant genome.

**IQ-TREE 3 would not have helped, and the workflow does not use it.** Everything is pinned to
IQ-TREE 2.2.6. IQ-TREE 3 has the same `+ASC` / `-fconst` semantics. Separately, `HANDOFF_TO_CLAUDE_CODE.md`
already settled that no phylogenetics tool here has usable GPU support — do not re-litigate either.
Useful corollary: **IQ-TREE 2.2.6 already embeds LSD2 (`--date`)**, so `tip_dates_lsd2.txt` works
against the version already installed.

**The workflow never roots the tree.** `graft_trees.py` manipulates the root structurally but makes no
biological rooting choice. Rooting is a free post-hoc decision — which matters, because the
colleague's midpoint rooting landed on `GCF_000756925_1`'s 1.44 subs/site terminal branch. That genome
sits in `cluster_73`, the only `Tier3_low_confidence_exploratory` cluster in the run (3 taxa, 27
filtered sites), and the report itself flags it as implausible for *B. pseudomallei sensu stricto*.
The entire "Australia basal → SE Asia derived" narrative rests on that root. **Use an explicit
outgroup.**

**LG1–LG9 are not a typing scheme.** No ST/MLST/cgMLST exists in the metadata. The lineages were
defined topologically by the colleague's `03_lineages.py`: midpoint-root, recursively split the eleven
deepest nodes with a floor of 40 genomes, number basal→derived. The geographic names were attached
afterwards from composition. They are a partition of one particular tree, will not map onto published
population-structure names, and will not carry over to a new run. They are carried in
`tip_annotation.tsv` as `prior_lineage` / `prior_lineage_label`, namespaced for exactly that reason.

**Two data problems in the colleague's run that the report does not resolve:**

- **30 genomes have no tip.** `clusters.tsv` holds 2,803 genomes in 87 clusters, but only the 63
  clusters with ≥3 members were built and grafted. Those 30 cannot appear in any downstream analysis.
- **7 genomes are mislabelled as in-house.** `SRR31608433`–`SRR31608440` are CDC public SRA
  submissions under `PRJNA908850` labelled "In-house prospective survey / Clinical" — inconsistently,
  since `SRR31608436` from the same BioProject and submitter was labelled public. The derived `source`
  column in `tip_annotation.tsv` corrects this; `prior_source` keeps the old call for comparison.

**`cluster_id` is not portable between runs.** `AUDIT_REPORT.md` §G.4: removing 0.3% of genomes
reorganised most clusters and changed the auto-selected threshold. Always report the threshold and the
collection alongside any cluster-based result.

---

## 6. Open work, in priority order

### A. Publish `backbone_alignment.fa` — one line, and it is losing data right now

`BUILD_BACKBONE_TREE` emits `backbone_alignment.fa` but `publishDir` only matches
`backbone.treefile` ([modules/local/build_backbone_tree/main.nf:15](modules/local/build_backbone_tree/main.nf)).
The alignment therefore exists **only in `work/`**. `RUNNING.md` tells users to delete `work/` between
runs — at which point the backbone can never be rebuilt.

This matters because the backbone is `parsnp --use-fasttree` (`-nt -gtr`, no support values), and it
is what produces the deep splits that any lineage definition would rest on. Publishing the alignment
is the cheap route to rebuilding it properly with IQ-TREE afterwards.

Change the pattern to match both files. Verify on a 10-genome smoke run.

### B. Mash-based divergent-genome screen

The `GCF_000756925` outlier was detectable from the Mash triangle **before any tree existed**, and
catching it there would have saved the root problem entirely. Write a small screen over
`Summaries/` — reuse `bin/mash_matrix_io.py` — that flags genomes whose minimum distance to any other
genome exceeds a threshold, i.e. candidates for sister species, contamination or chimeric assemblies.
Needs no pipeline re-run.

### C. Decide `--iqtree_support` on the live run

`iqtree_support = false` is the default ([conf/params.config:159](conf/params.config)); `README.md:418`
says turn it on for anything publishable. The colleague's run did not, so **no clade in that tree has
a confidence value** — a stated limitation of the report.

Support runs inside `IQTREE_ASC`, so adding it later means re-running that stage plus grafting for
every cluster. Changing the param *should* leave snippy and Gubbins cached under `-resume` — that is
the expensive 80% — but **verify on a smoke run before trusting it at full scale.** Heed the
`RUNNING.md:126` warning: a mid-flight container-option change once re-ran all 2,794 alignments and
filled the disk to 99%. Batch config changes; do not drip them in.

### D. Optional: an ML backbone method

`backbone_method` accepts only `parsnp` | `fasttree` ([conf/params.config:249](conf/params.config)).
Adding an IQ-TREE option would give the deep splits a model and support values. Larger change than A;
do A first regardless, since A is what makes D possible retroactively.

### E. Update `docs/TIP_ANNOTATION.md` §7 if the dating route changes

The doc currently recommends choosing the dating unit by depth and, if genome-scale branch lengths are
wanted directly from IQ-TREE, masking the `KEEP_INVARIANT_ATCG` full-length alignment with
`bin/mask_recombination.py` and inferring on that. That masking step is **not** wired into
`recombination_aware_snps.nf` — it exists only in `assembly_snps.nf` — so it currently has to be run
by hand. Wiring it in is a reasonable future change.

---

## 7. Running the scripts

```bash
python3 bin/build_tip_annotation.py \
  --tree       results/Graft/global_grafted.treefile \
  --metadata   metadata.tsv \
  --clusters   results/Summaries/clusters.tsv \
  --cluster-qc results/Summaries/cluster_phylogeny_summary.csv \
  --outdir     results/Tip_annotation
```

```bash
python3 bin/clock_signal_check.py --tree results/Graft/global_grafted.treefile --annotation results/Tip_annotation/tip_annotation.tsv --out results/Tip_annotation/clock_signal_by_cluster.tsv
```

Things that will bite on a new collection:

- **`REGION` and `COUNTRY_ALIASES` in `build_tip_annotation.py` are hard-coded.** An unmapped country
  is a hard error, not a blank — by design. Add new countries to the dict.
- **`--min-state-count`** (default 5) controls how aggressively rare countries fold into their region
  for `geo_state_pooled`. 40 countries → 25 states on this data.
- **`--carry FILE:COL,COL`** brings columns forward from a previous annotation, prefixed `prior_`, so
  a carried `Source` can never overwrite the derived `source`.

---

## 8. What is *not* verified

- **The colleague's report was not independently re-checked.** Its topology statistics, parsimony
  permutation tests and One Health distance analysis could not be reproduced — the tree and the
  `01_`–`07_` scripts are not on the machine. The assessment of that run is based on `Summaries/`,
  the report text, and the workflow source.
- **Figures were not regenerated** against the corrected countries. The reattribution materially
  changes the Americas lineage (LG3 goes from "USA 31" to "USA 19, Mexico 9, Guatemala 2, Aruba 2…"),
  so any figure or table showing country composition is stale.
- **`clock_signal_check.py` has not been run with `--tree`** on a real tree, only with `--rtt-col`
  against precomputed root-to-tip values. The Newick path is unit-tested but not exercised end to end.
- **The per-cluster variable-site counts come from `cluster_phylogeny_summary.csv`**, not from
  recounting the alignments — those were never copied out.
- **Nothing here has been run on the user's own in-progress results**, which did not exist yet.
