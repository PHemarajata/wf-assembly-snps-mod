# Mapping metadata onto tree tips

How to attach isolate metadata to the tips of a grafted tree so the result can
feed phylogeographic and molecular-clock analyses, and how to repeat it on a new
run. Two scripts do the work:

| script | what it does |
|---|---|
| `bin/build_tip_annotation.py` | joins metadata to tip labels, derives geography/date/source fields, writes TreeTime + LSD2 + mugration inputs |
| `bin/clock_signal_check.py` | root-to-tip regression, whole-tree and per-cluster — run this *before* any dating |

Neither touches the pipeline. They are post-processing on `Summaries/` plus the
grafted tree.

---

## 1. Run it

```bash
python3 bin/build_tip_annotation.py \
  --tree       results/Graft/global_grafted.treefile \
  --metadata   metadata.tsv \
  --clusters   results/Summaries/clusters.tsv \
  --cluster-qc results/Summaries/cluster_phylogeny_summary.csv \
  --outdir     results/Tip_annotation
```

Then:

```bash
python3 bin/clock_signal_check.py --tree results/Graft/global_grafted.treefile --annotation results/Tip_annotation/tip_annotation.tsv --out results/Tip_annotation/clock_signal_by_cluster.tsv
```

If the Newick file is not to hand, `--tips-from-csv FILE:COLUMN` reads the tip
labels from a table instead (that is how the 2026-08-08 run was re-annotated —
its tree was never copied out of the run directory).

Carry forward anything from a previous annotation pass so old and new sit side
by side and can be diffed:

```bash
--carry prev/annotated_tips.csv:Lineage,Lineage_label,Country_Final,Source,root_to_tip
```

Carried columns are prefixed `prior_`, so they can never silently overwrite a
column the script derives.

---

## 2. The join

Tip labels are not `sample_id`. In the 2,773-tip tree they take five shapes:

```
GCF_000756925_1_Australia_Townsville_Queensland   assembly + version + free-text geography
GCA_963561705_1                                   assembly + version
SRR33188703_SPAdes                                run accession + assembler suffix
IP-0009-1-R.fas                                   in-house ID + FASTA extension
IE-0051                                           in-house ID
```

Both sides are reduced to one key:

- `GC[AF]_(\d+)` → `ACC<digits>`. This drops the assembly version and any
  appended text, **and makes `GCF_` and `GCA_` interchangeable** — they denote
  the same assembly and the two files disagree about which prefix to use.
- `(SRR|ERR|DRR)\d+` → the run accession alone. Assembler suffixes are truncated
  differently in different files (`SRR33188703_SPAdes` in the tree,
  `SRR33188703_S` in the metadata).
- Otherwise the label minus any `.fa/.fas/.fasta/.fna` extension.

This reproduces the earlier hand-reconciled join exactly (2,773/2,773 identical
keys) without the four manual fixes it needed, and matches 2,773/2,773 tips to
metadata rows.

**The script aborts if any tip is unmatched.** That default is deliberate: a tip
silently dropped from a phylogeographic model is a country quietly losing a
genome, and nothing downstream will flag it. `--allow-unmatched` overrides, and
writes `unmatched_tips.tsv`.

It also aborts if two metadata `sample_id`s collapse to the same key, since the
join would then be ambiguous.

---

## 3. Geography — which country is the model's state

Two country columns are kept apart, and this distinction is the whole point of
the exercise:

| column | meaning |
|---|---|
| `country_acquired` | where the infection was **acquired** — `--country-col`, default `Country_Final` |
| `country_diagnosis` | where it was **diagnosed** — `--diagnosis-col`, default `Country_Diagnosis` |

A discrete phylogeographic model must use acquisition. Travel-associated
isolates entered under the country of diagnosis create migration events that
never happened: seventeen genomes here were diagnosed in the USA or Switzerland
but acquired in Mexico, Guatemala, Aruba, Costa Rica, Trinidad and Tobago,
Martinique, Panama/Peru, West Africa, Vietnam or Thailand. Attributing them to
the USA manufactures a USA→Latin America link for each one, and inflates the
apparent size of the New World population.

Every disagreement is written to `country_change_log.tsv`. In this dataset the
tip labels corroborate the reattribution independently — they carry the travel
history in the label itself (`GCF_002110945_1_USA_California_ex_Mexico`).

Derived fields:

- **`geo_region`** — ten regions, unchanged from the earlier interpretation so
  figures stay comparable. An unmapped country is a hard error, not a blank;
  add it to `REGION` in the script.
- **`geo_state`** — the trait a DTA/mugration model consumes. Equal to
  `country_acquired`, **except** where the value names more than one country
  ("Panama and Peru") or a continent ("Africa"). Those are not valid states, so
  they are blanked and flagged in `country_ambiguous` rather than guessed. Two
  tips here.
- **`geo_state_pooled`** — countries with fewer than `--min-state-count` tips
  (default 5) replaced by their region. 40 countries → 25 states at the default.
  Worth using: a discrete-trait matrix over 40 states with 15 singletons is
  poorly identified, and singleton states produce transition rates driven
  entirely by the prior.
- **`admin1` / `thai_province`** — ISO-3166-2 name preferred over the free-text
  `Subregion`, because `Subregion` mixes administrative levels (it holds Nakhon
  Phanom *districts* — Na Kae, Tha U-Then, Phon Sawan — where other rows hold
  provinces). ISO-3166-2 correctly aggregates those to the province.
  `admin1_source` records which was used. 2,232 of 2,773 tips resolve.

Country spellings are collapsed through `COUNTRY_ALIASES` (`Viet Nam` and
`Vietnam` are one state, not two).

---

## 4. Dates

`final_collection_dates` is mixed precision. The script parses each into a
precision class plus an interval:

| precision | n | `date_decimal` | `date_lower` / `date_upper` |
|---|---|---|---|
| `day` (YYYY-MM-DD) | 489 | exact | equal to the midpoint |
| `month` (YYYY-MM) | 85 | mid-month | first/last day of the month |
| `year` (YYYY) | 1,842 | YYYY.5 | YYYY.0 / YYYY.99999 |
| `unknown` | 357 | blank | blank |

2,416 of 2,773 tips (87.1%) are dated, spanning 1960–2025.

A year-only record is an **interval, not 1 January**. Collapsing it to the start
of the year biases every root-to-tip regression and every tip-dated clock
towards older tips. Both `tip_dates_treetime.csv` and `tip_dates_lsd2.txt`
therefore carry the bounds, so the tool integrates over the uncertainty:

```
# tip_dates_treetime.csv
name,date,date_lower_bound,date_upper_bound
GCA_015318835_1,2017.5,2017.0,2017.99999

# tip_dates_lsd2.txt
2416
GCA_015318835_1	b(2017.00000,2017.99999)
IE-0051.fa	2025.70959
```

---

## 5. Isolation source

Derived, in order: `EnviSampleID` present → Environmental; `Linked_PatientCase`
present → Clinical; `IE-` prefix → Environmental; `IP-` prefix → Clinical;
otherwise Unspecified. Result: 259 Clinical, 53 Environmental, 2,461
Unspecified.

This differs from the 2026-08-08 annotation, which called 266 Clinical / 319
in-house. The seven extra were `SRR31608433`–`SRR31608440`, CDC public SRA
submissions under `PRJNA908850`, labelled "In-house prospective survey" — and
inconsistently, since `SRR31608436` from the same BioProject and submitter was
labelled public. They are not part of the Nakhon Phanom prospective survey.
Pass `--carry ...:Source` to keep the old call as `prior_source` and compare.

---

## 6. Outputs

| file | use |
|---|---|
| `tip_annotation.tsv` | one row per tip; the master table |
| `tip_dates_treetime.csv` | `treetime --dates` |
| `tip_dates_lsd2.txt` | `lsd2 -d` (also IQ-TREE `--date`) |
| `traits_geography.tsv` | `treetime mugration --states`, or BEAST discrete traits |
| `country_change_log.tsv` | every acquisition≠diagnosis and every ambiguous country |
| `unmatched_tips.tsv` | written only if tips failed to join |
| `tip_annotation_QC.txt` | counts for everything above — read this before using the table |

---

## 7. Before you date anything: run the clock check

`clock_signal_check.py` regresses root-to-tip distance on collection date, whole
tree and per cluster. Read both numbers, because on a grafted tree they mean
different things:

- **Whole tree — not interpretable.** Each cluster's tree was inferred from its
  own alignment, so root-to-tip distance mixes incomparable branch-length
  scales. On the 2026-08-08 tree: r = −0.005, R² = 0.0000 over 2,416 dated
  tips. That is the expected consequence of the grafting, not evidence about
  the clock. **Do not date the grafted tree.**
- **Per cluster — valid.** Within one cluster every branch length comes from one
  alignment, and the shared backbone path from the root cancels out of both the
  correlation and the slope. On the 2026-08-08 tree, of 56 clusters with ≥10
  dated tips and ≥5 distinct dates: 13 have r ≥ 0.3, 24 have 0 < r < 0.3, and
  19 are zero or negative. Median r = +0.107.

**Get the units right before judging a rate.** `IQTREE_ASC` builds the final
per-cluster tree from `GUBBINS_CLUSTER.out.filtered_alignment`, i.e.
`*.filtered_polymorphic_sites.fasta`, under `GTR+ASC` with no genome-wide
`-fconst`. Lewis's ascertainment correction conditions the likelihood on the
sites being variable but cannot know the constant fraction of a 6.8 Mb genome,
so a branch length here is substitutions per **variable** site. Substitutions
per genome per year is `slope × n_variable_sites`, not `slope × genome_length`;
the script does this using `cluster_n_filtered_polymorphic_sites` and reports
both. That conversion is what supplying genome-wide `-fconst` counts would have
achieved at inference time, so **applying it is not a workaround to be replaced
later — the numbers below are already on the corrected scale.**

**Read the rate together with the variable-site count.** On the 2026-08-08 run
the implied rate tracks how much residual diversity the cluster carries:

| cluster | r | variable sites | subs/genome/year |
|---|---|---|---|
| cluster_25 | 0.32 | 2,631 | **1.6** |
| cluster_11 | 0.50 | 2,026 | **4.7** |
| cluster_4 | 0.36 | 20,891 | 43 |
| cluster_44 | 0.33 | 20,467 | 73 |
| cluster_38 | 0.37 | 28,670 | 132 |
| cluster_19 | 0.52 | 29,887 | 341 |
| cluster_2 | 0.41 | 43,682 | 708 |
| cluster_1 | 0.43 | 50,197 | 1,161 |
| cluster_46 | 0.64 | 106,459 | 1,364 |

Roughly 1–10 SNPs/genome/year is the expectation for *B. pseudomallei*.
`cluster_11` and `cluster_25` land in it. Everything above ~20,000 variable
sites does not, by two to three orders of magnitude, and the excess scales with
the site count. **A high r is not sufficient — `cluster_46` has the best
correlation in the set and the least believable rate.**

The reason is that the analysis unit is wrong, not the software. A cluster
holding 106,459 core variable sites across 22 genomes is not one population: it
mixes lineages whose divergence predates the 25-year sampling window by orders
of magnitude. Root-to-tip distance there is dominated by which lineage a genome
belongs to, and correlates with date only because lineages were sampled at
different times. No re-inference fixes that; a different tool will not either.

**Choose the dating unit by depth, not by `cluster_id`:**

1. Take clusters whose filtered-polymorphic-site count is of order 10³ — those
   are single populations on this data. Two of 56 here.
2. For a deep cluster, drop into it and take a subclade with comparable
   diversity rather than dating the whole cluster. Mash clustering at the
   collection-wide threshold (0.005455) is tuned for Gubbins' assumption of
   limited diversity sharing a recent ancestor, which is a much weaker
   requirement than a datable population.
3. Run TreeTime or LSD2 against `tip_dates_*` restricted to that clade. IQ-TREE
   2.2.6 already embeds LSD2 (`--date`), so no new tool is needed.
4. **Sanity-check the rate before believing anything downstream.** If it is not
   within roughly an order of magnitude of 1–10 SNPs/genome/year, the clade is
   still too deep — do not proceed to a dated analysis on it.
5. If you want branch lengths on a genome-wide scale directly out of IQ-TREE
   rather than rescaled after the fact, mask the recombinant positions from the
   `KEEP_INVARIANT_ATCG` full-length alignment with `bin/mask_recombination.py`
   and infer on that. Note that `--iqtree_asc_fallback fconst` does **not** do
   this: `asc_preflight.py` counts only the constant columns remaining inside
   the SNP alignment, not the ~6.8 Mb of invariant genome.

Treat the `usable` rows in `clock_signal_by_cluster.tsv` as candidates to test,
not as results.

---

## 8. Caveats that travel with the table

- **Cluster identity is not portable.** Cluster membership is a property of the
  sample set, not of the population — see `AUDIT_REPORT.md` §G.4, where removing
  0.3% of genomes reorganised most clusters and changed the auto-selected
  threshold. `cluster_id` is only meaningful alongside the threshold and the
  collection it came from. Always report both.
- **Genomes below the cluster-size floor are absent from the tree.** In the
  2026-08-08 run, `clusters.tsv` holds 2,803 genomes in 87 clusters, but only
  the 63 clusters with ≥3 members were built and grafted — 2,773 tips. The 30
  genomes in 24 sub-threshold clusters have no tip and cannot appear in any
  downstream analysis. `metadata rows unused` in the QC file counts them (31
  here: those 30 plus one row with no genome).
- **`prior_root_to_tip` and `prior_terminal_bl` are comparable within a cluster
  only.** Do not threshold on them across the tree.
- **`cluster_confidence_tier` is joined onto every tip.** Tier 1 requires only
  ≥6 isolates and ≥10 filtered polymorphic sites, so it is a floor, not a
  guarantee — `cluster_69` is Tier 1 on 6 isolates and 28 sites. The tier
  columns for filtered sites and recombination blocks are carried so this can be
  judged per tip.
- **`alignment_length` in `cluster_phylogeny_summary.csv` is total sequence
  characters, not columns.** Divide by `seq_count_in_alignment` for the
  per-sequence length (≈6.8 Mb here).
