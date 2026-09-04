# Provenance: which code produced which result

This mapping previously existed only in the analysis repository, which is the
wrong home for it: someone handed this pipeline alone could not establish what
produced the published numbers. It lives here now.

---

## The reported analysis

| | |
|---|---|
| **Release** | **`v1.0.5-mod`**, annotated tag at `79ab645` |
| **Commit** | `79ab645`, 2026-08-16 |
| **Nextflow** | 25.04.6, build 5954 |
| **Nextflow `scriptId`** | `e09a5c4eadba2c5984f6790095423ee4` |
| **Run name / session** | `agitated_coulomb` / `c90e1105-5b12-455e-9b31-4ecde888d559` |
| **Basis produced** | 85 units, 2,340 genomes; median in-window r/m 7.70 |
| **Cross-hardware control** | `insane_jennings`, NVIDIA DGX A100, 88 units, 2,342 genomes |

### How the commit was established, and what it does not prove

`nextflow run .` executes a **local directory**, so Nextflow records a `scriptId`,
a content hash of `main.nf`, rather than a git commit. **`e09a5c4ead…` is not a
git SHA**; `git cat-file -t e09a5c4ead` fails, and anyone treating it as one will
lose an afternoon.

The commit was established by bracketing: `79ab645` is the last commit before the
run, the next two (`f1167f4`, `19d764c`, both 2026-08-19) postdate it, and
**neither touches `main.nf`**. That is consistent with the recorded `scriptId` but
does not by itself discriminate between them, which is why the tag exists.

**Both the reported run and the cross-hardware control record the same
`scriptId`.** They therefore executed byte-identical pipeline code, which is a
stronger claim than a shared commit, since it hashes the workflow definition
itself. Containerised tool versions match across both; the runs differ in exactly
two respects, Nextflow version and resource profile.

### ⚠ The reported run is not seed-reproducible

`v1.0.5-mod` predates both `gubbins_seed` and `gubbins_deterministic`. Re-running
it under either produces a **different** run rather than validating the pinned
one. This is not recoverable and should be stated once, plainly, in any methods
section citing this release.

The tag's own message additionally asserts that the commit "predates the
`gubbins_seed` fix". **That wording is wrong**: on 2026-08-16 no such fix existed
anywhere, and none was added until 2026-09-04. The tag message cannot be corrected
without retagging a published release, so this file is the correction of record.

---

## Releases

| tag | commit | what changed | reproducible? |
|---|---|---|---|
| `v1.0.3-mod` | `2b9557b` | first release | |
| `v1.0.4-mod` | `0669624` | `graft_trees.py` | |
| **`v1.0.5-mod`** | **`79ab645`** | **the reported analysis** | no — see above |
| `v1.1.0-mod` | `git rev-list -n1 v1.1.0-mod` | determinism: `gubbins_seed`, `gubbins_deterministic`; first CI | yes, with `gubbins_deterministic = true` |

The determinism work itself is `0543892` (`--seed`), `4fd7b22`
(`gubbins_deterministic`) and `a49eac5` (the classic path); CI and this file are
`5089284`. The tag is placed after all four so the release carries the checks
that keep the seed from going missing again. Its commit is deliberately not
transcribed here: a tag cannot record the hash of the commit that records it,
and a hand-copied SHA in a table is exactly the kind of number this project has
had to correct too many times. Resolve it from the tag.

**The manifest at `79ab645` self-reports `version = '1.0.3-mod'`**, so run logs
from the reported analysis print that string while the `v1.0.3-mod` tag points at
a much older 2025 commit. Left uncorrected deliberately: bumping the manifest
would have changed the SHA the manuscript pins.

---

## Determinism, and what each parameter is for

The two parameters do different jobs and neither substitutes for the other.

| parameter | default | what it does |
|---|---|---|
| `gubbins_seed` | `20260904` | Removes a **silent unit loss**. Without it Gubbins draws RAxML's parsimony seed from an unseeded `randint(0, 10000)`; that is `0` about 1 time in 10,001, RAxML rejects it, and Gubbins reports only "Unable to fit model to data". With `errorStrategy 'ignore'` the unit is dropped and the run still exits 0. ~16% chance per panel. **Required always.** |
| `gubbins_deterministic` | `false` | Forces `--threads 1`, which is what actually gives byte-reproducibility. **Costs ~2x.** |

### Measured, 2026-09-04

Ten units, two runs each on the same alignment and the same code, comparing
`per_branch_statistics.csv`, `recombination_predictions.gff`,
`node_labelled.final_tree.tre` and `filtered_polymorphic_sites.fasta` byte for
byte:

| configuration | identical pairs |
|---|---|
| no seed, `--threads 4` | 4 / 10 |
| `--seed`, `--threads 4` | 5 / 10 |
| **`--seed`, `--threads 1`** | **10 / 10** |

**The seed alone is worth one unit in ten**, within noise of nothing. The five
that agreed under a seed are the same five that agreed without one; they are
stable regardless. Thread count, not the seed, is the dominant source of
run-to-run variation.

Cost: 1.28x at 8 taxa, 1.98x at 37. Units reach 159, so expect more.

This matches the result measured independently for IQ-TREE in the analysis
repository: `-seed` alone gave a different tree every run, `-seed` with `-T 1`
gave one identical tree across three. **Assume any threaded tree builder in this
stack is non-deterministic regardless of seed unless measured otherwise.**

Full write-up: `DETERMINISM_DEMONSTRATION_2026-09-04.md` in the analysis
repository.

### Both Gubbins paths are wired, and one was not used

`GUBBINS_CLUSTER` (clustered/scalable workflow) produced the reported analysis.
`RECOMBINATION_GUBBINS` (classic workflow, `subworkflows/local/recombination.nf`)
**produced no reported result**, and was wired for determinism anyway so that one
parameter governs both paths rather than two behaving differently for reasons
nobody documented.

Doing so surfaced a defect worth recording: that module passed **no `--threads` at
all**, so Gubbins used its own default of 1. The classic path had been
single-threaded, and after `gubbins_seed` accidentally deterministic, while
carrying a `process_medium` label and CPUs it never used. Implicit behaviour
arrived at by omission is the same class of defect as an argument default pointing
at a specific run. It is now explicit: allocated CPUs by default, 1 under
`gubbins_deterministic`.

---

## ⚠ Known incompatibility: Nextflow 26.x cannot parse this config

`nextflow config .` fails outright on Nextflow **26.04.6**:

```
Error nextflow.config:363:14: Unexpected input: '('
│ 363 | def check_max(obj, type) {
```

26.x ships a strict config parser that no longer allows function definitions in
`nextflow.config`. `check_max` is the nf-core resource-capping helper and is
used throughout the profile blocks, so this is not a one-line deletion; removing
it changes how every profile's resource ceilings are computed, and those
ceilings are already known to be sized for small units.

**Nothing in this repository has been run on 26.x, and nothing should be until
that is fixed.** Verified working: 24.10.5, 25.04.6 and 25.10.0, in both the
working tree and a fresh clone. CI tests 25.04.6 and 25.10.0 specifically,
because those are the two versions the reported results were produced on.

This was found by CI on its first run, which is the argument for having it.

---

## Pinned tool versions

Pins that alter scientific output, and must not be bumped casually:

| tool | pin | why it is pinned |
|---|---|---|
| **Gubbins** | `3.4.3--py310h5140242_0` | 3.4.2 made `--invariant-site-correction` **optional and defaulted it off**. Bumping from 3.3.5 without the flag silently drops a correction 3.3.5 always applied — a results change with no error. The flag and the version move together. |
| **RAxML** | via Gubbins | `gubbins_tree_builder` and `gubbins_first_tree_builder` are both pinned to `raxml`. Measured on 6 units x 2 replicons: rapidnj **systematically underestimates** r/m (median ratio 0.922, 11/12 low, sign test p = 0.0063, worst case 45.5% low). Do not reintroduce a distance-based builder. |
| IQ-TREE | `2.2.6--h21ec9f0_0` | |
| PopPUNK | `2.7.6--py310h4d0eb5b_0` | |
| SKA2 | `0.3.7--h4349ce8_2` | |
| parsnp | `1.7.4--hdcf5f25_2` | |
| snippy | `staphb/snippy:4.6.0` | |
| ClonalFrameML | `snads/clonalframeml@sha256:bc00db…` | |

`gubbins_filter_percentage` is also explicit at 25 rather than left to Gubbins'
invisible default: measured to silently exclude 7 of 30 genomes on a real
SKA-mapped set, and default versus 100 produced **different final trees**.

---

## Reproducing a result from this repository

```bash
git checkout v1.1.0-mod          # or v1.0.5-mod for the reported basis
nextflow run . -profile <profile> --gubbins_deterministic true \
        --input <samplesheet> --outdir <out>
```

Two things to check afterwards, because neither is implied by a zero exit code:

1. **Verify per unit, not by exit code.** `errorStrategy 'ignore'` means a dropped
   unit still returns 0. Read `gubbins_exit_code`, `iqtree_status` and
   `confidence_tier` from `Summaries/cluster_phylogeny_summary.csv` and confirm
   the replicon-unit count matches what was requested.
2. **Give every concurrent run its own working directory.** Gubbins writes scratch
   to the working directory regardless of `--prefix`, and concurrent runs sharing
   one collide. That failure is invisible in single-run testing, appears only under
   concurrency, and reports as a problem with the input rather than the
   invocation. It cost this project three wrong conclusions.
