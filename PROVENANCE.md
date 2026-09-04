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
`5089284`. The tag is placed after all of it so the release carries the checks
that keep the seed from going missing again. Its commit is deliberately not
transcribed here: a tag cannot record the hash of the commit that records it,
and a hand-copied SHA in a table is exactly the kind of number this project has
had to correct too many times. Resolve it from the tag.

`v1.1.0-mod` was moved forward more than once on the day it was created, first
from `e9e7753`. That placement was on a red build: the workflow there installed
`latest-stable`, which is 26.x, and 26.x could not then parse this config.
It can now, as of 2026-09-04; the CI matrix pins versions rather than
tracking `latest-stable`, so a future strict-parser change fails a named
job instead of every job at once. Nothing
referenced the tag while it stood, and **no pipeline code differs across any of
those commits**; only CI and this file changed. Recorded because a moved tag
nobody writes down is indistinguishable from a tag someone misread. Resolve the
final placement from the tag, not from this paragraph.

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

## Nextflow 26.x compatibility, and the six things that blocked it

**Resolved 2026-09-04.** `nextflow config .` now succeeds on **24.10.5, 25.04.6,
25.10.0 and 26.04.6**, across all eleven profiles: 44 combinations, all checked.
CI tests all four versions.

26.x ships a strict config parser. It rejects function definitions, variable
declarations, control flow, and bare environment variables in a config file.
This repository had all four, and each was hidden behind the one before it, so
the count went up as the work went on rather than down. In order:

| # | rejected | where | replaced with |
|---|---|---|---|
| 1 | `def check_max(obj, type)` | `nextflow.config` | the `resourceLimits` process directive, in `conf/base.config` and all six profiles |
| 2 | nested `if`/`each` building the include list | `nextflow.config` | the five `includeConfig` lines it always resolved to |
| 3 | `def trace_timestamp = ...` | `nextflow.config` | `params.trace_timestamp`, in `conf/params.config` |
| 4 | `"${LAB_HOME}/..."` | `nextflow.config`, `conf/profiles.config` | `System.getenv('LAB_HOME')`, which behaves identically on all four versions |
| 5 | `'parsnp' { description: ... }` | `conf/workflows.config` | unquoted scope names; `description:` was a Groovy label, never a key |
| 6 | `sge_options = { ...task... }` | `conf/profiles/rosalind_hpc.config` | the same closure on `process.clusterOptions`, where `task` is defined |

None of it changes what any profile allocates. That was measured, not assumed.

### The one way this could have changed the science, and how it was caught

`resourceLimits` written as a plain list is **bound once, when the file is
parsed**. `check_max()` read `params.max_*` per task. So the obvious
translation silently stops honouring a ceiling set in a later `-c` overlay:

| ceiling set in a later config | old `check_max` | `resourceLimits = [ ... ]` | `resourceLimits = { [ ... ] }` |
|---|---|---|---|
| `max_cpus = 3`, `max_memory = '5.GB'` | 3 cpus, 5 GB | **12 cpus, 16 GB — ignored** | 3 cpus, 5 GB |

The list form was written first and passed every check that had been run
against it, because a command-line `--max_cpus` *is* honoured either way and
that is what was tested. It was caught by a negative control: deliberately
lowering the ceiling and confirming the measurement moved. It did not.

**All seven declarations are closures, and CI asserts it.** This is the fourth
member of the family in `HANDOFF_2026-09-04.md`: a check that passes because it
was never able to fail.

### How the equivalence was verified

Not by reading the two implementations and finding them alike.

- **Config output.** `nextflow config` on 25.04.6, ten profiles, against a clean
  clone of `v1.1.0-mod`. Identical except the closure bodies, the added
  `trace_timestamp`, and the clone's own path.
- **Resolved resources.** A harness of seven processes, one per label, failing
  twice so each runs at `task.attempt` 1, 2 and 3, printing the cpus, memory and
  time it actually received. Five profiles, 21 measurements each, old against
  new: **identical in all 105**. The ceiling demonstrably binds — under
  `local_workstation`, `process_high` cpus go 6, 12, 12 against `max_cpus = 12`,
  and `process_high_memory` goes 32, 64, 64 GB against `max_memory = '64.GB'`.
- **Sensitivity.** The same harness with the ceiling lowered moves the numbers,
  so it is not measuring a constant.

Item 2 removed the only reader of `params.workflows`. The map is kept as the
written-down module list, and the static includes are checked against it.

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
