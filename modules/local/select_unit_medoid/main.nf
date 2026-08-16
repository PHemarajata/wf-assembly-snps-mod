#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * SELECT_UNIT_MEDOID
 *
 * Pick one representative genome per analysis unit: the MEDOID, the member with
 * the smallest mean SNP distance to the rest of its own unit.
 *
 * WHY THIS EXISTS. STEP 5/6 build a whole-collection backbone from per-cluster
 * medoids, and that path was skipped entirely in curated mode on the grounds
 * that "curated mode has no medoids -- its reference is external and not a
 * cluster member". The first half is wrong. The reference indeed cannot be a
 * backbone tip, but every unit still HAS a most-typical member; nothing about
 * curated mode removes it. What curated mode lacks is the MASH MATRIX the
 * default path used to find it.
 *
 * It does not need one. By this point each unit already has a filtered
 * polymorphic-sites alignment from Gubbins -- recombination removed, a few
 * thousand columns -- and pairwise distance on that is both cheaper than Mash
 * and a better measure of centrality, because it is computed on the clonal
 * frame rather than on whole-genome k-mer content.
 *
 * THE REFERENCE TAXON IS EXCLUDED. It sits outside the population by
 * construction, so it would win or lose the centrality contest for reasons that
 * have nothing to do with the unit. It is also not a member and must never
 * become a backbone tip.
 */

process SELECT_UNIT_MEDOID {

    tag "${cluster_id}"
    label 'process_low'
    container "quay.io/biocontainers/python:3.10"

    input:
    tuple val(cluster_id), path(filtered_aln), val(sample_ids), path(assemblies)

    output:
    tuple val(cluster_id), path("medoid/*.fa"), emit: medoid
    tuple val(cluster_id), path("${cluster_id}.medoid.txt"), emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail
    mkdir -p medoid

    python3 <<'PY'
import os, glob, sys, collections

aln = "${filtered_aln}"
cluster = "${cluster_id}"

seqs, name = {}, None
buf = []
with open(aln) as fh:
    for line in fh:
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = line[1:].strip().split()[0]
            buf = []
        else:
            buf.append(line.strip())
if name:
    seqs[name] = "".join(buf)

# The external reference is not a member of the population and must not be a
# backbone tip; drop it before measuring centrality.
seqs.pop("Reference", None)
if not seqs:
    sys.exit("ERROR: no sequences in %s after removing the Reference taxon" % aln)

names = sorted(seqs)
if len(names) == 1:
    medoid = names[0]
    mean_d = 0.0
else:
    # Pairwise mismatch counts over the filtered polymorphic sites. Ambiguous
    # positions in either sequence are skipped rather than counted as
    # differences, so a taxon with more missing data is not pushed to the edge.
    ACGT = set("ACGTacgt")
    cols = list(zip(*[seqs[n] for n in names]))
    tot = collections.Counter()
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            d = 0
            for c in cols:
                a, b = c[i], c[j]
                if a in ACGT and b in ACGT and a.upper() != b.upper():
                    d += 1
            tot[names[i]] += d
            tot[names[j]] += d
    denom = len(names) - 1
    medoid = min(names, key=lambda n: (tot[n] / denom, n))
    mean_d = tot[medoid] / denom

# Resolve the medoid back to its assembly. Names are matched on the stem so a
# staged '<id>.fasta' or '<id>.fa' both resolve.
staged = {}
for f in glob.glob("*"):
    if f.endswith((".fasta", ".fa", ".fna")) and not f.startswith("medoid"):
        staged.setdefault(os.path.basename(f).rsplit(".", 1)[0], f)

src = staged.get(medoid)
if src is None:
    for stem, f in staged.items():
        if stem.startswith(medoid) or medoid.startswith(stem):
            src = f
            break
if src is None:
    sys.exit("ERROR: medoid %s has no staged assembly (have %d: %s)"
             % (medoid, len(staged), sorted(staged)[:5]))

# The backbone tip is named for the BASE UNIT -- not the genome, and not the
# replicon-compound id. With --split_replicons the cluster_id arrives as
# "<unit>__<replicon>", and both replicons of a unit select the SAME genome, so
# naming the file after the compound id would put the replicon suffix (and a
# reference accession) into every tip label of a tree whose tips ARE units.
# Splitting on the first "__" is a no-op when replicons are not split.
tip = cluster.split("__")[0]
with open(os.path.join("medoid", tip + ".fa"), "w") as out:
    with open(src) as fh:
        for line in fh:
            out.write(line)

with open(cluster + ".medoid.txt", "w") as fh:
    fh.write("cluster_id\\tmedoid\\tn_members\\tmean_snp_distance\\n")
    fh.write("%s\\t%s\\t%d\\t%.3f\\n" % (cluster, medoid, len(names), mean_d))
print("%s: medoid %s of %d members, mean SNP distance %.3f"
      % (cluster, medoid, len(names), mean_d))
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
