process SPLIT_REFERENCE_REPLICONS {

    tag { "${cluster_id}" }
    label "process_low"
    container "ubuntu:jammy"

    input:
    tuple val(cluster_id), path(reference)

    output:
    tuple val(cluster_id), path("replicons/*.fa"), emit: replicons
    path("versions.yml")                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Split a COMPLETE reference into one FASTA per replicon (contig) so Gubbins
    // runs per replicon. Gubbins' 0.1-10 kb sliding window would otherwise scan
    // across contig junctions and call spurious recombination there, and
    // snp-sites hardcodes CHROM to "1" -- both are wrong on a multi-contig input.
    // This is the workflow equivalent of the manual pipeline's replicon split.
    //
    // GUARD: a reference with more than max_replicons contigs is a DRAFT, not a
    // complete multi-replicon genome. Splitting it would run Gubbins per contig
    // (meaningless) and is almost always a sign the wrong reference was supplied
    // (e.g. a 115-contig medoid draft). Fail loudly rather than produce plausible
    // nonsense. B. pseudomallei has 2 chromosomes; the default 4 leaves room for
    // plasmids.
    // SIZE GATE, separate from the count guard above. max_replicons is a COUNT
    // check and cannot catch a size problem. Measured: strain_12's reference
    // GCF_027856855.2 has 4 contigs -- the two real chromosomes at 3.95 Mb and
    // 3.09 Mb, plus two of 2,595 bp and 2,533 bp. Being 4 contigs it passed
    // max_replicons cleanly, then each 2.5 kb "replicon" spawned its own 37
    // snippy jobs (148 wasted tasks) and its own analysis unit. A 2.5 kb
    // alignment cannot yield a meaningful r/m, and units that cannot be
    // partitioned meaningfully before Gubbins must be dropped, not reported.
    // In that run both died at SNIPPY_CORE_GATHER anyway (the contig is present
    // in <2 of 37 genomes), so the only products were wasted compute and two
    // spurious failure rows in the summary.
    def max_replicons = (params.max_replicons ?: 4) as int
    def min_replicon_length = (params.min_replicon_length == null ? 100000 : params.min_replicon_length) as long
    """
    set -euo pipefail
    mkdir -p replicons

    # What Gubbins appends to the unit name before handing it to RAxML as -n.
    # Measured against a real run: a 99-character unit produced a 136-character
    # run id, which is exactly 99 + 37.
    RECON_SUFFIX=".core.full.iteration_1_reconstruction"

    n_contigs=\$(grep -c '^>' "${reference}")
    echo "Reference ${reference} has \${n_contigs} contig(s); max_replicons=${max_replicons}" >&2
    if [ "\${n_contigs}" -gt "${max_replicons}" ]; then
        echo "ERROR: reference for cluster ${cluster_id} has \${n_contigs} contigs (> ${max_replicons})." >&2
        echo "       This is a draft assembly, not a complete multi-replicon genome." >&2
        echo "       Replicon splitting requires a complete reference (supply one via" >&2
        echo "       --cluster_references, or raise --max_replicons only if you are certain" >&2
        echo "       every contig is a real replicon)." >&2
        exit 1
    fi

    # One file per contig, named for a sanitized version of the first header token
    # so the replicon id is filesystem-safe and stable (e.g. 'NC_006350.1' ->
    # 'NC_006350_1'). awk writes each record to its own file.
    awk -v OFS='' '
        /^>/ {
            id=substr(\$1,2)
            gsub(/[^A-Za-z0-9]/,"_",id)
            out="replicons/" id ".fa"
            print ">" id > out
            next
        }
        { print \$0 > out }
    ' "${reference}"

    # Drop anything below min_replicon_length. Report every drop with its length:
    # a silently discarded contig is indistinguishable from one that was never
    # there, and the whole point is that the operator can see what was excluded.
    dropped=0
    for f in replicons/*.fa; do
        len=\$(grep -v '^>' "\$f" | tr -d '\\n' | wc -c)
        if [ "\$len" -lt "${min_replicon_length}" ]; then
            echo "DROPPED replicon \$(basename "\$f" .fa): \${len} bp < min_replicon_length=${min_replicon_length}" >&2
            rm -f "\$f"
            dropped=\$(( dropped + 1 ))
        fi
    done
    if [ "\$dropped" -gt 0 ]; then
        echo "Dropped \${dropped} sub-threshold replicon(s) for cluster ${cluster_id}" >&2
    fi

    # LENGTH GATE. The workflow keys each analysis unit as
    # <cluster_id>__<replicon id>, and Gubbins hands RAxML
    # "<unit>.core.full.iteration_N_reconstruction" as its -n run id.
    #
    # raxmlHPC v8 SEGFAULTS (exit 139) on a -n run id of 128 characters or more.
    # Measured directly: identical inputs, only -n length varied -- 127 exits 0,
    # 128 exits 139 -- and the -w path length is irrelevant. RAxML contains the
    # string 'Error: run id after "-n" is too long' but crashes before printing
    # it, and Gubbins wraps the call in a bare `except` that reports only
    # "Unable to fit model to data". So the whole failure is silent and looks
    # like a model-fitting problem. It is not: it is the length of the
    # reference's FASTA defline, which is why no genome quality metric --
    # fastANI, contiguity, N50, ambiguous bases, GC, duplication ratio,
    # misassemblies -- separates the references that fail from those that work.
    #
    # Measured on a real 82-unit partition before the deflines were normalized:
    # 40 of 164 replicon-units (24%) were over the limit. Fail here, loudly,
    # naming the offender -- do not silently truncate, because that would change
    # unit identities the operator did not ask to change.
    max_unit=\$(( 127 - \${#RECON_SUFFIX} ))
    for f in replicons/*.fa; do
        rid=\$(basename "\$f" .fa)
        unit="${cluster_id}__\${rid}"
        if [ "\${#unit}" -gt "\$max_unit" ]; then
            echo "ERROR: analysis unit id is too long for RAxML." >&2
            echo "       unit    : \${unit}" >&2
            echo "       length  : \${#unit} characters; the maximum is \${max_unit}" >&2
            echo "       because : Gubbins appends '\${RECON_SUFFIX}' and raxmlHPC" >&2
            echo "                 segfaults at a -n run id of 128 characters or more." >&2
            echo "       fix     : shorten the reference's FASTA deflines. The replicon" >&2
            echo "                 id is the first header token, so a defline like" >&2
            echo "                 '>GCF_000954175_1_1' gives a short, stable unit id." >&2
            exit 1
        fi
    done

    if [ -z "\$(ls -A replicons/ 2>/dev/null)" ]; then
        echo "ERROR: every replicon of cluster ${cluster_id} fell below" >&2
        echo "       --min_replicon_length ${min_replicon_length}. Either the reference is" >&2
        echo "       not a genome, or the threshold is set too high." >&2
        exit 1
    fi

    echo "Wrote \$(ls replicons/ | wc -l) replicon file(s):" >&2
    ls -1 replicons/ >&2

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    awk: \$(awk --version 2>&1 | head -1 | sed 's/^/    /' || echo "    gawk")
END_VERSIONS
    """

    stub:
    // Run the REAL split on the staged reference: it is trivial awk, and a stub
    // that emitted fixed toy replicons made downstream snippy-core choke on 4 bp
    // input. Splitting the real reference keeps the stub chain faithful (correct
    // replicon count, real sequence) at negligible cost.
    """
    set -euo pipefail
    mkdir -p replicons
    awk -v OFS='' '
        /^>/ {
            id=substr(\$1,2)
            gsub(/[^A-Za-z0-9]/,"_",id)
            out="replicons/" id ".fa"
            print ">" id > out
            next
        }
        { print \$0 > out }
    ' "${reference}"
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    awk: stub
END_VERSIONS
    """
}
