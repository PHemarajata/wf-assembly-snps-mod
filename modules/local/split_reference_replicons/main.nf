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
    def max_replicons = (params.max_replicons ?: 4) as int
    """
    set -euo pipefail
    mkdir -p replicons

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
