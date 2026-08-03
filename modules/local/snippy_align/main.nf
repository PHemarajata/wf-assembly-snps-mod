#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Snippy per-cluster alignment, SCATTER + GATHER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    WHY THIS WAS SPLIT
    The previous SNIPPY_ALIGN was one Nextflow task per CLUSTER that ran snippy
    SERIALLY in a bash for-loop over every sample in that cluster. Cluster-level
    concurrency was then capped again by `maxForks = 2`, so with 2000 genomes in
    clusters of <=50 the wall clock was (samples per cluster) x (snippy runtime)
    x (clusters / 2) -- Nextflow could not overlap samples at all.

    SNIPPY_SCATTER emits one task per (cluster, sample) so the scheduler owns all
    the parallelism, and SNIPPY_CORE_GATHER does the snippy-core join per cluster.
    MODEL (not a measurement -- see below): ~3.5x at every scale, i.e. 2000 genomes
    ~49 h -> ~14 h at 3 min/sample. This has NOT been measured because snippy is
    only distributed as a container and Docker was unavailable in the development
    sandbox; the channel restructuring itself was executed end-to-end against a
    mock `snippy`/`snippy-core` that produces realistically shaped outputs.

    CORRECTNESS FIXES
    1. The old code located each sample's assembly with a substring test:
           for file in *.fa *.fasta *.fna; do if [[ "$file" == *"$sample"* ]] ...
       That mis-binds any sample id that is a prefix/substring of another (e.g.
       BP_12 matching BP_123.fasta), silently aligning the wrong genome under the
       right label. Sample id and assembly path now travel together as an explicit
       channel tuple, so no filename guessing happens at all.
    2. `|| { echo "Warning: Snippy failed"; continue; }` and the
       ">rep / N" one-base alignment fallbacks are gone. A failed sample now fails
       its task; a cluster that loses samples fails the gather rather than quietly
       producing a smaller (or 1 bp) alignment that flows into Gubbins.
*/

process SNIPPY_SCATTER {
    tag "cluster_${cluster_id}:${sample_id}"
    label 'process_medium'
    conda "bioconda::snippy=4.6.0 bioconda::samtools=1.9"
    container "staphb/snippy:4.6.0"

    input:
    tuple val(cluster_id), val(sample_id), path(assembly), val(representative_id), path(reference)

    output:
    tuple val(cluster_id), val(sample_id), path("snippy_out/${sample_id}"), emit: sample_dir
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    set -euo pipefail

    echo "Snippy: cluster=${cluster_id} sample=${sample_id} reference=${representative_id}"

    if [ ! -s "${assembly}" ]; then
        echo "ERROR: assembly for sample ${sample_id} is missing or empty: ${assembly}" >&2
        exit 1
    fi
    if [ ! -s "${reference}" ]; then
        echo "ERROR: reference for cluster ${cluster_id} is missing or empty: ${reference}" >&2
        exit 1
    fi

    mkdir -p snippy_out

    # No fallback: a snippy failure must fail the task. Silently continuing past a
    # failed sample is how a cluster ends up with an alignment of unknown membership.
    snippy \\
        --outdir "snippy_out/${sample_id}" \\
        --ref "${reference}" \\
        --ctgs "${assembly}" \\
        --cpus ${task.cpus} \\
        --force \\
        ${args}

    if [ ! -s "snippy_out/${sample_id}/snps.vcf" ]; then
        echo "ERROR: snippy produced no snps.vcf for ${sample_id}" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snippy: \$(snippy --version 2>&1 | head -n1 | sed 's/^.*snippy //')
    END_VERSIONS
    """
}

process SNIPPY_CORE_GATHER {
    tag "cluster_${cluster_id}"
    label 'process_medium'
    conda "bioconda::snippy=4.6.0 bioconda::samtools=1.9"
    container "staphb/snippy:4.6.0"

    publishDir "${params.outdir}/Clusters/cluster_${cluster_id}",
               mode: params.publish_dir_mode,
               pattern: "*.core.{full.aln,tab}"

    input:
    tuple val(cluster_id), val(sample_ids), path(sample_dirs, stageAs: 'snippy_in/*'), val(representative_id), path(reference)

    output:
    tuple val(cluster_id), path("${cluster_id}.core.full.aln"), emit: core_alignment
    tuple val(cluster_id), path("${cluster_id}.core.tab"),      emit: core_tab
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // See conf/params.config: with the default medoid reference, snippy-core's
    // "Reference" taxon duplicates a sample already in the cluster. Unset means
    // "drop it only when it is a duplicate"; an explicit value forces the choice,
    // and false restores the legacy alignment.
    def using_global_ref = (params.use_global_reference && params.ref) ? true : false
    def drop_ref = (params.drop_reference_taxon == null)
                     ? !using_global_ref
                     : (params.drop_reference_taxon.toString().toLowerCase() in ['true','1','yes'])
    """
    set -euo pipefail

    echo "snippy-core: cluster=${cluster_id} reference=${representative_id}"
    echo "Samples staged: ${sample_ids}"

    # -L: Nextflow stages the per-sample directories as SYMLINKS, so a bare
    # `find -type d` counts zero.
    n_dirs=\$(find -L snippy_in -mindepth 1 -maxdepth 1 -type d | wc -l)
    echo "Staged snippy sample directories: \$n_dirs"
    if [ "\$n_dirs" -lt 2 ]; then
        echo "ERROR: cluster ${cluster_id} has only \$n_dirs snippy sample directories; \\
snippy-core needs >=2 to produce a usable alignment" >&2
        exit 1
    fi

    snippy-core \\
        --ref "${reference}" \\
        --prefix "${cluster_id}" \\
        ${args} \\
        snippy_in/*/

    # snippy-core writes <prefix>.full.aln (whole-genome) and <prefix>.tab.
    # The whole-genome alignment is the one Gubbins needs -- never substitute
    # <prefix>.aln, which is SNP-only.
    if [ -s "${cluster_id}.full.aln" ]; then
        if ${drop_ref}; then
            # Drop the "Reference" record: with the medoid reference it is the same
            # genome as one of the samples below it in this very alignment.
            n_before=\$(grep -c "^>" "${cluster_id}.full.aln")
            awk '/^>/ { keep = (\$0 != ">Reference") } keep' \\
                "${cluster_id}.full.aln" > "${cluster_id}.core.full.aln"
            n_after=\$(grep -c "^>" "${cluster_id}.core.full.aln")
            echo "Dropped duplicate Reference taxon: \$n_before -> \$n_after sequences"
            if [ "\$n_after" -ne \$(( n_before - 1 )) ]; then
                echo "ERROR: expected to remove exactly one 'Reference' record, went \\
\$n_before -> \$n_after" >&2
                exit 1
            fi
            # CLUSTER_GENOMES only emits clusters of >=3 samples, so >=3 must survive.
            if [ "\$n_after" -lt 3 ]; then
                echo "ERROR: only \$n_after taxa remain in cluster ${cluster_id} after \\
dropping Reference; a tree needs >=3. Re-run with --drop_reference_taxon false to \\
restore the legacy alignment." >&2
                exit 1
            fi
        else
            mv "${cluster_id}.full.aln" "${cluster_id}.core.full.aln"
        fi
    else
        echo "ERROR: snippy-core did not produce ${cluster_id}.full.aln (the whole-genome \\
alignment). Gubbins cannot use a SNP-only alignment, so no substitution is made." >&2
        exit 1
    fi

    if [ -s "${cluster_id}.tab" ]; then
        mv "${cluster_id}.tab" "${cluster_id}.core.tab"
    else
        echo "ERROR: snippy-core did not produce ${cluster_id}.tab" >&2
        exit 1
    fi

    n_seqs=\$(grep -c '^>' "${cluster_id}.core.full.aln")
    echo "Core alignment sequences: \$n_seqs"
    if [ "\$n_seqs" -lt 3 ]; then
        echo "ERROR: cluster ${cluster_id} core alignment has \$n_seqs sequences (<3)" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snippy: \$(snippy --version 2>&1 | head -n1 | sed 's/^.*snippy //')
    END_VERSIONS
    """
}

// Retained only so that `include { SNIPPY_ALIGN }` in any other workflow keeps
// resolving. The recombination-aware workflow now uses SNIPPY_SCATTER +
// SNIPPY_CORE_GATHER; this monolithic version is deprecated and unparallelised.
process SNIPPY_ALIGN {
    tag "cluster_${cluster_id}"
    label 'process_high'
    conda "bioconda::snippy=4.6.0 bioconda::samtools=1.9"
    container "staphb/snippy:4.6.0"

    input:
    tuple val(cluster_id), val(sample_ids), path(assemblies), val(representative_id), path(reference)

    output:
    tuple val(cluster_id), path("${cluster_id}.core.full.aln"), emit: core_alignment
    tuple val(cluster_id), path("${cluster_id}.core.tab"), emit: core_tab
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bash <<'EOF'
    echo "Running Snippy alignment for cluster ${cluster_id}"
    echo "Reference: ${representative_id}"
    echo "Samples: ${sample_ids}"

    # Create working directory
    mkdir -p snippy_work

    for sample in \$(echo ${sample_ids} | sed 's/[][]//g; s/,/ /g'); do
        echo "Processing sample: \$sample"

        # Find the assembly file for this sample
        assembly_file=""
        for file in *.fa *.fasta *.fna; do
            if [[ "\$file" == *"\$sample"* ]] || [[ "\$(basename "\$file" .fa)" == "\$sample" ]] || [[ "\$(basename "\$file" .fasta)" == "\$sample" ]] || [[ "\$(basename "\$file" .fna)" == "\$sample" ]]; then
                assembly_file="\$file"
                break
            fi
        done

        if [[ -z "\$assembly_file" ]]; then
            echo "Warning: Could not find assembly file for sample \$sample"
            continue
        fi

        echo "Found assembly: \$assembly_file"

        # Skip if this is the reference sample
        if [[ "\$sample" == "${representative_id}" ]]; then
            echo "Skipping reference sample \$sample"
            continue
        fi

        # Run snippy
        snippy \\
            --outdir snippy_work/\$sample \\
            --ref ${reference} \\
            --ctgs \$assembly_file \\
            --cpus ${task.cpus} \\
            --force \\
            ${args} || {
            echo "Warning: Snippy failed for sample \$sample"
            continue
        }
    done

    # Collect all snippy output directories
    snippy_dirs=()
    for dir in snippy_work/*/; do
        if [[ -d "\$dir" ]] && [[ -f "\$dir/snps.vcf" ]]; then
            snippy_dirs+=("\$dir")
        fi
    done

    if [[ \${#snippy_dirs[@]} -eq 0 ]]; then
        echo "Warning: No successful Snippy runs found"
        # Create empty alignment
        echo ">${representative_id}" > ${cluster_id}.core.full.aln
        echo "N" >> ${cluster_id}.core.full.aln
        touch ${cluster_id}.core.tab
    else
        echo "Found \${#snippy_dirs[@]} successful Snippy runs"

        # Run snippy-core to generate core alignment
        snippy-core \\
            --ref ${reference} \\
            --prefix ${cluster_id} \\
            \${snippy_dirs[@]} || {
            echo "Warning: snippy-core failed, creating minimal alignment"
            echo ">${representative_id}" > ${cluster_id}.core.full.aln
            echo "N" >> ${cluster_id}.core.full.aln
            touch ${cluster_id}.core.tab
        }

        # Rename snippy-core output files to expected names
        if [[ -f "${cluster_id}.full.aln" ]]; then
            mv "${cluster_id}.full.aln" "${cluster_id}.core.full.aln"
        elif [[ -f "${cluster_id}.aln" ]]; then
            mv "${cluster_id}.aln" "${cluster_id}.core.full.aln"
        fi

        if [[ -f "${cluster_id}.tab" ]]; then
            mv "${cluster_id}.tab" "${cluster_id}.core.tab"
        fi

        # Ensure output files exist
        if [[ ! -f "${cluster_id}.core.full.aln" ]]; then
            echo "Warning: Core alignment not generated, creating minimal alignment"
            echo ">${representative_id}" > ${cluster_id}.core.full.aln
            echo "N" >> ${cluster_id}.core.full.aln
        fi

        if [[ ! -f "${cluster_id}.core.tab" ]]; then
            touch ${cluster_id}.core.tab
        fi
    fi

    echo "Snippy alignment completed for cluster ${cluster_id}"
    echo "Output alignment size: \$(wc -c < ${cluster_id}.core.full.aln) bytes"
    echo "Number of sequences: \$(grep -c '^>' ${cluster_id}.core.full.aln)"

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    snippy: \$(snippy --version 2>&1 | head -n1 | sed 's/^/    /')
END_VERSIONS
EOF
"""
}