#!/bin/bash -euo pipefail
set -euo pipefail

    echo "=== Gubbins Diagnostics for cluster cluster_2 ===" > .diagnostics.log
    seq_count=$(grep -c "^>" cluster_2.core.full.aln)
    aln_len=$(awk '/^[^>]/ {sum += length($0)} END {print sum}' cluster_2.core.full.aln)
    echo "Sequence count: $seq_count" >> .diagnostics.log
    echo "Alignment length: $aln_len" >> .diagnostics.log

    # Check if input files exist and are readable
    echo "Input file checks:" >> .diagnostics.log
    echo "Alignment file: cluster_2.core.full.aln - $(ls -la cluster_2.core.full.aln)" >> .diagnostics.log
    echo "Starting tree: cluster_2.treefile - $(ls -la cluster_2.treefile)" >> .diagnostics.log

    if [ $seq_count -lt 3 ]; then
        echo "WARNING: Alignment has only $seq_count sequences. Skipping Gubbins." >> .diagnostics.log
        touch cluster_2.filtered_polymorphic_sites.fasta
        touch cluster_2.recombination_predictions.gff
        touch cluster_2.node_labelled.final_tree.tre
        echo "Gubbins not run due to insufficient sequences." >> .diagnostics.log
    else
        echo "Running Gubbins with the following parameters:" >> .diagnostics.log
        echo "  iterations: 3" >> .diagnostics.log
        echo "  tree_builder: iqtree" >> .diagnostics.log
        echo "  first_tree_builder: rapidnj" >> .diagnostics.log
        echo "  min_snps: 5" >> .diagnostics.log
        echo "  use_hybrid: true" >> .diagnostics.log
        echo "  cpus: 8" >> .diagnostics.log

        if [ "true" = "true" ]; then
            echo "Running Gubbins with hybrid mode..." >> .diagnostics.log
            run_gubbins.py \
                --starting-tree cluster_2.treefile \
                --prefix cluster_2 \
                --first-tree-builder rapidnj \
                --tree-builder iqtree \
                --iterations 3 \
                --min-snps 5 \
                --threads 8 \
                cluster_2.core.full.aln >> .diagnostics.log 2>&1
            gubbins_exit_code=$?
        else
            echo "Running Gubbins without hybrid mode..." >> .diagnostics.log
            run_gubbins.py \
                --starting-tree cluster_2.treefile \
                --prefix cluster_2 \
                --tree-builder iqtree \
                --iterations 3 \
                --min-snps 5 \
                --threads 8 \
                cluster_2.core.full.aln >> .diagnostics.log 2>&1
            gubbins_exit_code=$?
        fi

        echo "Gubbins exit code: $gubbins_exit_code" >> .diagnostics.log

        if [ $gubbins_exit_code -ne 0 ]; then
            echo "ERROR: Gubbins failed for cluster cluster_2 with exit code $gubbins_exit_code" >> .diagnostics.log
            echo "Creating empty output files..." >> .diagnostics.log
            touch cluster_2.filtered_polymorphic_sites.fasta
            touch cluster_2.recombination_predictions.gff
            touch cluster_2.node_labelled.final_tree.tre
        fi

        # Check output file sizes and log
        for f in cluster_2.filtered_polymorphic_sites.fasta cluster_2.recombination_predictions.gff cluster_2.node_labelled.final_tree.tre; do
            if [ ! -s "$f" ]; then
                echo "WARNING: Output file $f is empty." >> .diagnostics.log
                touch "$f"
            else
                echo "Output file $f size: $(stat -c%s "$f") bytes" >> .diagnostics.log
            fi
        done
    fi

    cat <<-END_VERSIONS > versions.yml
"RECOMBINATION_AWARE_SNPS:GUBBINS_CLUSTER_DIAGNOSTIC":
    gubbins: $(run_gubbins.py --version | sed 's/^/    /')
END_VERSIONS
