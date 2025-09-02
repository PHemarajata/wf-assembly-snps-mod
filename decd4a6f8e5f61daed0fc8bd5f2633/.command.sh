#!/bin/bash -euo pipefail
bash <<'EOF'
    echo "Running Snippy alignment for cluster cluster_2"
    echo "Reference: /home/cdcadmin/wf-assembly-snps-mod/work/73/b0327fbf73b6cd89fc5ddc3daa9725/representative_id.txt"
    echo "Samples: [IE-0003_S8_L001-SPAdes, IE-0004_S9_L001-SPAdes, IE-0005_S10_L001-SPAdes, IE-0006_S11_L001-SPAdes, IE-0007_S1_L001-SPAdes, IE-0009_S3_L001-SPAdes, IE-0010_S4_L001-SPAdes, IE-0011_S5_L001-SPAdes, IE-0012_S6_L001-SPAdes, IE-0013_S7_L001-SPAdes, IE-0014_S8_L001-SPAdes, IE-0015_S9_L001-SPAdes, IE-0016_S10_L001-SPAdes, IE-0017_S11_L001-SPAdes, IE-0018_S12_L001-SPAdes, IE-0019_S5-SPAdes, IE-0021_S7-SPAdes, IE-0022_S8-SPAdes, IE-0024_S10-SPAdes, IE-0026_S8-SPAdes, IP-0001-8_S1_L001-SPAdes, IP-0002-6_S2_L001-SPAdes, IP-0003-4_S3_L001-SPAdes, IP-0004-2_S4_L001-SPAdes, IP-0005-9_S5_L001-SPAdes, IP-0006-7_S6_L001-SPAdes, IP-0007-5_S7_L001-SPAdes, IP-0008-3_S8_L001-SPAdes, IP-0009-1_S9_L001-SPAdes, IP-0010-9_S10_L001-SPAdes, IP-0011-7_S2_L001-SPAdes, IP-0012-5_S3_L001-SPAdes, IP-0013-3_S4_L001-SPAdes, IP-0014-1_S5_L001-SPAdes, IP-0015-8_S6_L001-SPAdes, IP-0016-6_S7_L001-SPAdes, IP-0017-4_S8_L001-SPAdes, IP-0018-2_S9_L001-SPAdes, IP-0019-0_S10_L001-SPAdes, IP-0020-8_S11_L001-SPAdes, IP-0021-6_S1_L001-SPAdes, IP-0022-4_S2_L001-SPAdes, IP-0023-2_S3_L001-SPAdes, IP-0024-0_S4_L001-SPAdes, IP-0025-7_S5_L001-SPAdes, IP-0026-5_S6_L001-SPAdes, IP-0027-3_S7_L001-SPAdes, IP-0028-1_S8_L001-SPAdes, IP-0029-9_S9_L001-SPAdes, IP-0030-7_S10_L001-SPAdes]"

    # Create working directory
    mkdir -p snippy_work

    for sample in $(echo [IE-0003_S8_L001-SPAdes, IE-0004_S9_L001-SPAdes, IE-0005_S10_L001-SPAdes, IE-0006_S11_L001-SPAdes, IE-0007_S1_L001-SPAdes, IE-0009_S3_L001-SPAdes, IE-0010_S4_L001-SPAdes, IE-0011_S5_L001-SPAdes, IE-0012_S6_L001-SPAdes, IE-0013_S7_L001-SPAdes, IE-0014_S8_L001-SPAdes, IE-0015_S9_L001-SPAdes, IE-0016_S10_L001-SPAdes, IE-0017_S11_L001-SPAdes, IE-0018_S12_L001-SPAdes, IE-0019_S5-SPAdes, IE-0021_S7-SPAdes, IE-0022_S8-SPAdes, IE-0024_S10-SPAdes, IE-0026_S8-SPAdes, IP-0001-8_S1_L001-SPAdes, IP-0002-6_S2_L001-SPAdes, IP-0003-4_S3_L001-SPAdes, IP-0004-2_S4_L001-SPAdes, IP-0005-9_S5_L001-SPAdes, IP-0006-7_S6_L001-SPAdes, IP-0007-5_S7_L001-SPAdes, IP-0008-3_S8_L001-SPAdes, IP-0009-1_S9_L001-SPAdes, IP-0010-9_S10_L001-SPAdes, IP-0011-7_S2_L001-SPAdes, IP-0012-5_S3_L001-SPAdes, IP-0013-3_S4_L001-SPAdes, IP-0014-1_S5_L001-SPAdes, IP-0015-8_S6_L001-SPAdes, IP-0016-6_S7_L001-SPAdes, IP-0017-4_S8_L001-SPAdes, IP-0018-2_S9_L001-SPAdes, IP-0019-0_S10_L001-SPAdes, IP-0020-8_S11_L001-SPAdes, IP-0021-6_S1_L001-SPAdes, IP-0022-4_S2_L001-SPAdes, IP-0023-2_S3_L001-SPAdes, IP-0024-0_S4_L001-SPAdes, IP-0025-7_S5_L001-SPAdes, IP-0026-5_S6_L001-SPAdes, IP-0027-3_S7_L001-SPAdes, IP-0028-1_S8_L001-SPAdes, IP-0029-9_S9_L001-SPAdes, IP-0030-7_S10_L001-SPAdes] | sed 's/[][]//g; s/,/ /g'); do
        echo "Processing sample: $sample"

        # Find the assembly file for this sample
        assembly_file=""
        for file in *.fa *.fasta *.fna; do
            if [[ "$file" == *"$sample"* ]] || [[ "$(basename "$file" .fa)" == "$sample" ]] || [[ "$(basename "$file" .fasta)" == "$sample" ]] || [[ "$(basename "$file" .fna)" == "$sample" ]]; then
                assembly_file="$file"
                break
            fi
        done

        if [[ -z "$assembly_file" ]]; then
            echo "Warning: Could not find assembly file for sample $sample"
            continue
        fi

        echo "Found assembly: $assembly_file"

        # Skip if this is the reference sample
        if [[ "$sample" == "/home/cdcadmin/wf-assembly-snps-mod/work/73/b0327fbf73b6cd89fc5ddc3daa9725/representative_id.txt" ]]; then
            echo "Skipping reference sample $sample"
            continue
        fi

        # Run snippy
        snippy             --outdir snippy_work/$sample             --ref IP-0030-7_S10_L001-SPAdes.fa             --ctgs $assembly_file             --cpus 16             --force              || {
            echo "Warning: Snippy failed for sample $sample"
            continue
        }
    done

    # Collect all snippy output directories
    snippy_dirs=()
    for dir in snippy_work/*/; do
        if [[ -d "$dir" ]] && [[ -f "$dir/snps.vcf" ]]; then
            snippy_dirs+=("$dir")
        fi
    done

    if [[ ${#snippy_dirs[@]} -eq 0 ]]; then
        echo "Warning: No successful Snippy runs found"
        # Create empty alignment
        echo ">/home/cdcadmin/wf-assembly-snps-mod/work/73/b0327fbf73b6cd89fc5ddc3daa9725/representative_id.txt" > cluster_2.core.full.aln
        echo "N" >> cluster_2.core.full.aln
        touch cluster_2.core.tab
    else
        echo "Found ${#snippy_dirs[@]} successful Snippy runs"

        # Run snippy-core to generate core alignment
        snippy-core             --ref IP-0030-7_S10_L001-SPAdes.fa             --prefix cluster_2             ${snippy_dirs[@]} || {
            echo "Warning: snippy-core failed, creating minimal alignment"
            echo ">/home/cdcadmin/wf-assembly-snps-mod/work/73/b0327fbf73b6cd89fc5ddc3daa9725/representative_id.txt" > cluster_2.core.full.aln
            echo "N" >> cluster_2.core.full.aln
            touch cluster_2.core.tab
        }

        # Ensure output files exist
        if [[ ! -f "cluster_2.core.full.aln" ]]; then
            echo "Warning: Core alignment not generated, creating minimal alignment"
            echo ">/home/cdcadmin/wf-assembly-snps-mod/work/73/b0327fbf73b6cd89fc5ddc3daa9725/representative_id.txt" > cluster_2.core.full.aln
            echo "N" >> cluster_2.core.full.aln
        fi

        if [[ ! -f "cluster_2.core.tab" ]]; then
            touch cluster_2.core.tab
        fi
    fi

    echo "Snippy alignment completed for cluster cluster_2"
    echo "Output alignment size: $(wc -c < cluster_2.core.full.aln) bytes"
    echo "Number of sequences: $(grep -c '^>' cluster_2.core.full.aln)"

    cat <<-END_VERSIONS > versions.yml
"RECOMBINATION_AWARE_SNPS:SNIPPY_ALIGN":
        snippy: $(snippy --version 2>&1 | head -n1 | sed 's/^/    /')
    END_VERSIONS
    EOF
