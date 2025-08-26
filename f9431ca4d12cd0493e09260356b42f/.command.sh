#!/bin/bash -euo pipefail
# Create input file for SKA
    cat > cluster_18_input.tsv <<'EOFSKA'
GCA_963562965_1	GCA_963562965_1.fasta
IP-0006-7_S6_L001-SPAdes	IP-0006-7_S6_L001-SPAdes.fasta
IP-0037-2_S7_L001-SPAdes	IP-0037-2_S7_L001-SPAdes.fasta
IP-0096-8_S8_L001-SPAdes	IP-0096-8_S8_L001-SPAdes.fasta
EOFSKA

    # Verify input file was created correctly
    echo "Input file contents:"
    cat cluster_18_input.tsv
    
    # Verify that all files exist
    echo "Checking file existence:"
    for file in *.fasta *.fa *.fna *.fas *.fsa; do
        if [ -f "$file" ]; then
            echo "Found: $file"
        fi
    done 2>/dev/null || echo "No FASTA files found with standard extensions"
    
    # Build SKA file
    ska build \
        -o cluster_18 \
         \
        cluster_18_input.tsv

    cat <<-END_VERSIONS > versions.yml
    "ASSEMBLY_SNPS_SCALABLE:CLUSTERED_SNP_TREE:SKA_BUILD":
        ska: $(ska --version 2>&1 | head -n1 | sed 's/^/    /')
    END_VERSIONS
