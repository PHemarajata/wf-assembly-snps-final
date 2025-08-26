#!/bin/bash -euo pipefail
# Create input file for SKA
    cat > cluster_13_input.tsv <<'EOFSKA'
ERS013364	ERS013364.fasta
IP-0087-7_S10_L001-SPAdes	IP-0087-7_S10_L001-SPAdes.fasta
IP-0103-2_S3_L001-SPAdes	IP-0103-2_S3_L001-SPAdes.fasta
IP-0134-7_S4_L001-SPAdes	IP-0134-7_S4_L001-SPAdes.fasta
EOFSKA

    # Verify input file was created correctly
    echo "Input file contents:"
    cat cluster_13_input.tsv
    
    # Build SKA file
    ska build \
        -o cluster_13 \
         \
        cluster_13_input.tsv

    cat <<-END_VERSIONS > versions.yml
    "ASSEMBLY_SNPS_SCALABLE:CLUSTERED_SNP_TREE:SKA_BUILD":
        ska: $(ska --version 2>&1 | head -n1 | sed 's/^/    /')
    END_VERSIONS
