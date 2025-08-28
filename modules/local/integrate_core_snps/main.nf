process INTEGRATE_CORE_SNPS {
    tag "global_integration"
    label 'process_medium'
    container "quay.io/biocontainers/python:3.9--1"
    
    publishDir "${params.outdir}/Integrated_Results", mode: params.publish_dir_mode, pattern: "*.{fa,tsv,txt}"

    input:
    path core_snp_files
    path snp_position_files
    path clusters_file

    output:
    path "integrated_core_snps.fa", emit: integrated_alignment
    path "integrated_snp_positions.tsv", emit: integrated_positions
    path "core_snp_summary.txt", emit: summary
    path "sample_cluster_mapping.tsv", emit: sample_mapping
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Install compatible versions of numpy and pandas
    pip install --upgrade numpy>=1.15.4
    pip install pandas biopython

    python3 << 'EOF'
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict
import glob

def integrate_core_snps():
    # Integrate core SNPs from all clusters into a global alignment
    
    # Read cluster assignments
    clusters_df = pd.read_csv("${clusters_file}", sep='\\t')
    sample_to_cluster = dict(zip(clusters_df['sample_id'], clusters_df['cluster_id']))
    
    # Collect all core SNP sequences
    all_sequences = {}
    cluster_snp_counts = {}
    
    # Process each core SNP file
    core_snp_files = glob.glob("*_core_snps.fa")
    
    for snp_file in core_snp_files:
        cluster_id = snp_file.replace("_core_snps.fa", "")
        snp_count = 0
        
        try:
            for record in SeqIO.parse(snp_file, "fasta"):
                sample_id = record.id
                snp_sequence = str(record.seq)
                
                if sample_id not in all_sequences:
                    all_sequences[sample_id] = []
                
                all_sequences[sample_id].append({
                    'cluster': cluster_id,
                    'sequence': snp_sequence,
                    'length': len(snp_sequence)
                })
                snp_count = len(snp_sequence)
            
            cluster_snp_counts[cluster_id] = snp_count
            print(f"Processed cluster {cluster_id}: {snp_count} SNPs")
            
        except Exception as e:
            print(f"Warning: Could not process {snp_file}: {e}")
            cluster_snp_counts[cluster_id] = 0
    
    # Create integrated alignment by concatenating SNPs from all clusters
    integrated_sequences = {}
    
    for sample_id, cluster_data in all_sequences.items():
        # Sort by cluster to ensure consistent order
        cluster_data.sort(key=lambda x: x['cluster'])
        
        # Concatenate sequences from all clusters
        concatenated_seq = ''.join([data['sequence'] for data in cluster_data])
        integrated_sequences[sample_id] = concatenated_seq
    
    # Write integrated core SNPs alignment
    with open("integrated_core_snps.fa", 'w') as f:
        for sample_id, sequence in integrated_sequences.items():
            if sequence:  # Only write non-empty sequences
                f.write(f">{sample_id}\\n{sequence}\\n")
    
    # Integrate SNP position information
    all_positions = []
    position_offset = 0
    
    snp_position_files = glob.glob("*_snp_positions.tsv")
    
    for pos_file in snp_position_files:
        try:
            cluster_positions = pd.read_csv(pos_file, sep='\\t')
            if not cluster_positions.empty:
                # Adjust positions by offset
                cluster_positions['global_position'] = cluster_positions['position'] + position_offset
                cluster_positions['original_position'] = cluster_positions['position']
                all_positions.append(cluster_positions)
                
                # Update offset for next cluster
                position_offset += cluster_positions['position'].max()
                
        except Exception as e:
            print(f"Warning: Could not process {pos_file}: {e}")
    
    # Combine all position data
    if all_positions:
        integrated_positions = pd.concat(all_positions, ignore_index=True)
        integrated_positions.to_csv("integrated_snp_positions.tsv", sep='\\t', index=False)
    else:
        # Create empty file
        pd.DataFrame(columns=['position', 'ref_base', 'alt_bases', 'cluster_id', 'global_position', 'original_position']).to_csv("integrated_snp_positions.tsv", sep='\\t', index=False)
    
    # Create sample-cluster mapping
    sample_mapping = []
    for sample_id in integrated_sequences.keys():
        cluster_id = sample_to_cluster.get(sample_id, 'unknown')
        sample_mapping.append({
            'sample_id': sample_id,
            'cluster_id': cluster_id,
            'total_snps': len(integrated_sequences[sample_id])
        })
    
    mapping_df = pd.DataFrame(sample_mapping)
    mapping_df.to_csv("sample_cluster_mapping.tsv", sep='\\t', index=False)
    
    # Create summary report
    with open("core_snp_summary.txt", 'w') as f:
        f.write("INTEGRATED CORE SNP ANALYSIS SUMMARY\\n")
        f.write("=" * 50 + "\\n\\n")
        
        f.write(f"Total samples: {len(integrated_sequences)}\\n")
        f.write(f"Total clusters processed: {len(cluster_snp_counts)}\\n")
        f.write(f"Total integrated SNP positions: {sum(cluster_snp_counts.values())}\\n\\n")
        
        f.write("SNPs per cluster:\\n")
        f.write("-" * 20 + "\\n")
        for cluster_id, count in sorted(cluster_snp_counts.items()):
            f.write(f"Cluster {cluster_id}: {count} SNPs\\n")
        
        f.write(f"\\nIntegrated alignment length: {max([len(seq) for seq in integrated_sequences.values()]) if integrated_sequences else 0}\\n")
        
        if integrated_sequences:
            avg_snps = sum([len(seq) for seq in integrated_sequences.values()]) / len(integrated_sequences)
            f.write(f"Average SNPs per sample: {avg_snps:.2f}\\n")
    
    print(f"Integration complete: {len(integrated_sequences)} samples, {sum(cluster_snp_counts.values())} total SNP positions")

# Run integration
integrate_core_snps()
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}