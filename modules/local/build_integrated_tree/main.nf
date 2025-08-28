process BUILD_INTEGRATED_TREE {
    tag "integrated_phylogeny"
    label 'process_high'
    container "quay.io/biocontainers/iqtree:2.2.6--h21ec9f0_0"
    
    publishDir "${params.outdir}/Integrated_Results", mode: params.publish_dir_mode, pattern: "*.{treefile,iqtree,log}"

    input:
    path integrated_alignment
    path sample_mapping

    output:
    path "integrated_core_snps.treefile", emit: tree
    path "integrated_core_snps.iqtree", emit: log
    path "integrated_phylogeny_report.txt", emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def model = params.iqtree_model ?: 'GTR+ASC'
    """
    # Check if alignment has sequences
    seq_count=\$(grep -c "^>" ${integrated_alignment} || echo "0")
    
    if [ \$seq_count -lt 3 ]; then
        echo "WARNING: Integrated alignment has only \$seq_count sequences. Cannot build phylogenetic tree."
        
        # Create empty output files
        touch integrated_core_snps.treefile
        touch integrated_core_snps.iqtree
        
        echo "Insufficient sequences for phylogenetic analysis" > integrated_phylogeny_report.txt
        
    else
        echo "Building phylogenetic tree from integrated core SNPs alignment (\$seq_count sequences)"
        
        # Build tree with IQ-TREE
        iqtree2 \\
            -s ${integrated_alignment} \\
            -st DNA \\
            -m MFP \\
            -bb 1000 \\
            -alrt 1000 \\
            -nt AUTO \\
            --prefix integrated_core_snps \\
            || {
            echo "WARNING: IQ-TREE failed. Creating empty output files."
            touch integrated_core_snps.treefile
            touch integrated_core_snps.iqtree
        }
        
        # Create phylogeny report
        # Install compatible versions of numpy and pandas
        pip install --upgrade numpy>=1.15.4
        pip install pandas || echo "pandas installation failed, continuing without it"
        
        python3 << 'EOF'
try:
    import pandas as pd
    pandas_available = True
except ImportError:
    pandas_available = False
    print("Warning: pandas not available, creating basic report")

import os

# Read sample mapping
try:
    if pandas_available:
        mapping_df = pd.read_csv("${sample_mapping}", sep='\\t')
    else:
        # Basic file reading without pandas
        with open("${sample_mapping}", 'r') as f:
            lines = f.readlines()
        mapping_data = []
        for line in lines[1:]:  # Skip header
            parts = line.strip().split('\\t')
            if len(parts) >= 2:
                mapping_data.append({'sample_id': parts[0], 'cluster_id': parts[1]})
    
    with open("integrated_phylogeny_report.txt", 'w') as f:
        f.write("INTEGRATED PHYLOGENETIC ANALYSIS REPORT\\n")
        f.write("=" * 50 + "\\n\\n")
        
        if pandas_available:
            f.write(f"Total samples in phylogeny: {len(mapping_df)}\\n")
            f.write(f"Clusters represented: {mapping_df['cluster_id'].nunique()}\\n\\n")
            
            # Cluster representation
            cluster_counts = mapping_df['cluster_id'].value_counts()
            f.write("Samples per cluster in integrated tree:\\n")
            f.write("-" * 40 + "\\n")
            for cluster_id, count in cluster_counts.items():
                f.write(f"Cluster {cluster_id}: {count} samples\\n")
            
            # SNP statistics
            if 'total_snps' in mapping_df.columns:
                avg_snps = mapping_df['total_snps'].mean()
                max_snps = mapping_df['total_snps'].max()
                min_snps = mapping_df['total_snps'].min()
                
                f.write(f"\\nSNP statistics:\\n")
                f.write(f"Average SNPs per sample: {avg_snps:.2f}\\n")
                f.write(f"Maximum SNPs per sample: {max_snps}\\n")
                f.write(f"Minimum SNPs per sample: {min_snps}\\n")
        else:
            f.write(f"Total samples in phylogeny: {len(mapping_data)}\\n")
            clusters = set([item['cluster_id'] for item in mapping_data])
            f.write(f"Clusters represented: {len(clusters)}\\n\\n")
        
        f.write(f"\\nPhylogenetic method: IQ-TREE with model selection\\n")
        f.write(f"Bootstrap support: 1000 replicates\\n")
        f.write(f"SH-aLRT support: 1000 replicates\\n")
        
        # Check if tree was successfully built
        if os.path.exists("integrated_core_snps.treefile") and os.path.getsize("integrated_core_snps.treefile") > 0:
            f.write(f"\\nTree construction: SUCCESSFUL\\n")
        else:
            f.write(f"\\nTree construction: FAILED\\n")

except Exception as e:
    with open("integrated_phylogeny_report.txt", 'w') as f:
        f.write(f"Error creating phylogeny report: {e}\\n")
EOF
    fi

    # Ensure all output files exist
    for file in integrated_core_snps.treefile integrated_core_snps.iqtree integrated_phylogeny_report.txt; do
        if [ ! -f "\$file" ]; then
            touch "\$file"
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(iqtree2 --version 2>&1 | head -n1 | sed 's/^/    /')
        pandas: \$(python -c "try: import pandas; print(pandas.__version__); except: print('not available')")
        numpy: \$(python -c "try: import numpy; print(numpy.__version__); except: print('not available')")
    END_VERSIONS
    """
}