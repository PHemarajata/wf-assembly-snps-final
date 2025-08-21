#!/usr/bin/env python3
"""
Cluster genomes based on Mash distances using single-linkage clustering.
"""

import argparse
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import sys

def parse_mash_output(mash_file):
    """Parse Mash distance output file."""
    df = pd.read_csv(mash_file, sep='\t', header=None, 
                     names=['ref', 'query', 'distance', 'p_value', 'shared_hashes'])
    return df

def build_distance_matrix(df, threshold=0.03):
    """Build a sparse adjacency matrix from Mash distances."""
    # Get unique sample names
    samples = sorted(set(df['ref'].tolist() + df['query'].tolist()))
    sample_to_idx = {sample: idx for idx, sample in enumerate(samples)}
    
    # Filter by threshold and build adjacency matrix
    filtered_df = df[df['distance'] <= threshold]
    
    rows = []
    cols = []
    
    for _, row in filtered_df.iterrows():
        ref_idx = sample_to_idx[row['ref']]
        query_idx = sample_to_idx[row['query']]
        
        # Add both directions for undirected graph
        rows.extend([ref_idx, query_idx])
        cols.extend([query_idx, ref_idx])
    
    # Create sparse adjacency matrix
    data = [1] * len(rows)
    adj_matrix = csr_matrix((data, (rows, cols)), shape=(len(samples), len(samples)))
    
    return adj_matrix, samples

def cluster_samples(adj_matrix, samples):
    """Perform connected components clustering."""
    n_components, labels = connected_components(adj_matrix, directed=False)
    
    # Create cluster assignments
    clusters = {}
    for idx, label in enumerate(labels):
        sample = samples[idx]
        if label not in clusters:
            clusters[label] = []
        clusters[label].append(sample)
    
    return clusters

def write_clusters(clusters, output_file):
    """Write cluster assignments to TSV file."""
    with open(output_file, 'w') as f:
        f.write("cluster_id\tsample_id\n")
        for cluster_id, samples in clusters.items():
            for sample in samples:
                f.write(f"cluster_{cluster_id}\t{sample}\n")

def main():
    parser = argparse.ArgumentParser(description='Cluster genomes based on Mash distances')
    parser.add_argument('mash_file', help='Mash distance output file')
    parser.add_argument('output_file', help='Output cluster assignments file')
    parser.add_argument('--threshold', type=float, default=0.03, 
                       help='Distance threshold for clustering (default: 0.03)')
    parser.add_argument('--max-cluster-size', type=int, default=100,
                       help='Maximum cluster size (default: 100)')
    
    args = parser.parse_args()
    
    # Parse Mash output
    print(f"Reading Mash distances from {args.mash_file}")
    df = parse_mash_output(args.mash_file)
    
    # Build distance matrix and cluster
    print(f"Clustering with threshold {args.threshold}")
    adj_matrix, samples = build_distance_matrix(df, args.threshold)
    clusters = cluster_samples(adj_matrix, samples)
    
    # Split large clusters if needed
    final_clusters = {}
    cluster_counter = 0
    
    for cluster_id, samples_in_cluster in clusters.items():
        if len(samples_in_cluster) <= args.max_cluster_size:
            final_clusters[cluster_counter] = samples_in_cluster
            cluster_counter += 1
        else:
            # Split large clusters into smaller chunks
            chunk_size = args.max_cluster_size
            for i in range(0, len(samples_in_cluster), chunk_size):
                chunk = samples_in_cluster[i:i + chunk_size]
                final_clusters[cluster_counter] = chunk
                cluster_counter += 1
    
    # Write results
    print(f"Found {len(final_clusters)} clusters")
    for cluster_id, samples_in_cluster in final_clusters.items():
        print(f"  Cluster {cluster_id}: {len(samples_in_cluster)} samples")
    
    write_clusters(final_clusters, args.output_file)
    print(f"Cluster assignments written to {args.output_file}")

if __name__ == '__main__':
    main()