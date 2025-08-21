//
// Pre-clustering using Mash distances
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Local modules
//
include { MASH_SKETCH      } from "../../modules/local/mash_sketch/main"
include { MASH_DIST        } from "../../modules/local/mash_dist/main"
include { CLUSTER_GENOMES  } from "../../modules/local/cluster_genomes/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN CLUSTERING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CLUSTERING {

    take:
    ch_assemblies // channel: [ val(sample_id), path(assembly) ]

    main:
    ch_versions = Channel.empty()

    // PROCESS: Create Mash sketches for each assembly
    MASH_SKETCH (
        ch_assemblies
    )
    ch_versions = ch_versions.mix(MASH_SKETCH.out.versions)

    // PROCESS: Calculate pairwise distances
    MASH_DIST (
        MASH_SKETCH.out.sketch.map{ sample_id, sketch -> sketch }.collect()
    )
    ch_versions = ch_versions.mix(MASH_DIST.out.versions)

    // PROCESS: Cluster genomes based on distances
    CLUSTER_GENOMES (
        MASH_DIST.out.distances
    )
    ch_versions = ch_versions.mix(CLUSTER_GENOMES.out.versions)

    // Create channel with cluster assignments
    ch_cluster_assignments = CLUSTER_GENOMES.out.clusters
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            tuple(row.cluster_id, row.sample_id)
        }

    // Join with original assemblies to create clustered channel
    ch_clustered_assemblies = ch_cluster_assignments
        .map { cluster_id, sample_id -> tuple(sample_id, cluster_id) }
        .join(ch_assemblies)
        .map { sample_id, cluster_id, assembly -> tuple(cluster_id, sample_id, assembly) }
        .groupTuple(by: 0)
        .branch { cluster_id, sample_ids, assemblies ->
            multi_sample: sample_ids.size() > 1
            singleton: sample_ids.size() == 1
        }

    // Handle singleton clusters based on parameter
    if (params.merge_singletons) {
        // Merge all singletons into one large cluster
        ch_merged_singletons = ch_clustered_assemblies.singleton
            .map { cluster_id, sample_ids, assemblies -> [sample_ids[0], assemblies[0]] }
            .collect()
            .map { singleton_pairs ->
                if (singleton_pairs.size() > 1) {
                    tuple("merged_singletons", singleton_pairs)
                } else {
                    // If only one singleton, skip it
                    null
                }
            }
            .filter { it != null }

        // Combine multi-sample clusters with merged singletons
        ch_clustered_assemblies = ch_clustered_assemblies.multi_sample
            .mix(ch_merged_singletons)
    } else {
        // Log and skip singleton clusters
        ch_clustered_assemblies.singleton
            .subscribe { cluster_id, sample_ids, assemblies ->
                log.info "Skipping singleton cluster ${cluster_id} (sample: ${sample_ids[0]}) - phylogenetic analysis requires multiple samples"
            }

        // Process only multi-sample clusters
        ch_clustered_assemblies = ch_clustered_assemblies.multi_sample
    }
        .map { cluster_id, sample_ids, assemblies ->
            def sample_assembly_pairs = [sample_ids, assemblies].transpose()
            tuple(cluster_id, sample_assembly_pairs)
        }

    // Count clusters for summary
    ch_clustered_assemblies
        .count()
        .subscribe { count ->
            if (count == 0) {
                log.warn "No multi-sample clusters found for phylogenetic analysis."
                log.warn "All samples appear to be singletons (each in its own cluster)."
                log.warn "Consider:"
                log.warn "  - Adjusting --mash_threshold (current: ${params.mash_threshold ?: 0.03}) to allow more similar samples to cluster together"
                log.warn "  - Using --merge_singletons to combine all singletons into one large cluster"
                log.warn "Phylogenetic analysis requires at least 2 samples per cluster."
            } else {
                log.info "Found ${count} clusters for phylogenetic analysis"
            }
        }

    emit:
    versions             = ch_versions
    clusters             = CLUSTER_GENOMES.out.clusters
    cluster_summary      = CLUSTER_GENOMES.out.summary
    clustered_assemblies = ch_clustered_assemblies // channel: [ val(cluster_id), [ [sample_id, assembly], ... ] ]
}