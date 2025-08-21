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
        .combine(ch_assemblies, by: 1)
        .map { sample_id, cluster_id, assembly ->
            tuple(cluster_id, sample_id, assembly)
        }
        .groupTuple(by: 0)
        .map { cluster_id, sample_ids, assemblies ->
            def sample_assembly_pairs = [sample_ids, assemblies].transpose()
            tuple(cluster_id, sample_assembly_pairs)
        }

    emit:
    versions             = ch_versions
    clusters             = CLUSTER_GENOMES.out.clusters
    cluster_summary      = CLUSTER_GENOMES.out.summary
    clustered_assemblies = ch_clustered_assemblies // channel: [ val(cluster_id), [ [sample_id, assembly], ... ] ]
}