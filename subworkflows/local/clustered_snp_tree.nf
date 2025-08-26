//
// Per-cluster SNP analysis and tree building
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Local modules
//
include { SKA_BUILD        } from "../../modules/local/ska_build/main"
include { SKA_ALIGN        } from "../../modules/local/ska_align/main"
include { IQTREE_FAST      } from "../../modules/local/iqtree_fast/main"
include { GUBBINS_CLUSTER  } from "../../modules/local/gubbins_cluster/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN CLUSTERED SNP TREE WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CLUSTERED_SNP_TREE {

    take:
    ch_clustered_assemblies // channel: [ val(cluster_id), val(sample_ids), path(assemblies) ]

    main:
    ch_versions = Channel.empty()

    // PROCESS: Build SKA files for each cluster
    SKA_BUILD (
        ch_clustered_assemblies
    )
    ch_versions = ch_versions.mix(SKA_BUILD.out.versions)

    // PROCESS: Create alignments from SKA files
    SKA_ALIGN (
        SKA_BUILD.out.ska_file
    )
    ch_versions = ch_versions.mix(SKA_ALIGN.out.versions)

    // PROCESS: Build fast ML trees
    IQTREE_FAST (
        SKA_ALIGN.out.alignment
    )
    ch_versions = ch_versions.mix(IQTREE_FAST.out.versions)

    // PROCESS: Run Gubbins on clusters if requested
    if (params.run_gubbins) {
        // Combine alignment and tree for Gubbins input
        ch_gubbins_input = SKA_ALIGN.out.alignment
            .join(IQTREE_FAST.out.tree, by: 0)

        GUBBINS_CLUSTER (
            ch_gubbins_input
        )
        ch_versions = ch_versions.mix(GUBBINS_CLUSTER.out.versions)

        // Use Gubbins filtered alignment for downstream analysis
        ch_final_alignments = GUBBINS_CLUSTER.out.filtered_alignment
        ch_final_trees = GUBBINS_CLUSTER.out.final_tree
        ch_recombination_gff = GUBBINS_CLUSTER.out.recombination_gff
    } else {
        // Use original SKA alignment and IQ-TREE trees
        ch_final_alignments = SKA_ALIGN.out.alignment
        ch_final_trees = IQTREE_FAST.out.tree
        ch_recombination_gff = Channel.empty()
    }

    emit:
    versions           = ch_versions
    ska_files          = SKA_BUILD.out.ska_file
    alignments         = ch_final_alignments
    trees              = ch_final_trees
    recombination_gff  = ch_recombination_gff
}