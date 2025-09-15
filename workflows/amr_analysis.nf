#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process AMRFINDER_RUN {
    container 'quay.io/biocontainers/ncbi-amrfinderplus:3.12.8--h283d18e_0'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_amr.tsv"), emit: results
    tuple val(meta), path("${meta.id}_mutations.tsv"), emit: mutations, optional: true

    script:
    def organism_flag = meta.organism ? "-O ${meta.organism}" : ""
    """
    amrfinder \\
        -n ${fasta} \\
        ${organism_flag} \\
        --plus \\
        --threads ${task.cpus} \\
        -d /fastscratch/tylerdoe/ks_samples_project/shared_amrfinder_db/latest \\
        -o ${meta.id}_amr.tsv \\
        --mutation_all ${meta.id}_mutations.tsv
    """
}

workflow AMR_ANALYSIS {
    take:
    samples_ch

    main:
    // Run AMRFinder on each sample with absolute database path
    AMRFINDER_RUN(samples_ch)

    emit:
    results = AMRFINDER_RUN.out.results
    mutations = AMRFINDER_RUN.out.mutations
}
