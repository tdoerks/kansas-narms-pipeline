#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AMR_ANALYSIS } from './workflows/amr_analysis'

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.input)
        .splitCsv(header:true)
        .map { row -> 
            def meta = [:]
            meta.id = row.sample
            meta.organism = row.organism
            return [meta, file(row.fasta)]
        }
        .set { ch_samples }

    AMR_ANALYSIS(ch_samples)
}
