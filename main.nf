/*
 * Proof of concept of a RNAseq pipeline implemented with Nextflow
 *
 * Author:
 * - Saurin Parikh <saurinbp@gmail.com>  
 */

/* 
 * enables modules 
 */
nextflow.enable.dsl = 2

/*
 * Default pipeline parameters.
 */
params.reads = "$baseDir/data/ggal/ggal_gut_{1,2}.fq"
params.sample_names = ["sample1", "sample2", "sample3", "sample4", "sample5", "sample6"]
params.replicate_info = ["1", "1", "1", "2", "2", "2"]
params.outdir = "results"

// import modules
include { QUANT } from './modules/quant/'

/* 
 * main script flow
 */
workflow {
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true ) 
    QUANT(read_pairs_ch)
}