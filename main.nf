/*
 * RNAseq pipeline implemented with Nextflow
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
params.indir = "$baseDir/fasta_files"
params.index = "$baseDir/index"
params.outdir = "$baseDir/results"
params.sample_info = "$baseDir/fasta_files/sample_info.csv"
params.ref_sample = "WT"
params.orfids = "$baseDir/translatome_orfid.csv"


// import modules
include { TRIM } from './modules/trim/'
include { QUANT } from './modules/quant/'
include { GENERATE_READS } from './modules/generate_reads/'
include { DIFF_EXP } from './modules/diff_exp/'
include { GO_KEGG } from './modules/go_kegg/'

/* 
 * main script flow
 */
workflow {
    TRIM(params.indir)
    QUANT(params.indir, params.index, TRIM.out)
    GENERATE_READS(params.indir, params.outdir, QUANT.out)
    DIFF_EXP(params.outdir, params.sample_info, params.ref_sample, params.orfids, GENERATE_READS.out)
    GO_KEGG(params.outdir, DIFF_EXP.out)
}