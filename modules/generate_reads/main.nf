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


params.input_dir = "/fasta_files"
params.sample_names = ["sample1", "sample2", "sample3", "sample4", "sample5", "sample6"]
params.replicate_info = ["1", "1", "1", "2", "2", "2"]

input_files = files("${params.input_dir}/*.fastq.gz")

process read_map {
    input:
        tuple val(input), file("${params.input_dir}/*.fastq.gz")

    output:
        tuple val(input), file("${params.input_dir}/${input}_salmon")

    shell:
        """
        salmon quant -i /index -p 20 --validateMappings --writeUnmappedNames -l A -r ${input} -o ${output}
        """
}

workflow {
    read_map(input_files)
}
