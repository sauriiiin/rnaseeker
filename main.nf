params.input_dir = "/fasta_files"
params.sample_names = ["sample1", "sample2", "sample3", "sample4", "sample5", "sample6"]
params.replicate_info = ["1", "1", "1", "2", "2", "2"]

input:
  set val(input_files), file("${params.input_dir}/*.fasta") into input_fasta_files

process map {
    input:
        set val(input), file(fasta_file) from input_fasta_files

    output:
        file("${params.input_dir}/${fasta_file.baseName}_salmon")

    shell:
        "salmon quant -i /index -p 20 --validateMappings --writeUnmappedNames -l A -r ${input} -o ${output}"
}

process trim {
    input:
        set val(input), file(salmon_output_file) from "${params.input_dir}/*_salmon/quant.sf"
    
    output:
        file("${params.input_dir}/${salmon_output_file.baseName}_trimmed")

    shell:
        "trim_galore -q 20 --length 50 --fastqc --illumina --stringency 3 ${input} -o ${output}"
}

process generate_reads {
    input:
        set val(input), file(trimmed_file) from "${params.input_dir}/*_trimmed"

    output:
        file("${params.input_dir}/${trimmed_file.baseName}_reads.csv")

    script:
        """
        Rscript generate_reads.R ${input} ${output} ${params.sample_names} ${params.replicate_info}
        """
}

process diff_exp {
    input:
        set val(input), file(reads_file) from "${params.input_dir}/*_reads.csv"

    output:
        file("${params.input_dir}/${reads_file.baseName}_diff_exp.csv")

    script:
        """
        Rscript diff_exp.R ${input} ${output} ${params.sample_names} ${params.replicate_info}
        """
}

process go_kegg {
    input:
        set val(input), file(diff_exp_file) from "${params.input_dir}/*_diff_exp.csv"

    output:
        file("${params.input_dir}/${diff_exp_file.baseName}_go_kegg.csv")

    script:
        """
        Rscript go_kegg.R ${input} ${output}
        """
}

output:
    set file("${params.input_dir}/*_trimmed/*"),
        file("${params.input_dir}/*_reads.csv"),
        file("${params.input_dir}/*_diff_exp.csv"),
        file("${params.input_dir}/*_go_kegg.csv") into result_files
    result_files.collect().each { f -> f.baseName = f.baseName.replaceAll('_.*', '') }
    file("${params.input_dir}") into output_dir
