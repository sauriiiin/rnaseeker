params {
  genome_fasta = "sacCer3.fa"
  gtf_file = "sacCer3.gtf"
  salmon_index = "sacCer3_salmon_index"
  r1_trim = 10
  r2_trim = 10
}

process generate_reads {
  input:
  set pair_id, file(fasta) from input_fasta_pairs
  
  output:
  set pair_id, file("R1.fastq.gz"), file("R2.fastq.gz") into reads_ch
  
  script:
  """
  Rscript generate_reads.R --fasta ${fasta} --pair_id ${pair_id} --output_dir reads
  """
}

process salmon_quant {
  input:
  set pair_id, file("R1.fastq.gz"), file("R2.fastq.gz") from reads_ch
  
  output:
  file("${pair_id}.tximport.tsv") into tximport_ch
  
  script:
  """
  salmon quant -i ${params.salmon_index} -l A -1 ${"R1.fastq.gz"} -2 ${"R2.fastq.gz"} -o salmon_output/${pair_id}
  """
}

process tximport {
  input:
  set pair_id, file(tximport_tsv) from tximport_ch
  
  output:
  file("${pair_id}_counts.tsv") into counts_ch
  
  script:
  """
  Rscript tximport.R --tximport_tsv ${tximport_tsv} --pair_id ${pair_id} --gtf_file ${params.gtf_file} --output_dir counts
  """
}

process differential_expression {
  input:
  set pair_id, file(counts_tsv) from counts_ch
  
  output:
  file("${pair_id}_DESeq2_results.tsv") into deseq2_ch
  
  script:
  """
  Rscript diff_exp.R --counts_tsv ${counts_tsv} --pair_id ${pair_id} --output_dir differential_expression
  """
}

process go_kegg {
  input:
  set pair_id, file(deseq2_tsv) from deseq2_ch
  
  output:
  file("${pair_id}_GO_enrichment.tsv") into go_kegg_ch
  
  script:
  """
  Rscript go_kegg.R --deseq2_tsv ${deseq2_tsv} --pair_id ${pair_id} --output_dir go_kegg
  """
}

workflow {
  input_fasta_pairs = Channel.fromPath("/input/*.fa").pair().map{ pair -> [ pair.left.baseName, pair.left, pair.right ] }
  
  generate_reads()
    .into(salmon_quant)
    .into(tximport)
    .into(differential_expression)
    .into(go_kegg)
}
