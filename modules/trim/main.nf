
process TRIM {
    input:
        path fasta_files
    
    output:
        path "${fasta_files}/trimmed"
    
    shell:
        """
        for fn in ${fasta_files}/rawdata/*.fastq;
        do
          samp=`basename \${fn}`
          dir=`dirname \${fn}`
          mkdir -p ${fasta_files}/trimmed
          
          trim_galore -q 20 --fastqc --stringency 3 ${fasta_files}/rawdata/\${samp} -o ${fasta_files}/trimmed/\${samp:0:-6}.trimmed.fastq
        done
        """
}
