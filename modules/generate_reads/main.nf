
process GENERATE_READS {
    input:
        path indir
        path outdir
        path salmon_quant

    output:
        path "$outdir/raw_counts.csv"

    shell:
      """
      mkdir -p $baseDir/results
      Rscript $baseDir/generate_reads.R $indir $outdir
      """
}