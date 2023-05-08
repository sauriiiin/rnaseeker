
process DIFF_EXP {
    input:
        path outdir
        path sample_info
        val ref_sample
        path orfids
        path read_cnts

    output:
        path "$outdir/dea.csv"

    shell:
      """
      Rscript $baseDir/diff_exp.R $outdir/raw_counts.csv $outdir/dea.csv $sample_info $ref_sample $orfids
      """
}