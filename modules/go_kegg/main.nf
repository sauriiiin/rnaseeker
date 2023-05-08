
process GO_KEGG {
    input:
        path outdir
        path dea

    output:
        path "$outdir/goe.csv"

    shell:
      """
      Rscript $baseDir/go_kegg.R $outdir/dea.csv $outdir/goe.csv
      """
}