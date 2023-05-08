
process QUANT {
    input:
        path indir
        path index
        path trimmed
        
    output:
        path "${indir}/salmon_quant"

    shell:
        """
        rm -f $indir/foo.bar
        touch $indir/foo.bar
        for fn in $indir/rawdata/*.fastq;
        do
          samp=`basename \${fn}`
          dir=`dirname \${fn}`
          grep -qxF "\${samp:0:7}" $indir/foo.bar || echo "\${samp:0:7}" >> $indir/foo.bar
        done
        
        
        input="$indir/foo.bar"
        mkdir -p $indir/salmon_quant
        
        while IFS= read -r samp
        do
          salmon quant -i $index -p 20 --validateMappings --writeUnmappedNames -l A -1 $indir/trimmed/\${samp}_R1_001.trimmed.fastq/\${samp}_R1_001_trimmed.fq -2 $indir/trimmed/\${samp}_R2_001.trimmed.fastq/\${samp}_R2_001_trimmed.fq  -o "$indir/salmon_quant/\${samp}_quant"
        done < "\$input"
        
        """
}