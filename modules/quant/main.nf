params.index = "./index"

process QUANT {
    tag "$pair_id"

    input:
        tuple val(pair_id), path(reads)

    output:
        path pair_id

    shell:
        """
        salmon quant -i ${params.index} -p 20 --validateMappings --writeUnmappedNames -l A -1 ${reads[0]} -2 ${reads[1]} -o salmon_quants/$pair_id
        """
}