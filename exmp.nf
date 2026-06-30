#!/usr/bin/env nextflow

nextflow.enable.dsl=2

# run using
# $ nextflow run main.nf --input input.txt --run_optional true

#params.input = "input.txt"
#params.run_optional = true

Channel
    .fromPath(params.config)
    .splitCsv(header: true)
    .set { params_ch }

params_ch = Channel
    .fromPath(params.config)
    .splitCsv(header:true)
    .map { row ->
        tuple(
            row.prefix,
            file(row.csv),
            file(row.genotypes),
            file(row.pgs),
            file(row.mapping)
        )
    }

process csv2plink.format {

    tag "$prefix"

    input:
    tuple val(prefix), path(csv)

    output:
    tuple val(prefix), path("${prefix}.intnpx")

    script:
    """
    bash csv2plink.sh $csv $prefix
    """

}

#process SCRIPT2 {
#
#    input:
#    path infile
#
#    output:
#    path "b.txt"
#
#    script:
#    """
#    bash script2.sh $infile
#    """
#}

workflow {
    runc = Channel
        .fromPath(params.config)
        .splitCsv(header:true)
        .map { row ->
            tuple(
                row.prefix,
                file(row.csv),
                file(row.genotypes),
                file(row.pgs),
                file(row.mapping)
            )
        }

    step1_out = csv2plink.format(runc)

}
