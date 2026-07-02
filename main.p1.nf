#!/usr/bin/env nextflow



process PHENO2PL {

    tag "$prefix"
    
    input:
    tuple val(prefix),
          path(pheno_input),
          path(pheno_colmap),
          path(omicsmap),
          path(pgs)

    output:
    tuple val(prefix),
          path("${prefix}.assay.tsv"),
          path("${prefix}.samples.tsv"),
          path("${prefix}.npx_plink.txt"),
          path("${prefix}.idmap.txt"),
          path(pgs)

    script:
    """
    csv2plink.sh \
        $pheno_input \
        $pheno_colmap \
        $omicsmap \
        $prefix
    """
}


workflow {
    println "CONFIG CHECK: ${params}"

    runx = Channel
        .fromPath(params.settings)
        .splitCsv(header:true)
        .view()
        .map { row ->
            tuple(
                row.prefix,
                file(row.pheno),
                file(row.colmap),
                file(params.omicsmap),
                file(params.pgs)
            )
        }

    PHENO2PL(runx)
}

