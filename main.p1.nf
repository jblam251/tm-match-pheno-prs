#!/usr/bin/env nextflow



process PHENO2PL {

    input:
    tuple val(prefix),
          path(pheno_input),
          path(pheno_colmap),
          path(omicsmap),
          path(pgs)

    output:
    tuple val(prefix),
          path("assay.tsv"),
          path("samples.tsv"),
          path("npx_plink.txt"),
          path("idmap.txt"),
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
        .map { row ->
            tuple(
                row.prefix,
                file(row.pheno),
                file(row.colmap)
            )
        }

    PHENO2PL(runx)
}

