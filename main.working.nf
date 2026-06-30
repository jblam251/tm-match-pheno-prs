#!/usr/bin/env nextflow

nextflow.enable.dsl=2

#PGS = params.pgs
#OMICSMAP = params.omicsmap
#THREADS = params.threads
#DISTSET = params.distset
#DISTLAMBDA = params.distlambda

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


#process PHENO2PL {
#
#    input:
#    tuple path(pheno_input), path(pheno_colmap), path(OMICSMAP), path(PGS), value(prefix)
#
#    output:
#    tuple path("assay.tsv"), path("samples.tsv"), path("npx_plink.txt"), path("idmap.txt"), path(PGS)
#
#    script:
#    """
#    csv2plink.sh \
#        $pheno_input \
#        $pheno_colmap \
#        $OMICSMAP \
#        $prefix
#    """
#}

process PHENO_PREPROCESS {

    input:
    tuple path(npx_plink),
          path(IDMAP),
          path(PGS)

    output:
    tuple path("npx_plink_preprocessed.txt")

    script:
    """
    npx.preprocess.sh \
        --input $npx_plink \
        --out npx_plink_preprocessed.txt
    """
}

process MATCH_PRS_PHENO {

    input:
    tuple path(npx_plink),
          path(npx_plink_preprocessed),
          path(IDMAP),
          path(PGS)

    script:
    """
    pqtl.match.sh \
        --pheno $preprocessed \
        --pgs ${prefix}.pgs \
        --mapping ${prefix}.map
    """
}


workflow {

    analyses = Channel
        .fromPath(params.settings)
        .splitCsv(header:true)
        .map { row ->
            tuple(
                row.prefix,
                file(row.prefix),
                file(row.pheno),
                file(row.colmap)
            )
        }

    step1 = CSV2REGENIE(analyses)

    step2 = NPX_PREPROCESS(step1)

    PQTL_MATCH(step2)
}

