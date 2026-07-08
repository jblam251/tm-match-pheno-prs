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
          path("${prefix}.npx.wide.tsv"),
          path(omicsmap),
          path(pgs)

    script:
    """
    csv2plink.R \
        $pheno_input \
        $pheno_colmap \
        $prefix
    """
}


process GENIDMAP {

    tag "$prefix"
    
    input:
    tuple val(prefix),
          path(pheno_wide),
          path(omicsmap),
          path(pgs)
    
    output:
    tuple val(prefix),
          path(pheno_wide),
          path("${prefix}.idmap.tsv"),
          path(pgs)
    
    script:
    """
    Rscript ${params.scripts}/gen.idmap.R $omicsmap $pheno_wide ${prefix}.idmap.tsv
    """          
}

process SUBSET_PRS {

    tag "$prefix"
    
    input:
    tuple val(prefix),
          path(pheno_wide),
          path(idmap),
          path(pgs)
    
    output:
    tuple val(prefix),
          path(pheno_wide),
          path(idmap),
          path("${prefix}.pgs.tsv")
    
    script:
    """
    ${params.scripts}/subset.prs.sh $pgs $idmap ${prefix}.pgs.tsv
    """


}

process MATCH_PRS {

    tag "$prefix"
    
    input:
    tuple val(prefix),
          path(pheno_wide),
          path(idmap),
          path(pgs_subset)
    
    script:
    """
    /net/fantasia/home/hmkang/code/working/qpgen/bin/qpgentools match-prs-pheno \
	--pheno $pheno_wide \
	--prs $pgs_subset \
	--sample-tsv $idmap \
	--out ${prefix}.results \
        --threads $task.cpus \
	--lambda 1 \
	--mahalanobis
    """
}



workflow {
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

    MATCH_PRS(SUBSET_PRS(GENIDMAP(PHENO2PL(runx))))
}

