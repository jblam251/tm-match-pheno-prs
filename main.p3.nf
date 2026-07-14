#!/usr/bin/env nextflow


process PHENO2PL_PROTEIN {

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
    ${params.scripts}/pheno2pl_protein.R\
        $pheno_input \
        $pheno_colmap \
        $prefix
    """
}

process PHENO2PL_METABOLITE {

    tag "$prefix"
    
    input:
    tuple val(prefix),
          path(pheno_input),
          path(metabol_annotation),
          path(traits),
          path(omicsmap),
          path(pgs)
    
    output:
    tuple val(prefix),
          path("${prefix}.peakareas.wide.tsv"),
          path(omicsmap),
          path(pgs)

    script:
    """
    Rscript ${params.scripts}/pheno2pl_metabolite.R \
	--peak_areas $pheno_input \
	--annotation $metabol_annotation \
	--traits $traits \
	--out $prefix
    """
}

process PHENO2PL_METHYLATION {

    tag "$prefix"
    
    input:
    tuple val(prefix),
          path(data_dir),
          path(traits),
          path(omicsmap),
          path(pgs)

    output:
    tuple val(prefix),
          path("${prefix}.beta.noob.tsv"),
          path(omicsmap),
          path(pgs)

    script:
    """
    ${params.scripts}/pheno2pl_methylation.sh \
        $data_dir \
        $prefix \
        $traits
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

    if (params.type == "proteomics") {
        runx = Channel
            .fromPath(params.settings)
            .splitCsv(header:true)
            .map { row ->
                tuple(
                    row.prefix,
                    file(row.pheno),
                    file(row.colmap),
                    file(params.omicsmap),
                    file(params.pgs)
                )
            }
        step1 = PHENO2PL_PROTEIN(runx)
    }

    else if (params.type == "metabolomics") {
        runx = Channel
            .fromPath(params.settings)
            .splitCsv(header:true)
            .map { row ->
                tuple(
                    row.prefix,
                    file(row.pheno),
                    file(row.metabol_annotation),
                    file(params.omicsmap),
                    file(params.pgs)
                )
            }
        step1 = PHENO2PL_METABOLITE(runx)
    }

    else if (params.type == "methylation") {
        runx = Channel
            .fromPath(params.settings)
            .splitCsv(header:true)
            .map { row ->
                tuple(
                    row.prefix,
                    file(row.pheno),
                    file(params.traits),
                    file(params.omicsmap),
                    file(params.pgs)
                )
            }
        step1 = PHENO2PL_METHYLATION(runx)
    }

    MATCH_PRS(SUBSET_PRS(GENIDMAP(step1)))
}

