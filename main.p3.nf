#!/usr/bin/env nextflow

nextflow.enable.dsl=2


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
    ${params.scripts}/pheno2pl_protein.R \
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
        $traits \
        ${params.scripts}/bind.betas.R
    """
}

process PHENO2PL_RNASEQ {

    tag "$prefix"
    
    input:
    tuple val(prefix),
          path(pheno_input),
          path(omicsmap),
          path(pgs)

    output:
    tuple val(prefix),
          path("${prefix}.reads.wide.tsv"),
          path(omicsmap),
          path(pgs)

    script:
    """
    ${params.scripts}/pheno2pl_rnaseq.R \
        $pheno_input \
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
	--out ${params.type}.${prefix}.results \
        --threads $task.cpus \
	--lambda 1 \
	--mahalanobis
    """
}



workflow {

    if (params.type == "proteomics") {
            params.pgs = "/net/topmed11/working/jblamer/qc.xqtl/match/pqtl/pgs/ukb_ppp_v2_sentinel_GeneNames.topmed_freeze12c_minDP0_cis.prs.tsv.gz"
    }
    if (params.type == "methylation") {
            params.pgs = "/net/topmed11/working/jblamer/qc.xqtl/match/methqtl.0318/pgs/methyl.qtl.freeze12c.sentinel.p6e14.rank1.top10k.prs.tsv.gz"
            params.traits = "/net/topmed11/working/jblamer/qc.xqtl/match/methqtl.0318/traits/TOPMed_mQTL_freeze1_traits_p6e14_rank1.top10k.tsv"
    }
    if (params.type == "metabolomics") {
            params.pgs = "/net/topmed11/working/jblamer/qc.xqtl/match/metab.phase2/pgs/metab.qtl.2026.phase2.conditional.freeze12c.prs.tsv.gz"
            params.traits = "/net/topmed11/working/jblamer/qc.xqtl/match/metab.phase2/traits/metab.qtl.2026.phase2.traits.txt"
    }
    if (params.type == "rnaseq") {
            params.pgs = "/net/topmed11/working/jblamer/qc.xqtl/match/rna.bulk/pgs/gtex.v9.whole.blood.susie.tm.eqtl.freeze12c.prs.tsv.gz"
    }




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
    
    else if (params.type == "rnaseq") {
        runx = Channel
            .fromPath(params.settings)
            .splitCsv(header:true)
            .map { row ->
                tuple(
                    row.prefix,
                    file(row.pheno),
                    file(params.omicsmap),
                    file(params.pgs)
                )
            }
        step1 = PHENO2PL_RNASEQ(runx)
    }

    MATCH_PRS(SUBSET_PRS(GENIDMAP(step1)))
}

