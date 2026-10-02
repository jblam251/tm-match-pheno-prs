#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process PHENO2PL_PROTEIN {

    tag "$prefix"
    
    input:
    tuple val(type),
          val(prefix),
          path(pheno_input),
          path(pgs),
          path(omicsmap),
          path(pheno_colmap)

    output:
    tuple val(type),
          val(prefix),
          path("${prefix}.wide.tsv"),
          path(omicsmap),
          path(pgs)

    script:
    """
    pheno2pl_protein.R \
        $pheno_input \
        $pheno_colmap \
        $prefix
    """
}

process PHENO2PL_METABOLITE {

    tag "$prefix"
    
    input:
    tuple val(type),
          val(prefix),
          path(pheno_input),
          path(pgs),
          path(omicsmap),
          path(traits),
          path(metabol_annotation)
    
    output:
    tuple val(type),
          val(prefix),
          path("${prefix}.wide.tsv"),
          path(omicsmap),
          path(pgs)

    script:
    """
    pheno2pl_metabolite.R \
	$pheno_input \
	$metabol_annotation \
	$traits \
	$prefix
    """
}

process PHENO2PL_METHYLATION {

    tag "$prefix"
    
    input:
    tuple val(type),
          val(prefix),
          path(data_dir),
          path(pgs),
          path(omicsmap),
          path(traits)

    output:
    tuple val(type),
          val(prefix),
          path("${prefix}.wide.tsv"),
          path(omicsmap),
          path(pgs)

    script:
    """
    pheno2pl_methylation.sh \
        $data_dir \
        $prefix \
        $traits \
        bind.betas.R
    """
}

process PHENO2PL_RNASEQ {

    tag "$prefix"
    
    input:
    tuple val(type),
          val(prefix),
          path(pheno_input),
          path(pgs),
          path(omicsmap)

    output:
    tuple val(type),
          val(prefix),
          path("${prefix}.wide.tsv"),
          path(omicsmap),
          path(pgs)

    script:
    """
    pheno2pl_rnaseq.R \
        $pheno_input \
        $prefix
    """
}

process GENIDMAP {

    tag "$prefix"
    
    input:
    tuple val(type),
          val(prefix),
          path(pheno_wide),
          path(omicsmap),
          path(pgs)
    
    output:
    tuple val(type),
          val(prefix),
          path(pheno_wide),
          path("${prefix}.idmap.tsv"),
          path(pgs)
    
    script:
    """
    gen.idmap.R $omicsmap $pheno_wide ${prefix}.idmap.tsv
    """          
}

process SUBSET_PRS {

    tag "$prefix"
    
    input:
    tuple val(type),
          val(prefix),
          path(pheno_wide),
          path(idmap),
          path(pgs)
    
    output:
    tuple val(type),
          val(prefix),
          path(pheno_wide),
          path(idmap),
          path("${prefix}.pgs.tsv")
    
    script:
    """
    subset.prs.sh $pgs $idmap ${prefix}.pgs.tsv
    """


}

process MATCH_PRS {

    tag "$prefix"
    
    input:
    tuple val(type),
          val(prefix),
          path(pheno_wide),
          path(idmap),
          path(pgs_subset)

    output:
    tuple val(type),
          val(prefix),
          path("${type}.${prefix}.results.match.all.tsv.gz")
    
    
    script:
    """
    qpgentools match-prs-pheno \
	--pheno $pheno_wide \
	--prs $pgs_subset \
	--sample-tsv $idmap \
	--out ${type}.${prefix}.results \
        --threads $task.cpus \
	--z-threshold 1.96 \
        --rint \
	--lambda 1 \
	--mahalanobis
    """
}


process PLOT_DIAGNOSTIC {

    tag "$prefix"
    
    input:
    tuple val(type),
          val(prefix),
          path(result)

    script:
    """
    plot.diagnostics.R $result
    """          

}

workflow {


    samples = Channel
        .fromPath(params.settings)
        .splitCsv(header: true)
    
    protein = samples
        .filter { it.type == "proteomics" }
        .map { row ->
            tuple(
                row.type,
                row.prefix,
                file(row.pheno),
                file(row.pgs),
                file(row.omicsmap),
                file(row.protein_colmap)
            )
        }
    s1_protein = PHENO2PL_PROTEIN(protein)

    metabolite = samples
        .filter { it.type == "metabolomics" }
        .map { row ->
            tuple(
                row.type,
                row.prefix,
                file(row.pheno),
                file(row.pgs),
                file(row.omicsmap),
                file(row.traits),
                file(row.metabolite_annotation)
            )
        }
    s1_metabolite = PHENO2PL_METABOLITE(metabolite)


    methyl = samples
        .filter { it.type == "methylation" }
        .map { row ->
            tuple(
                row.type,
                row.prefix,
                file(row.pheno),
                file(row.pgs),
                file(row.omicsmap),
                file(row.traits)
            )
        }
    s1_methyl = PHENO2PL_METHYLATION(methyl)


    rna = samples
        .filter { it.type == "rnaseq" }
        .map { row ->
            tuple(
                row.type,
                row.prefix,
                file(row.pheno),
                file(row.pgs),
                file(row.omicsmap)
            )
        }
    s1_rna = PHENO2PL_RNASEQ(rna)
    
    step1 = s1_protein
        .mix(s1_metabolite)
        .mix(s1_methyl)
        .mix(s1_rna)
    
    PLOT_DIAGNOSTIC(MATCH_PRS(SUBSET_PRS(GENIDMAP(step1))))

}

