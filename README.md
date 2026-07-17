TOPMed Multi-Omics Sample Identity Checking
Summer 2026

Sample identity checking for omics data generated through the TOPMed program is performed at the TOPMed IRC.  While all sample checks are executed via genotype-matching with the qpgentools software, the stratagy depends on the omics type: RNAseq data relies on pipeups from RNAseq reads on a specified set of variants to match to genotypes.  Methylation, Metabolomics, and Proteomics data uses published QTLs to calculate genotype-based scores, which are then compared to measured assay values to determine sample identity. Details on each approach can be found below

PGS-Based Sample Identity Checking (Methylation, Metabolomics, Proteomics)
leverages published molecular quantitative trait locus (xQTL) studies for comprehensive QC of genotyped multi-omics datasets. anQChor uses summary statistics from cis- and trans-xQTL analyses to compute polygenic scores (PGS) for thousands of molecular traits, including gene expression and protein abundance. Although individual PGS explain only small fractions of trait variance, their aggregated signal provides robust and powerful QC metrics when combined across multiple individuals or multiple molecular traits.

Pileup-Based Sample Identity Checking (RNAseq)
In brief, the RNA sample identity checking follows the following procedure
1. Input data we need include (a) VCF file containing WGS-based genotypes, and (b) BAM files from RNA-seq reads
2. Genotype VCF file is converted into PGEN format using `plink2`, subsetted to relevant samples and variants within exons. 
3. Using the `bcdseq plp-dump` tool, create pileups from RNA-seq reads on the specified variant sites. Typical output file size could be ~5MB.
4. Using the `qpgentools match-plp-geno` tool, identify the best-matching sample IDs from the VCF with the pileups from RNA-seq



### MISC NOTES TO SELF ###

qc identity checking from genotyped multi-omics data using 
anqchor workflow.  three similar but slightly different 
pipelines: protein match, metabolite match, and methylation 
match. while the exact scripts and inputs/outputs may vary 
depending on data type, the general workflow steps are 
generally the same: 
	1) raw phenotype data reformat to plink format
	2) generate 2 column idmap file (NWD-to-TO*)
	3) subset full freeze12 PRS file
	4) run prs-pheno match
these steps are preceeded by tasks which are not currently
included in the nextflow pipeline: 1) QTL summary statistic 
generation, 2) QTL summary statistic re-formatting, and 3)
calculation of genotype-anchored PRS for the TOPMed freeze 12
samples. there will also probably be some downstream analysis
and aggregation which hasn't been configured yet

TO DO
- implement import of covariate file
- (X) fix detection of duplicate TO* IDs in phenotypes
- (X) fix idmap issue for CARDIA study
- rename scripts and variables for consistancy
- re-organize the three pipelines into one (somehow)

omics "type" goes in config file, assign in workflow block
workflow {

    if (params.type == "olink")
        step1 = PHENO2PL_OLINK(samples)

    else if (params.type == "metabolon")
        step1 = PHENO2PL_METABOLON(samples)

    else if (params.type == "epic")
        step1 = PHENO2PL_EPIC(samples)

    step2 = GENIDMAP(step1)
    step3 = SUBSET_PRS(step2)
    MATCH_PRS(step3)
}
