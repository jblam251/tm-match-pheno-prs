# TOPMed Multi-Omics Sample Identity Checking
Multi-omics data generated through the TOPMed program is evaluated for sample identity at the TOPMed Informatics Research Center (IRC).  The IRC leverages multiple genotyope-guided approaches to peform this task.  While these stratagies were developed for application to TOPmed data, some are assay- and technology-independent and can thus be deployed on arbitrary multi-omics studies with matched genotypes


## PGS-Based Sample Identity Checking
***RNAseq, Methylation, Metabolomics, Proteomics***
The IRC has developed a software, anQChor, that leverages published molecular quantitative trait locus (xQTL) studies for comprehensive QC of genotyped multi-omics datasets. anQChor uses summary statistics from cis- and trans-xQTL analyses to compute polygenic scores (PGS) for thousands of molecular traits.  Although individual PGS explain only small fractions of trait variance, their aggregated signal provides powerful QC metrics when combined across multiple individuals or multiple molecular traits
### How to Run
'''
nextflow run main.nf -c config.runx --type [data type] --settings [local parameters]
'''

The *type* argument can take one of four values: rnaseq, methylation, metabolomics, or proteomics.  The *--settings* argument is a CSV file with the first column specifying the desired prefix for the output files and the second column the full path location to the phenotype files.  For RNAseq, metabolomics, and proteomics, this must specify the location of the gene expression summary table, metabolite peak area table, and the proteomics NPX table respectively.  For methylation, this should be the LEVEL3 directory which contains the noob-adjusted beta values. Finally for metabolomics and proteomics, the settings file must contain a third column: for metabolomics this is the full path location to the chemical annotation file, whereas for proteomics this is a column mapping file which specifies which columns in the NPX data file correspond to the NPX, Sample ID, and Assay ID.  Prior to any run, a preliminary step must be taken to calculate PGS based on known QTLs (see below).  


## Pileup-Based Sample Identity Checking
***RNAseq***
In brief, the RNA sample identity checking follows the following procedure
1. Input data we need include (a) VCF file containing WGS-based genotypes, and (b) BAM files from RNA-seq reads
2. Genotype VCF file is converted into PGEN format using `plink2`, subsetted to relevant samples and variants within exons. 
3. Using the `bcdseq plp-dump` tool, create pileups from RNA-seq reads on the specified variant sites. Typical output file size could be ~5MB.
4. Using the `qpgentools match-plp-geno` tool, identify the best-matching sample IDs from the VCF with the pileups from RNA-seq

## Methylation Finger-print SNP Identity Checking (Metylation)
***Methylation***


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
