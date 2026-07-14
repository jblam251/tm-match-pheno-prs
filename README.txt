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





