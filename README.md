# PGS-Based Sample Identity Checking TOPMed Multi-Omics Data
Multi-omics data generated through the TOPMed program is evaluated for sample identity at the TOPMed Informatics Research Center (IRC).  The IRC leverages summary statistics from published cis- and trans- molecular quantitative trait locus (xQTL) studies to compute polygenic scores (PGS) for thousands of molecular traits.  Although individual PGS explain only small fractions of trait variance, their aggregated signal provides useful QC metrics when combined across multiple individuals. While this pipeline was developed for application to TOPmed data, this stratagy is assay- and technology-independent and can thus be deployed on arbitrary multi-omics studies with matched genotypes


## Pipeline workflow
![Alt Text](images/pgs.schematic.1.png)


## How to run
```
nextflow run main.nf -c config.runx --type [data type] --settings [local parameters]
```

```
| Argument | Required | Description |
|----------|----------|-------------|
| `--type` | Yes | Omics data type to process. Accepted values: `rnseq`,  `proteomics`, `methylation`, `metabolomics`. |
| `--settings` | Yes | Path to the sample settings CSV containing one row per dataset to process. |
```

## Results

## Creating a Settings file
The first column specifying the desired prefix for output files.  The second column should contain the full path location to the molecular phenotype files.  For RNAseq, metabolomics, and proteomics, this must specify the location of the gene expression summary table, metabolite peak area table, and the proteomics NPX table respectively.  For methylation, this should be the LEVEL3 directory which contains the noob-adjusted beta values. Finally for metabolomics and proteomics, the settings file must contain a third column: for metabolomics this is the full path location to the chemical annotation file, whereas for proteomics this is a column mapping file which specifies which columns in the NPX data file correspond to the NPX, Sample ID, and Assay ID.  Prior to any run, a preliminary step must be taken to calculate PGS based on known QTLs (see below).  

