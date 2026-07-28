# Sample Identity Checking TOPMed Multi-Omics Data
Multi-omics data generated through the TOPMed program is evaluated for sample identity at the TOPMed Informatics Research Center (IRC).  The IRC utilizes multiple approaches for this task. The information below describes a recently developed method that leverages summary statistics from published cis- and trans- molecular quantitative trait locus (xQTL) studies to compute polygenic scores (PGS) for thousands of molecular traits.  Although individual PGS explain only small fractions of trait variance, their aggregated signal provides useful QC metrics when combined across multiple individuals. While this pipeline was developed for application to TOPmed data, this stratagy is assay- and technology-independent and can thus be deployed on arbitrary multi-omics studies with matched genotypes


## Pipeline workflow
![Alt Text](images/pgs.schematic.1.png)


## How to run
```
nextflow run main.nf -c config.runx --settings [sample settings CSV]
```

## Results
to do

## Creating a Settings file
A settings file is required to set run-specific parameters for analysis.  This is a comma-separated (CSV) file containing one row per dataset and 8 columns (see below).  Not all columns are required for each omics data type, this is handled automatically in the workflow.  This file can include multiple rows for batch processing

| Column | Required? | Description |
|--------|-----------|-------------|
| `type` | Yes | The omics data type specification. It must take one of the following: `rnaseq`, `methylation`, `metabolomics`, or `proteomics`. |
| `prefix` | Yes | Prefix used for naming output files. |
| `pheno` | Yes | Input file of molecular phenotypes. For RNAseq, metabolomics, and proteomics, this must specify the location of the gene expression summary table, metabolite peak area table, and the proteomics NPX table respectively.  For methylation, this should be the LEVEL3 directory which contains the noob-adjusted beta values.|
| `pgs` |  Yes | Genotype-derrived polygenic scores for each molecular trait.  If scores have yet to be generated, `qpgentools pair-prs` can calculate PGS when provided genotypes and a set of known QTL summary statistics.  See `pair-prs` below for more detail. |
| `traits` | Metabolomics, Methylation | Single-column file of molecular trait labels. |
| `metabolite_annotation` | Metabolomics | The metabolite annotation file which was provided during data generation.  This is sometimes called the Chemical annotation file. |
| `protein_colmap` | Proteomics | Column mapping file which specifies which columns in the NPX data file correspond to the NPX, Sample ID, and Assay ID. |

```text
$ head -n6 traits.txt
cg00002190
cg00002646
cg00004883
cg00006735
cg00008795
cg00010187

$ head -n6 chemical.annotation.tsv | cut -f1-5
CHEM_ID	LIB_ID	COMP_ID	SUPER_PATHWAY	SUB_PATHWAY	
35	400	42370	Amino Acid	Glutamate Metabolism
50	400	485	Amino Acid	Polyamine Metabolism
55	400	27665	Cofactors and Vitamins	Nicotinate Metabolism
62	209	38395	Lipid	Fatty Acid, Dihydroxy	
93	305	528	Energy	TCA Cycle	

$ cat protein.column.map.csv
Field,ColumnNumber
TOP_ID,1
GeneName,6
NPX,17
```

## Calculating PGS using `pair-prs`
