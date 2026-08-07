# Sample Identity Checking TOPMed Multi-Omics Data
This pipeline automates the execution of a recently developed identity check method that leverages summary statistics from published cis- and trans- molecular quantitative trait locus (xQTL) studies to compute polygenic scores (PGS) for thousands of molecular traits.  PGS is then compared to measured trait values from multi-omics data to evaluate sample assignment.


## Pipeline workflow
![Alt Text](images/pgs.schematic.1.png)


## How to run
```
nextflow run main.nf -c config.runx --settings [sample settings CSV]
```

## Results
Identity check results are written to a file [prefix].match.all.tsv.gz which contains the identifier for the omics sample, the identifier for the assigned corresponding genotype, and the z-score between them.  It also includes the genotype identifiers and z-scores for the top five genotype matches 

Based on the above results, each sample is assigned a status in the `MatchStatus` column which can be interpreted as follows: if the top match for the omics sample is its assigned corresponding genotype, it is labeled `SELF_BEST`.  If it is not the top match, it may be considered a lenient match (`SELF_LENIENT`) or a non-match (`UNCLEAR`) depending on whether the z-score is greater or less than the lenient match threshold defined by the `--z-threshold` argument.  Omics samples that fail to match their assigned genotype but exhibit a strong match to a different genotype are flagged as `SINGLE_NEW_BEST` (or `MULTI_NEW_BEST` if it matches more than one non-assigned genotype).  If an omics sample does not have an assigned genotype, the status is automatically set to `SELF_BEST` (but may be alternatively tagged as `NO_ASSIGNED_GT` for diagnostic plots)

![Alt Text](images/diagnostic.plot.1.png)



## Creating a Settings file
This is a comma-separated (CSV) file containing one row per dataset and 8 columns (see below).  Not all columns are required for each omics data type, this is handled automatically in the workflow.  This file can include multiple rows for batch processing

| Column | Required | Description |
|--------|----------|-------------|
| `type` | Yes | The omics data type specification. It must take one of the following: `rnaseq`, `methylation`, `metabolomics`, or `proteomics`. |
| `prefix` | Yes | The prefix used for naming output files. |
| `pheno` | Yes | The input file of molecular phenotypes. For RNAseq, metabolomics, and proteomics, this must specify the location of the gene expression summary table, metabolite peak area table, and the proteomics NPX table respectively.  For methylation, this should be the LEVEL3 directory which contains the noob-adjusted beta values.|
| `pgs` |  Yes | The genotype-derived polygenic scores for each molecular trait.  If scores have yet to be generated, `qpgentools pair-prs` can calculate PGS when provided genotypes and a set of known QTL summary statistics.  See `pair-prs` below for more detail. |
| `omicsmap` | Yes | A file for mapping genotype identifiers to omics identifiers for all tested participants.  Genotype identifiers must appear in a column named `NWD_ID` while omics identifiers in column `SAMPLE_ID`. The cross-cohort omics sample attributes files available on the TOPMed combined exchange area fulfills these requirements.| 
| `traits` | Metabolomics, Methylation | A single-column file of molecular trait labels. |
| `metabolite_annotation` | Metabolomics | The metabolite annotation file which was provided during data generation.  This is sometimes called the Chemical annotation file. |
| `protein_colmap` | Proteomics | A column mapping file which specifies which columns in the NPX data file correspond to the NPX, Sample ID, and Assay ID. |

## Calculating PGS using `pair-prs`
Polygenic scores must be in PLINK format with the genotype identifiers in the first two columns (FID, IID) and the molecular traits in the remaining columns.  It's important to ensure the trait labels in the PGS file match those found in the molecular data set.  While calculating PGS can be done in any number of ways, one method to accomplish this is `pair-prs` from the qpgentools software.  It requires two arguments: the genotypes and the xQTL summary statistics.  Genotypes are passed using a "pairs" file which contains the paths to the per-chromosome VCF files.  The xQTL summary statistics must contain the trait, variant, beta, standard error, and log10 p-value. With the genoypes and correctly-formatted summary statistics, PGS can be calculated using

```
qpgentools pair-prs \
        --pgen-list genotype.pfiles.tsv \
        --pairs xqtl.summary.stats.tsv \
        --out pgs.tsv
```


## Example input files
```text
$ head -n6 pheno.tsv | cut -f1-5
Name	rna_1	rna_2	rna_3	rna_4
ENSG00000268903	7	5	26	20
ENSG00000241860	316	389	147	222
ENSG00000308579	0	0	0	0
ENSG00000278267	0	0	1	2
ENSG00000248901	0	0	0	0

$ head -n6 pgs.tsv | cut -f1-5
FID	IID	ENSG00000268903	ENSG00000241860	ENSG00000308579	
geno_a	geno_a	1.1592	-0.8913	-0.1258	
geno_b	geno_b	0.0995	-0.3332	1.1983	
geno_c	geno_c	-0.0892	-1.239	-0.1258	
geno_d	geno_d	-0.0995	-0.8913	1.1983	
geno_e	geno_e	1.1592	-0.3332	-0.1258	

$ head -n6 omicsmap.tsv
NWD_ID	SAMPLE_ID
geno_a	rna_1
geno_a	rna_2
geno_c	rna_3
geno_d	rna_4
geno_e	rna_5

$ head -n6 traits.txt
cg00002190
cg00002646
cg00004883
cg00006735
cg00008795
cg00010187

$ head -n6 metabolite_annotation.tsv | cut -f1-5
CHEM_ID	LIB_ID	COMP_ID	SUPER_PATHWAY	SUB_PATHWAY	
35	400	42370	Amino Acid	Glutamate Metabolism
50	400	485	Amino Acid	Polyamine Metabolism
55	400	27665	Cofactors and Vitamins	Nicotinate Metabolism
62	209	38395	Lipid	Fatty Acid, Dihydroxy	
93	305	528	Energy	TCA Cycle	

$ cat protein_colmap.csv
Field,ColumnNumber
TOP_ID,1
GeneName,6
NPX,17

$ head genotype.pfiles.tsv
1       1       1000000000      /path/to/genos.chr1
2       1       1000000000      /path/to/genos.chr2
3       1       1000000000      /path/to/genos.chr3
4       1       1000000000      /path/to/genos.chr4
5       1       1000000000      /path/to/genos.chr5

$ head -n6 xqtl.summary.stats.tsv
trait	variant		beta	se	log10p
NOC2L	1:953778:G:C	-0.7512	0.0321	115.00
KLHL17	1:959193:G:A	0.6864	0.0223	193.21
HES4	1:1000112:G:T	-0.8752	0.0178	445.97
AGRN	1:1010481:T:A	0.6347	0.0283	106.22
ISG15	1:1013490:C:G	1.3448	0.0273	447.27

```

