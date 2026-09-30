# PRS-Based Sample Identity Matching for Multi-Omics Data

This Nextflow pipeline automates an identity-check method that leverages summary
statistics from published cis- and trans-molecular quantitative trait locus
(xQTL) studies to compute polygenic scores (PGS) for thousands of molecular
traits. These PGS are then compared against measured trait values from
multi-omics data to evaluate whether omics samples are correctly assigned to
their corresponding genotype data.

The method was originally developed and validated using data from the
[TOPMed](https://topmed.nhlbi.nih.gov/) program, but it is generalizable to
any study with paired genotype and multi-omics data (RNA-seq, DNA methylation,
metabolomics, or proteomics).

## Table of Contents
- [Overview](#prs-based-sample-identity-matching-for-multi-omics-data)
- [Citation](#citation)
- [Requirements](#requirements)
- [Quick Start](#quick-start)
- [Pipeline Workflow](#pipeline-workflow)
- [Configuration](#configuration)
- [Creating a Settings File](#creating-a-settings-file)
- [Calculating PGS Using `pair-prs`](#calculating-pgs-using-pair-prs)
- [Example Input Files](#example-input-files)
- [Running the Pipeline](#running-the-pipeline)
- [Understanding the Results](#understanding-the-results)
- [Troubleshooting](#troubleshooting)
- [License](#license)
- [Support](#support)

---

## Citation

If you use this pipeline in published work, please cite:

> [Author names, year, title, journal/preprint venue, DOI]
> *(Add citation details here prior to public release.)*

Please also cite the underlying `qpgentools` software (see [Requirements](#requirements)).

---

## Requirements

| Dependency | Notes |
|---|---|
| [Nextflow](https://www.nextflow.io/) | Tested with version `>=XX.XX.X` — *(fill in tested version)* |
| Java | Required by Nextflow; version `>=11` recommended |
| R | Version `>=X.X` — *(fill in tested version)*, with the following packages installed: `data.table`, `dplyr`, *(add full list)* |
| [`qpgentools`](#) | Required for PGS calculation (`pair-prs`) and identity matching (`match-prs-pheno`). *(Add public repository link, installation instructions, and license here.)* |
| PLINK-format genotype data | Per-chromosome VCF or PLINK2 pgen/pvar/psam files, depending on your `pair-prs` input |

> **Note:** This pipeline calls `qpgentools` as an external binary. Ensure it
> is installed and available on your system, and set its location using the
> `--qpgentools` parameter (see [Configuration](#configuration)) or by adding
> it to your `PATH`. This pipeline does **not** bundle `qpgentools`.

### Compute Requirements
Resource needs scale with the number of samples and molecular traits being
tested. As a general guideline:
- *(Fill in typical memory/CPU/runtime for a dataset of N samples)*
- The pipeline can be run on a local workstation, an HPC cluster (via Slurm,
  SGE, etc.), or in a containerized environment (Docker/Singularity), depending
  on the Nextflow profile/config you provide.

---

## Quick Start

A minimal test dataset is provided in `test_data/` to verify your installation
works end-to-end. *(Add a small example dataset and settings CSV if planning
public release — this significantly lowers the barrier to adoption.)*

```bash
git clone https://github.com/jblam251/tm-match-pheno-prs.git
cd tm-match-pheno-prs
nextflow run main.nf -c config.runx --settings test_data/settings.example.csv
```

## Pipeline Workflow
The pipeline:
![Alt Text](images/pgs.schematic.1.png)

1. Converts omics-specific molecular phenotype files (RNA-seq, methylation,metabolomics, or proteomics) into a common "wide" format.
2. Generates a genotype-to-omics identifier map (GENIDMAP).
3. Subsets the polygenic score file to the relevant traits and samples(SUBSET_PRS).
4. Compares polygenic scores to observed molecular phenotypes to assesssample identity (MATCH_PRS).

## Configuration


## Creating a Settings File
This is a comma-separated (CSV) file containing one row per dataset and 8 columns (see below).  Not all columns are required for each omics data type, this is handled automatically in the workflow.  This file can include multiple rows for batch processing

| Column | Required | Description |
|--------|----------|-------------|
| `type` | Yes | The omics data type specification. It must take one of the following: `rnaseq`, `methylation`, `metabolomics`, or `proteomics`. |
| `prefix` | Yes | The prefix used for naming output files. |
| `pheno` | Yes | The input file of molecular phenotypes. For RNAseq, metabolomics, and proteomics, this must specify the location of the gene expression summary table, metabolite peak area table, and the proteomics NPX table respectively.  For methylation, this should be the LEVEL3 directory which contains the noob-adjusted beta values.|
| `pgs` |  Yes | The genotype-derived polygenic scores for each molecular trait.  If scores have yet to be generated, `qpgentools pair-prs` can calculate PGS when provided genotypes and a set of known QTL summary statistics.  See `pair-prs` below for more detail. |
| `omicsmap` | Yes | A file for mapping genotype identifiers to omics identifiers for all tested participants.  Genotype identifiers must appear in a column named `NWD_ID` while omics identifiers in column `SAMPLE_ID`. | 
| `traits` | Metabolomics, Methylation | A single-column file of molecular trait labels. |
| `metabolite_annotation` | Metabolomics | The metabolite annotation file which was provided during data generation.  This is sometimes called the Chemical annotation file. |
| `protein_colmap` | Proteomics | A column mapping file which specifies which columns in the NPX data file correspond to the NPX, Sample ID, and Assay ID. |



## Calculating PGS Using `pair-prs`
Polygenic scores must be in PLINK format with the genotype identifiers in the first two columns (FID, IID) and the molecular traits in the remaining columns.  It's important to ensure the trait labels in the PGS file match those found in the molecular data set.

While calculating PGS can be done in any number of ways, one method to accomplish this is `pair-prs` from the `qpgentools` software.  It requires two inputs: 
1. Genotypes — passed as a "pairs" file listing paths to per-chromosome VCF files. See the example format for the expected columns.
2. xQTL summary statistics — must contain columns for trait, variant,beta, standard error, and log10 p-value.

```
qpgentools pair-prs \
        --pgen-list genotype.pfiles.tsv \
        --pairs xqtl.summary.stats.tsv \
        --out pgs.tsv
```


## Example Input Files
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



## Running the Pipeline
```
nextflow run main.nf -c config.runx --settings [sample settings CSV]
```
* -c config.runx — a custom Nextflow config file specifying your executionenvironment (see Configuration).
* --settings — path to your settings CSV file.


## Understanding the Result

### Output File

Identity check results are written to a compressed tab-separated file named
[prefix].match.all.tsv.gz

This file contains one row per omics sample and includes the following information:

| Column | Description |
|---|---|
| `SAMPLE_ID` | Identifier for the omics sample. |
| `NWD_ID` (assigned) | Identifier for the genotype originally assigned to this omics sample (per the `omicsmap` input). |
| `Z_SELF` | Z-score between the omics sample and its assigned genotype. |
| `TOP1_ID` ... `TOP5_ID` | Genotype identifiers for the top five matching genotypes, ranked by z-score. |
| `TOP1_Z` ... `TOP5_Z` | Corresponding z-scores for each of the top five matches. |
| `MatchStatus` | Categorical assignment describing the identity match outcome (see below). |

> *(Confirm exact column names/order against your actual output — update this
> table to match precisely, including whether ranks 1–5 are inclusive of the
> assigned genotype or reported separately.)*

### Interpreting `MatchStatus`

Each omics sample is assigned one of the following statuses, based on
whether its top genotype match corresponds to its assigned genotype (from
`omicsmap`), and how strong that match is relative to a lenient-match
z-score threshold:

| Status | Condition |
|---|---|
| `SELF_BEST` | The assigned genotype **is** the top match (highest z-score) for the omics sample. |
| `SELF_LENIENT` | The assigned genotype is **not** the top match, but its z-score still exceeds the lenient match threshold. |
| `UNCLEAR` | The assigned genotype is **not** the top match, and its z-score falls below the lenient match threshold. |
| `SINGLE_NEW_BEST` | The omics sample fails to match its assigned genotype, but shows a strong match (above threshold) to exactly **one** other, non-assigned genotype. |
| `MULTI_NEW_BEST` | The omics sample fails to match its assigned genotype, but shows a strong match to **more than one** other, non-assigned genotype. |
| `NO_ASSIGNED_GT` | The omics sample has no assigned genotype in `omicsmap`. By default this is reported as `SELF_BEST` in the output, but is separately labeled `NO_ASSIGNED_GT` when generating diagnostic plots, so these samples can be visually distinguished. |


### The Lenient Match Threshold

The z-score cutoff used to distinguish `SELF_LENIENT`/strong matches from
`UNCLEAR` results is referred to throughout this document as the **lenient 
match threshold**. This parameter is not yet configurable 

Note: The lenient match threshold referenced above is intended to beconfigurable (e.g., via a --z-threshold argument). (Confirm where thisthreshold is applied in the current pipeline code — it does not currentlyappear as an argument to MATCH_PRS in main.nf and may need to be addedor clarified.)


### Diagnostic Plots

![Alt Text](images/diagnostic.plot.1.png)

Diagnostic plots visualize the distribution of z-scores across samples and
highlight `MatchStatus` categories, making it easier to spot systematic
mismatches (e.g., sample swaps, mislabeling, or batch effects).





## Troubleshooting


## License


## Support


