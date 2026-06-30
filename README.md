# topmed-pqtl-id-check
jake blamer
05/29/2026

a few scripts for the pqtl match workflow

if plan is to use nextflow, then
 * 4-5 process blocks
 * connected by channels
 * with 1-2 conditional branches 

1) csv2regenie.sh*
	- input  : NPX csv, column mapping file
	- output : assay tab, sample tab, NPX in plink format

2) npx.preprocess.sh
	- input  : NPX in plink format (from step 1)
	- output : NPX files with imputation, PC adjustment, and RINT

3) pqtl.pair.sh (if PGS needs to be calculated)
	- input  : file of paths to genotypes split by chr, known QTL summary stats
	- output : merged PGS sample table

4) pqtl.match.sh
	- input  : NPX in plink format (from step 1) NPX preprocessed (from step 2), PGS, map file
	- output : qc match results

5) pqtl.match.res.analysis.R
	- input  : qc match results
	- output : post-process results/plots at several thresholds
