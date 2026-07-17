#!/usr/bin/bash

# USAGE
# $ ./beta.merge.sh [BETA DIR] [PREFIX] [PROBE LIST]

#BINDR=/net/topmed11/working/jblamer/id-checking/tm-match-pheno-prs/bin/bind.betas.R
BETADIR=$1
PREFIX=$2
TRAITS=$3
BINDR=$4

mkdir -p tmpd.batch
mkdir -p tmpd.batch/tmpd.stage
OUT=tmpd.batch

echo "[$(date +%T)] extracting beta noob files from LEVEL3 zipped directories"
for zip in $BETADIR/*LEVEL3*; do
	echo "[$(date +%T)] processing $zip"
	unzip -qq -j $zip "*beta.noob.csv" -d ${OUT}
	for file in ${OUT}/*.beta.noob.csv; do
		head -n1 $file  > ${OUT}/tmpd.stage/tmpf.01.tsv
		grep -f $TRAITS $file | sed -E 's/_.{4}//g' >> ${OUT}/tmpd.stage/tmpf.01.tsv
		mv ${OUT}/tmpd.stage/tmpf.01.tsv $file
	done
done


echo "[$(date +%T)] merging per-sample beta noob results"
Rscript $BINDR $OUT $PREFIX

echo "[$(date +%T)] pre-processing complete"
rm -r tmpd.batch

