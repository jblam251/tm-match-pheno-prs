#!/usr/bin/bash

# args
BASE=$(pwd)
CSV=$1
OUT=${BASE}/$2.intnpx
# scripts
NPX2WIDE=/net/1000g/hmkang/topmed/pQTL/scripts/npx_sorted_to_wide.py
WIDE2QPGEN=/net/topmed11/working/jblamer/qc.xqtl/scripts/wide2qpgen.R
QPGEN2REGENIE=/net/1000g/hmkang/topmed/pQTL/scripts/topmed_reformat_pheno_for_regenie.py
# columns in the npx csv - the order must match: topID olinkID npx uniprot gene plateID wellID assayType
COLS=(1 6 14 7 8 4 3 9)

# get assay information
echo "extracting assay information"
CUTCOL_ASSAYS=${COLS[1]},${COLS[3]},${COLS[4]},${COLS[7]}
cat ${CSV} | grep -v ^Sample | grep -v ^SAMPLE | grep -v -w ext_ctrl | grep -v -w inc_ctrl | grep -v -w amp_ctrl | cut -d "," -f ${CUTCOL_ASSAYS} | tr , '\t' | sort | uniq > ${OUT}.assays.tsv

# get sample information
echo "extracting sample information"
CUTCOL_SAMPLES=${COLS[0]},${COLS[6]},${COLS[5]}
cat ${CSV} | grep -v ^Sample | grep -v ^SAMPLE | grep -v -w ext_ctrl | grep -v -w inc_ctrl | grep -v -w amp_ctrl | grep -v _CONTROL | cut -d "," -f ${CUTCOL_SAMPLES} | tr , '\t' | sort | uniq > ${OUT}.samples.tsv

# csv to long-format
echo "converting CSV to long-format NPX"
CUTCOL_NPX=${COLS[0]},${COLS[1]},${COLS[2]}
cat ${CSV} | grep -v ^Sample | grep -v ^SAMPLE | grep -v -w ext_ctrl | grep -v -w inc_ctrl | grep -v -w amp_ctrl | grep -v _CONTROL | cut -d "," -f ${CUTCOL_NPX} | tr , '\t' | grep -v -w NA | sort -k 1,2 | pigz -p32 -c > ${OUT}.npx.tsv.gz

# long-format to wide
echo "converting long-format NPX to wide-format"
python3 ${NPX2WIDE} --assays ${OUT}.assays.tsv --samples ${OUT}.samples.tsv --npx ${OUT}.npx.tsv.gz --out ${OUT}

# wide to qpgen
echo "converting wide-format npx to qpgen format"
Rscript ${WIDE2QPGEN} ${OUT}.npx.wide.tsv.gz ${OUT}.assays.tsv ${OUT}.qpgen.tsv
gzip ${OUT}.qpgen.tsv

# qpgen to regenie
echo "converting qpgen format to REGENIE format"
python3 ${QPGEN2REGENIE} --npx ${OUT}.qpgen.tsv.gz --out ${OUT}.regenie.tsv.gz
