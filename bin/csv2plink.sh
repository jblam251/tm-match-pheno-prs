#!/usr/bin/bash

# scripts
NPX2WIDE=/net/topmed11/working/jblamer/id-checking/proteomics/bin/npx_sorted_to_plink.py
# input/output
CSV=/net/topmed11/incoming/topmed/proteomics/bidmc/2024.0903.FHS.OlinkHT/FHS_OlinkHT_11182024.csv
OUT=/net/topmed11/working/jblamer/id-checking/proteomics/out/tmpf.npx
# map columns
COLMAP=/net/topmed11/working/jblamer/id-checking/proteomics/npx.column.mapping.csv
mapfile -t COLS < <(cut -d"," -f2 $COLMAP | tail -n +2)


# get assay information
echo "extracting assay information"
CUTCOL_ASSAYS=${COLS[1]},${COLS[3]},${COLS[4]},${COLS[7]}
cat ${CSV} | grep -v ^Sample | grep -v ^SAMPLE | grep -v -w ext_ctrl | grep -v -w inc_ctrl | grep -v -w amp_ctrl | cut -d "," -f ${CUTCOL_ASSAYS} | tr , '\t' | sort | uniq > ${OUT}.assays.tsv

# get sample information
echo "extracting sample information"
CUTCOL_SAMPLES=${COLS[0]},${COLS[6]},${COLS[5]}
cat ${CSV} | grep -v ^Sample | grep -v ^SAMPLE | grep -v -w ext_ctrl | grep -v -w inc_ctrl | grep -v -w amp_ctrl | grep -v CONTROL | cut -d "," -f ${CUTCOL_SAMPLES} | tr , '\t' | sort | uniq > ${OUT}.samples.tsv

# csv to long-format
echo "converting CSV to long-format NPX"
CUTCOL_NPX=${COLS[0]},${COLS[1]},${COLS[2]}
cat ${CSV} | grep -v ^Sample | grep -v ^SAMPLE | grep -v -w ext_ctrl | grep -v -w inc_ctrl | grep -v -w amp_ctrl | grep -v _CONTROL | cut -d "," -f ${CUTCOL_NPX} | tr , '\t' | grep -v -w NA | sort -k 1,2 | pigz -p32 -c > ${OUT}.npx.tsv.gz

# long-format to wide plink format
echo "converting long-format NPX to wide-format"
python3 ${NPX2WIDE} --assays ${OUT}.assays.tsv --samples ${OUT}.samples.tsv --npx ${OUT}.npx.tsv.gz --out ${OUT}
