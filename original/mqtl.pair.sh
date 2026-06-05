#!/bin/bash
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=%j.out


QPGEN=/net/fantasia/home/hmkang/code/working/qpgen/bin/qpgentools
BASE=/net/topmed11/working/jblamer/qc.xqtl/match/metab.phase2.afr
LIST=/net/1000g/hmkang/topmed/pQTL/analysis/2025_12/pqtl/geno/topmed_freeze12c.minDP0.autosomes.pfiles.tsv
PAIRS=${BASE}/stats/TOPMed_mQTL_phase2_afr.sentinel.p1e06.merged.sorted.pair.format.tsv
OUT=${BASE}/pgs/metab.qtl.2026.phase2.sentinel.1e06.freeze12c

${QPGEN} pair-prs --pgen-list ${LIST} \
	--pairs ${PAIRS} \
	--out ${OUT}


