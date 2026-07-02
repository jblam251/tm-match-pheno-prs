#!/bin/bash
#SBATCH --mem=48G
#SBATCH --cpus-per-task=12
#SBATCH --time=01:00:00
#SBATCH --output=%j.out

#when comparing PGS to full freeze12, set to 8 threads! 

STUDY=ltrc
QPGEN=/net/fantasia/home/hmkang/code/working/qpgen/bin/qpgentools
BASE=/net/topmed11/working/jblamer/qc.xqtl
PRS=${BASE}/match/pgs/ukb_ppp_v2_sentinel_HT.topmed_freeze12c_minDP0_cis.prs.tsv.gz
NPX=${BASE}/data/${STUDY}.intnpx.pc${PC}.regenie.tsv.gz
IDMAP=${BASE}/match/idmaps/pqtl.combined.idmap.cardia.unmatched.nwds.rm.mesa.fails.added.tsv
SAMPLES=${BASE}/data/${STUDY}.intnpx.samples.tsv
OUT=${BASE}/match/scratch/${STUDY}.intnpx.pc${PC}.ukb.freeze12.cis.${ADJ}

## run match-prs-pheno (use --mahalanobis to adj for correlation)
${QPGEN} match-prs-pheno \
	--pheno ${NPX} \
	--prs $PRS \
	--sample-tsv ${IDMAP} \
	--out ${OUT} \
	--threads ${SLURM_CPUS_PER_TASK} \
	--lambda 1 \
	--mahalanobis

