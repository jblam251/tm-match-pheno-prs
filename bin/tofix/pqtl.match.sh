#!/bin/bash
#SBATCH --mem=48G
#SBATCH --cpus-per-task=12
#SBATCH --time=01:00:00
#SBATCH --output=%j.out

#study prefixes: *cardia,*chs,*mesa1-4,*ltrc,*cra,*camp,*fhs,*fhs-bcm,*copdgene,*jhs-batch1

#when comparing PGS to full freeze12, set to 8 threads! 

STUDY=ltrc
PC=50
ADJ=default #one of default or mahalanobis
QPGEN=/net/fantasia/home/hmkang/code/working/qpgen/bin/qpgentools
BASE=/net/topmed11/working/jblamer/qc.xqtl
PRS=${BASE}/match/pgs/ukb_ppp_v2_sentinel_HT.topmed_freeze12c_minDP0_cis.prs.tsv.gz
NPX=${BASE}/data/${STUDY}.intnpx.pc${PC}.regenie.tsv.gz
IDMAP=${BASE}/match/idmaps/pqtl.combined.idmap.cardia.unmatched.nwds.rm.mesa.fails.added.tsv
SAMPLES=${BASE}/data/${STUDY}.intnpx.samples.tsv
OUT=${BASE}/match/scratch/${STUDY}.intnpx.pc${PC}.ukb.freeze12.cis.${ADJ}

# if we want to check all NWDs in IDMAP match PGS...
#cut -f1 $IDMAP | sort | uniq > pqtls.nwds.all
#$(wc -l pqtls.nwds.all)
#$(zcat $PRS | grep -f pqtl.nwds.all | wc -l)

# subset the idmap file
cut -f1 $SAMPLES > ${BASE}/match/tmpf.sample.list.${SLURM_JOB_ID}
grep -f ${BASE}/match/tmpf.sample.list.${SLURM_JOB_ID} $IDMAP > ${BASE}/match/tmpf.idmap.${SLURM_JOB_ID}

## run match-prs-pheno (use --mahalanobis to adj for correlation)
if [[ "$ADJ" == "default" ]]; then
	${QPGEN} match-prs-pheno --pheno ${NPX} \
		--prs $PRS \
		--sample-tsv ${BASE}/match/tmpf.idmap.${SLURM_JOB_ID} \
		--out ${OUT} \
		--threads ${SLURM_CPUS_PER_TASK}
fi
if [[ "$ADJ" == "mahalanobis" ]]; then
	${QPGEN} match-prs-pheno --pheno ${NPX} \
		--prs $PRS \
		--sample-tsv ${BASE}/match/tmpf.idmap.${SLURM_JOB_ID} \
		--out ${OUT} \
		--threads ${SLURM_CPUS_PER_TASK} \
		--mahalanobis
fi



## rm tmp files
rm ${BASE}/match/tmpf.sample.list.${SLURM_JOB_ID}
rm ${BASE}/match/tmpf.idmap.${SLURM_JOB_ID}
