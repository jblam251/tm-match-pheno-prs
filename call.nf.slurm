#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --output=slurm.out/job.log.out

BASE=/net/topmed11/working/jblamer/id-checking/tm-match-pheno-prs
NF=${BASE}/main.p3.nf
CONFIG=${BASE}/etc/nf.config.pqtl

nextflow $NF -c $CONFIG -resume 
