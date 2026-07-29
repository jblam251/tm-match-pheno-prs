#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=8:00:00
#SBATCH --output=slurm.out/job.log.out

BASE=/net/topmed11/working/jblamer/id-checking/tm-match-pheno-prs
TYPE=methylation  #one of proteomics, metabolomcs, rnaseq, or methylation
NF=${BASE}/main.p3.nf
CONFIG=${BASE}/etc/nf.config.runx
SETTINGS=${BASE}/etc/match.settings.run.csv

nextflow run $NF -c $CONFIG -resume \
	--type $TYPE \
	--settings $SETTINGS
