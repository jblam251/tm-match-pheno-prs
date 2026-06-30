#!/usr/bin/Rscript

# define empty arg vector
args=commandArgs(trailingOnly = TRUE)

## FOR TESTING BEG ##
#args[1]="/net/topmed3/working/exchange.area.mirror/freeze.12.sources/TOPMed_Combined_Omics_SampleAttributes_DS_20260212.txt"
#args[2]="/net/topmed11/working/proteomics_subgroup/jblamer/npx.data/20250422/samples/FHS_OlinkHT_11182024.samples.tsv"
#args[3]="/net/topmed11/working/jblamer/id-checking/tm-match-pheno-prs/tmp.output.files/temp.output.testing"
## FOR TESTING END ##

# import files
md6=read.table(args[1], sep="\t", header=T)
smp=read.table(args[2], sep="\t", header=F)

# generate ID map file
idmap=md6[match(smp$V1, md6$SAMPLE_ID),c(7,11)]

# write
write.table(args[3], sep="\t", col.names=F, row.names=F, quote=F)

