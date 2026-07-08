#!/usr/bin/Rscript

## USAGE
## $ Rscript gen.idmap.R [combined omics table] [phenotypes] [output]

# define empty arg vector
args=commandArgs(trailingOnly = TRUE)
print("script called...")

# import files
md6=read.table(args[1], sep="\t", header=T)
smp=read.table(args[2], sep="\t", header=F)[-1,1]
gc(); Sys.sleep(5)
print("files imported...")

# generate ID map file
idmap=md6[match(smp, md6$SAMPLE_ID),c(7,11)]

# remove phenotypes which are lacking NWDs
if(any(is.na(idmap[,1]))) {idmap=idmap[-which(is.na(idmap[,1])),]}
if(any(is.na(idmap[,2]))) {idmap=idmap[-which(is.na(idmap[,2])),]}
if(any(idmap[,1]=="")) {idmap=idmap[-which(idmap[,1]==""),]}

# write
write.table(idmap, args[3], sep="\t", col.names=F, row.names=F, quote=F)
print(paste0("idmap written to ", args[3]))
