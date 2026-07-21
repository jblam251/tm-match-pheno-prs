#!/usr/bin/Rscript

## USAGE
## $ pheno2pl_rnaseq [PHENO CSV] [PREFIX]

## initialize arg vector
args=commandArgs(trailingOnly = T)
print(paste0("script called"))

## import phenotype files
dat=as.data.frame(data.table::fread(args[1], sep="\t", header=T))
print(paste0("expression data imported"))

## find/rm any duplicate traits
if(any(duplicated(dat[,2]))) {dat=dat[-which(duplicated(dat[,2])),]}

## store traits
kp.trait=dat[,2]

## store samples
kp.smp=colnames(dat)[-c(1,2)]

## transpose data set
dat[,c(1,2)]=NULL
dat=as.data.frame(t(dat))

## now get to PLINK format
dat=cbind(kp.smp, kp.smp, dat)
colnames(dat)=c("FID","IID",kp.trait)
rownames(dat)=1:nrow(dat)

## any duplicate samples in phenotypes?
dup.samples=c()
if(any(duplicated(dat[,1]))) {
  dup.samples=dat[,1][which(duplicated(dat[,1]))]
  dat=dat[-which(dat[,1]%in%dup.samples),]
}

# print some stats
writeLines(paste0("number of samples in phenotype file : ", nrow(dat)))
writeLines(paste0("number of assays in phenotype file : ", (ncol(dat)-2)))
writeLines(paste0("duplicate samples identified and removed : ", length(dup.samples)))
writeLines(paste0(dup.samples, collapse=","))

# write new data
write.table(dat, paste0(args[2], ".reads.wide.tsv"), sep = "\t", col.names=T, row.names=F, quote=F)
gc(); Sys.sleep(5)
print(paste0("gene expression data written"))
