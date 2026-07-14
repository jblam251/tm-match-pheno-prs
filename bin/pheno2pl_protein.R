#!/usr/bin/Rscript

## USAGE
## $ csv2plink.$ [NPX CSV] [COL MAP FILE] [PREFIX]

## initialize arg vector
args=commandArgs(trailingOnly = T)
#library(data.table)
print(paste0("script called"))

## import npx and column-mapping files
dat=as.data.frame(data.table::fread(args[1], sep=",", header=T))
#dat=read.csv(args[1], sep=",", header=T)
map=read.csv(args[2], sep=",", header=T)
print(paste0("npx data imported"))

## remove sample controls
idx.smp.ctrl=which(colnames(dat)=="SampleType")
if(length(idx.smp.ctrl)!=0) {
  dat=dat[-grep("CONTROL", dat[[idx.smp.ctrl]]),]
}

## remove assay controls
idx.assay.ctrl=which(colnames(dat)=="AssayType")
if(length(idx.assay.ctrl)!=0) {
  dat=dat[-grep("ctrl", dat[[idx.assay.ctrl]]),]
}

# get unique assays and samples
assays=unique(dat[,map$ColumnNumber[which(map$Field=="GeneName")]])
smps=unique(dat[,map$ColumnNumber[which(map$Field=="TOP_ID")]])

## subset only sample, assay, and npx
dat=dat[,map$ColumnNumber]
gc(); Sys.sleep(5)

## initialize matrix for wide format mapping
mat=matrix(NA, nrow=length(smps), ncol=length(assays), dimnames=list(smps, assays))

# now fill matrix
mat[cbind(match(dat[[1]], smps),match(dat[[2]], assays))]=dat[[3]]
print(paste0("npx data reformatted to wide"))

# now get to PLINK format
mat=data.frame(FID=smps, IID=smps, mat)

# any duplicate samples in phenotypes?
dup.samples=c()
if(any(duplicated(mat[,1]))) {
	dup.samples=mat[,1][which(duplicated(mat[,1]))]
	mat=mat[-which(mat[,1]%in%dup.samples),]
}

# print some stats
writeLines(paste0("number of samples in phenotype file : ", nrow(mat)))
writeLines(paste0("number of assays in phenotype file : ", (ncol(mat)-2)))
writeLines(paste0("duplicate samples identified and removed : ", length(dup.samples)))
writeLines(paste0(dup.samples, collapse=","))

# write new data
write.table(mat, paste0(args[3], ".npx.wide.tsv"), sep = "\t", col.names=T, row.names=F, quote=F)
gc(); Sys.sleep(5)
print(paste0("npx data written"))
