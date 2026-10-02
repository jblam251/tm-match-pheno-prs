#!/usr/bin/Rscript

## 09-09-26 : a slightly different implementation of RNAseq pre-processing which
##            better mimics the GTEx workflow. This is primarily to ensure there's
##            compatability when comparing the two (GTEx vs TOPMed for construction
##            of PGS). Although it'll likely preform better than the current pre-
##            processing stratagy which uses raw gene counts
##
## USAGE
## $ pheno2pl_rnaseq [TPM_CSV] [PREFIX]
##
## DEPENDENCIES
library(edgeR)



## initialize arg vector
args=commandArgs(trailingOnly = T)
print(paste0("script called"))

## import phenotype files
dat=as.data.frame(data.table::fread(args[1], sep="\t", header=T))
print(paste0("expression data imported"))

# if phenotypes already appear in PLINK format, write copy + exit
if(identical(colnames(dat)[1:2],c("FID","IID"))) {
  print(paste0("gene expression already appears to be in PLINK format, exiting.."))
  write.table(dat,paste0(args[2], ".wide.tsv"), sep = "\t", col.names=T, row.names=F, quote=F)
  quit(save="no", status=0)
}

## find/rm any duplicate traits
if(any(duplicated(dat[,2]))) {dat=dat[-which(duplicated(dat[,2])),]}

## store traits and samples
kp.smp=colnames(dat)[-c(1,2)]
kp.trait=dat[,2]
dat[,c(1,2)]=NULL

## TMM normalization
print(paste0("performing TMM normalization"))
dat.norm=DGEList(counts = dat)
dat.norm=calcNormFactors(dat.norm, method = "TMM")
dat.norm=cpm(dat.norm, normalized.lib.sizes = TRUE)
rownames(dat.norm)=kp.trait

## filters on low abundence genes
keep=rowSums(dat.norm>0.1) >= (0.20 * ncol(dat))
kp.trait=names(which(keep==TRUE))
dat.norm=dat.norm[keep, , drop = FALSE]

## perform RINT
print(paste0("performing RINT"))
rint <- function(x) {
  n <- length(x)
  qnorm((rank(x, ties.method="average") - 0.5) / n)
}
dat.norm=t(apply(dat.norm, 1, rint))
rownames(dat.norm)=kp.trait
colnames(dat.norm)=kp.smp


## transpose data set
dat.norm=as.data.frame(t(dat.norm))

## now get to PLINK format
dat.norm=cbind(kp.smp, kp.smp, dat.norm)
colnames(dat.norm)=c("FID","IID",kp.trait)
rownames(dat.norm)=1:nrow(dat.norm)

## any duplicate samples in phenotypes?
dup.samples=c()
if(any(duplicated(dat.norm[,1]))) {
  dup.samples=dat.norm[,1][which(duplicated(dat.norm[,1]))]
  dat.norm=dat.norm[-which(dat.norm[,1]%in%dup.samples),]
}

# print some stats
writeLines(paste0("number of samples in phenotype file : ", nrow(dat.norm)))
writeLines(paste0("number of assays in phenotype file : ", (ncol(dat.norm)-2)))
writeLines(paste0("duplicate samples identified and removed : ", length(dup.samples)))
writeLines(paste0(dup.samples, collapse=","))

# write new data
write.table(dat.norm, paste0(args[2], ".wide.tsv"), sep = "\t", col.names=T, row.names=F, quote=F)
gc(); Sys.sleep(5)
print(paste0("tpm expression data written"))
