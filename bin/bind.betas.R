#!/usr/bin/Rscript

## USAGE
## Rscript bind.betas.R [dir with beta noob csvs] [prefix]

## set args
args=commandArgs(trailingOnly = T)

## import beta files and bind into one table
lf=list.files(args[1], pattern="*.beta.noob.csv", full.names = T)
for (i in 1:length(lf)) {
  if(i==1) {
    beta.out=read.csv(lf[i], sep="\t", header=T)
    if(ncol(beta.out)==1) {beta.out=read.csv(lf[i], sep=",", header=T)}
    colnames(beta.out)[2]=regmatches(colnames(beta.out)[2], regexpr("TOE[[:alnum:]]{6}", colnames(beta.out)[2]))
  }
  if(i!=1) {
    tmpf.in=read.csv(lf[i], sep="\t", header=T)
    if(ncol(tmpf.in)==1) {tmpf.in=read.csv(lf[i], sep=",", header=T)}
    tmpf.in=tmpf.in[match(beta.out[,1], tmpf.in[,1]),]
    beta.out=cbind(beta.out, tmpf.in[,2])
    colnames(beta.out)[i+1]=regmatches(colnames(tmpf.in)[2], regexpr("TOE[[:alnum:]]{6}", colnames(tmpf.in)[2]))
    rm(tmpf.in)
  }
}


## get to plink format
traits=beta.out[,1]
beta.out[,1]=NULL
smps=colnames(beta.out)
beta.out=t(as.matrix(beta.out))
beta.out=data.frame(FID=smps, IID=smps, beta.out)
rownames(beta.out)=1:nrow(beta.out)
colnames(beta.out)=c("FID","IID",traits)

# any duplicate samples in phenotypes?
dup.samples=c()
if(any(duplicated(beta.out[,1]))) {
	dup.samples=beta.out[,1][which(duplicated(beta.out[,1]))]
	beta.out=beta.out[-which(beta.out[,1]%in%dup.samples),]
}

# print some stats
writeLines(paste0("number of samples in phenotype file : ", nrow(beta.out)))
writeLines(paste0("number of assays in phenotype file : ", (ncol(beta.out)-2)))
writeLines(paste0("duplicate samples identified and removed : ", length(dup.samples)))
writeLines(paste0(dup.samples, collapse=","))

## now write table
write.table(beta.out, paste0(args[2], ".beta.noob.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
