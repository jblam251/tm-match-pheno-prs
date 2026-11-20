#!/usr/bin/Rscript
##
## Peak areas from BCM
##
## USAGE
## Rscript pheno2pl_metabolite.R [PEAK AREAS] [METABOLITE ANNOTATION] [TRAIT FILE] [OUT]

## initialize arg vector
args=commandArgs(trailingOnly = T)
print(paste0("script called"))


#args[1]="/net/topmed11/incoming/topmed/metabolomics/baylor/HCHSSOL.2025.0515/BatchNormData_DS.txt"
#args[2]="/net/topmed11/incoming/topmed/metabolomics/baylor/HCHSSOL.2025.0515/ChemicalAnnotation_DS.txt"
#args[3]="/net/topmed11/working/jblamer/qc.xqtl/match/metab.phase2/analysis/TOPMed_SOP_Phase2_Catalog_Final_01152025.Suggested.match.final.tsv"
#args[4]="/net/topmed11/working/jblamer/qc.xqtl/match/scratch/metab.catalog/data/hchssol.metabolon.peakarea.pc0.noimpute.norint"


# read in metabolite traits (ie from pgs file)
traits=read.csv(args[3], header=T, sep="\t")


# read in metabolon peak areas
dat=read.csv(args[1], header=T, sep="\t")
dat=as.data.frame(dat)


# if phenotypes already appear in PLINK format, write copy + exit
if(identical(colnames(dat)[1:2],c("FID","IID"))) {
  print(paste0("phenotypes already appear to be in PLINK format, exiting.."))
  write.table(dat, paste0(args[4], ".wide.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
  quit(save="no", status=0)
}

# read in metabolon metabolite annotation file
chemanno=read.csv(args[2], header=T, sep="\t")

# Use Bing's catalog to match metabolite names to QTL nomenclature
chemanno$Match.Catalog.Name=traits$Suggested.match.final[match(chemanno$CHEMICAL_NAME, traits$Metabolon_Name)]

# rename assay column names in peak area file to match prs
colnames(dat)=chemanno$Match.Catalog.Name[match(colnames(dat), paste0("X",chemanno$CHEM_ID))]
colnames(dat)[1]="SAMPLE_ID"

# remove missing or duplicate columns
if(any(is.na(colnames(dat)))) {dat=dat[,-which(is.na(colnames(dat)))]}
if(length(colnames(dat)[grep("[.]1", colnames(dat))])!=0) {dat=dat[,-grep("[.]1", colnames(dat))]}

# create FID and IID columns
dat$FID=dat$SAMPLE_ID
dat$IID=dat$SAMPLE_ID
dat$SAMPLE_ID=NULL
dat=dat[,c((ncol(dat)-1), ncol(dat), 1:(ncol(dat)-2))]

# any duplicate samples in phenotypes?
dup.samples=c()
if(any(duplicated(dat[,1]))) {
  dup.samples=dat[,1][which(duplicated(dat[,1]))]
  dat=dat[-which(dat[,1]%in%dup.samples),]
}

# print some stats
writeLines(paste0("traits detected in PRS trait file : ", nrow(traits)))
writeLines(paste0("traits matched in peak areas : ", length(intersect(traits$Suggested.match.final, chemanno$Match.Catalog.Name))))
writeLines(paste0("resulting number of missing data points in peak area file : ", length(which(is.na(dat))), " (", round(100*((length(which(is.na(dat))))/(nrow(dat)*(ncol(dat)-2))),2), "%)"))
writeLines(paste0("duplicate samples identified and removed : ", length(dup.samples)))
writeLines(paste0(dup.samples, collapse=","))


# write reforamtted pa data in regenie format
write.table(dat, paste0(args[4], ".wide.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
writeLines(paste0("reforamtted peak area data written to : ", args[4]))
