#!/usr/bin/Rscript

## USAGE
## Rscript pheno2pl_metabolite.R [PEAK AREAS] [METABOLITE ANNOTATION] [TRAIT FILE] [OUT]

## initialize arg vector
args=commandArgs(trailingOnly = T)
print(paste0("script called"))

# read in metabolite traits (ie from pgs file)
traits=as.character(read.table(args[3], header=F, sep="\n")[,1])

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

#fix the metabolites with duplicated HMDBs (Bing resolved this by adding a suffex of ".2")
chemanno$HMDB[which(duplicated(chemanno$HMDB))]=paste0(chemanno$HMDB[which(duplicated(chemanno$HMDB))], ".2")
chemanno$HMDB[which(chemanno$HMDB==".2")]=""
#fix the metabolites that have multiple HMDBs seperated by a comma
chemanno$HMDB=gsub(",", ".", chemanno$HMDB)
#fix the metabolites that are assigned to multiple HMDBs
chemanno.split.hmdb=data.frame(
  chemanno[rep(seq_len(nrow(chemanno)), lengths(strsplit(chemanno$HMDB, ","))), ],
  HMDB.fix = unlist(strsplit(chemanno$HMDB, ","))
)
chemanno.split.hmdb$HMDB=chemanno.split.hmdb$HMDB.fix
chemanno.split.hmdb$HMDB.fix=NULL
chemanno=rbind(chemanno.split.hmdb, chemanno[which(chemanno$HMDB==""),])
rm(chemanno.split.hmdb)


# column HMDB should now match traits in PRS, now run the overlap
chemanno$chem.match=chemanno$HMDB
chemanno$chem.match[which(chemanno$chem.match%in%setdiff(chemanno$chem.match, traits))]=NA
chemanno$chem.match[which(is.na(chemanno$chem.match))]=paste0("HMDB.", tolower(gsub("[^[:alnum:]]", "", chemanno$CHEMICAL_NAME[which(is.na(chemanno$chem.match))])))
chemanno$chem.match[which(chemanno$chem.match%in%setdiff(chemanno$chem.match, traits))]=NA

# rename assay column names in peak area file to match prs
colnames(dat)=chemanno$chem.match[match(colnames(dat), paste0("X",chemanno$CHEM_ID))]
colnames(dat)[1]="SAMPLE_ID"
if(any(is.na(colnames(dat)))) {dat=dat[,-which(is.na(colnames(dat)))]}

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
writeLines(paste0("traits detected in PRS trait file : ", length(traits)))
writeLines(paste0("traits matched in peak areas : ", length(intersect(traits, chemanno$chem.match))))
writeLines(paste0("resulting number of missing data points in peak area file : ", length(which(is.na(dat))), " (", round(100*((length(which(is.na(dat))))/(nrow(dat)*(ncol(dat)-2))),2), "%)"))
writeLines(paste0("duplicate samples identified and removed : ", length(dup.samples)))
writeLines(paste0(dup.samples, collapse=","))


# write reforamtted pa data in regenie format
write.table(dat, paste0(args[4], ".wide.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
writeLines(paste0("reforamtted peak area data written to : ", args[4]))
