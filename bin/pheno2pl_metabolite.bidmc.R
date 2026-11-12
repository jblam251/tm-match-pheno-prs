#!/usr/bin/Rscript

## USAGE
## Rscript pheno2pl_metabolite.R [PEAK AREAS] [METABOLITE ANNOTATION] [TRAIT FILE] [OUT]

## initialize arg vector
args=commandArgs(trailingOnly = T)
#args[1]="/net/topmed11/incoming/topmed/metabolomics/broad/CAMP.2024.0125.modified.data/CAMP_metabolomics_peakArea_20240125.txt.gz"
#args[2]="/net/topmed11/incoming/topmed/metabolomics/broad/CAMP.2024.0125.modified.data/CAMP_metabolomics_chemicalAnnotation_20240125.txt"
#args[3]="/net/topmed11/working/jblamer/id-checking/tm-match-pheno-prs/resources/metab.qtl.2026.phase2.traits.txt"
#args[4]="/net/topmed11/working/jblamer/qc.xqtl/match/metab.bidmc/data/camp.broad.peakarea.pc0.noimpute.norint"
print(paste0("script called"))

# read in metabolite traits (ie from summary statistics)
traits=as.character(read.table(args[3], header=F, sep="\n")[,1])
##*    a 'trait' file is needed to harmonize 
##*    the metabolite labels in the peak area
##*    file to the metabolite labels featured
##*    in the QTL summary results. this trait
##*    file should contain the metabolites as
##*    they appear in summary stats. this is
##*    less intuitive than say RNAseq or prot
##*    data where gene names or UniProt IDs
##*    are more portable


# read in metabolon peak areas
dat=read.csv(args[1], header=T, sep="\t")

# determine the structure of the peak areas
if(length(grep("^TOM", colnames(dat)))/ncol(dat)>0.10) {
  print("Metabolites detected in rows, samples in columns -- transposing peak areas...")
  dat=as.data.frame(t(dat))
  colnames(dat)=dat[1,]
  dat=dat[-1,]
  dat=cbind(SampleID=rownames(dat), dat)
}

# ensure metabolites in peak areas follow the format [ionmode].[compoundID]
colnames(dat)=gsub("-",".", colnames(dat))

# if phenotypes already appear in PLINK format, write copy + exit
if(identical(colnames(dat)[1:2],c("FID","IID"))) {
  print(paste0("phenotypes already appear to be in PLINK format, exiting.."))
  write.table(dat, paste0(args[4], ".wide.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
  quit(save="no", status=0)
}

# read in Broad metabolite annotation file
chemanno=read.csv(args[2], header=T, sep="\t")

# add key to match the notation in the peak area data
chemanno$Method.Compound_ID=paste0(gsub("-",".",chemanno$Method), ".", chemanno$Compound_ID)

# add final column which contains either HMDB (if available) or chem name -- this will be primary key for matching to "traits" file
chemanno$chem.match=chemanno$HMDB_ID

# if there's no HMDB, let's pull the chemical name over into primary match key column
chemanno$chem.match[-grep("HMDB", chemanno$chem.match)]=chemanno$Metabolite[-grep("HMDB", chemanno$chem.match)]

# remove metabolites that don't have either an HMDB or chemical name
chemanno=chemanno[-which(chemanno$chem.match==""),]

# for metabolites that don't have HMDBs, remove punctuation and add prefix "HMDB." to chemical name
chemanno$chem.match[-grep("HMDB", chemanno$chem.match)]=paste0("HMDB.", tolower(gsub("[^[:alnum:]]", "", chemanno$chem.match[-grep("HMDB", chemanno$chem.match)])))

# remove duplicated metabolites (usually due to same metabolite in multiple ionization modes)
chemanno=chemanno[-which(duplicated(chemanno$chem.match)),]

# 549 metabolites detected match
length(intersect(traits, chemanno$chem.match))
# 531 metabolites overlap between "traits" and HMDB ID
length(intersect(traits, chemanno$HMDB_ID))
# 18 metabolites overlap between "traits" and chemical name
length(intersect(traits, paste0("HMDB.",tolower(gsub("[^[:alnum:]]", "", chemanno$Metabolite)))))

# re-format the peak areas -- start by removing the TOPMed IDs in column 1
tmp.header=data.frame(FID=dat[,1], IID=dat[,1])
dat[,1]=NULL

# remove all metabolites in peak areas that don't have labels that match "traits" file
dat=dat[,which(colnames(dat)%in%chemanno$Method.Compound_ID)]

# now re-map the metabolite labels in peak areas into HMDB IDs
colnames(dat)=chemanno$chem.match[match(colnames(dat), chemanno$Method.Compound_ID)]

# add back the sample ids
dat=cbind(tmp.header, dat)

# any duplicate samples in phenotypes?
dup.samples=c()
if(any(duplicated(dat[,1]))) {
  dup.samples=dat[,1][which(duplicated(dat[,1]))]
  dat=dat[-which(dat[,1]%in%dup.samples),]
}

# print some stats
writeLines(paste0("traits detected in summary stats trait file : ", length(traits)))
writeLines(paste0("traits matched in peak areas : ", length(intersect(traits, chemanno$chem.match))))
writeLines(paste0("resulting number of missing data points in peak area file : ", length(which(is.na(dat))), " (", round(100*((length(which(is.na(dat))))/(nrow(dat)*(ncol(dat)-2))),2), "%)"))
writeLines(paste0("duplicate samples identified and removed : ", length(dup.samples)))
writeLines(paste0(dup.samples, collapse=","))


# write reforamtted pa data in regenie format
write.table(dat, paste0(args[4], ".wide.tsv"), sep="\t", col.names=T, row.names=F, quote=F)
writeLines(paste0("reforamtted peak area data written to : ", args[4]))


