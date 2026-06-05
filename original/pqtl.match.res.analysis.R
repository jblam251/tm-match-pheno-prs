


## import  *assigned results files
res.list=list()
res.files=list.files("/net/topmed11/working/jblamer/qc.xqtl/match/results/match.freeze12c.pgs", pattern="*.pc50.*match.assigned.tsv.gz", full.names = T)
res.files=res.files[-grep("mesa1-3", res.files)]
# read in each file, pull data
for (i in 1:length(res.files)) {
  res.list[[i]]=read.table(res.files[i], sep="\t", header=T)
  res.list[[i]]$MatchStatus=factor(res.list[[i]]$MatchStatus, levels=c("BEST_MATCH","LENIENT_MATCH","NO_MATCH"))
}
res.df=do.call("rbind.data.frame", res.list)
res.df$file=rep(sapply(strsplit(res.files, "/"), '[', 10), unlist(lapply(res.list, nrow)))



## import *all results files
res.all.list=list()
res.all.files=list.files("/net/topmed11/working/jblamer/qc.xqtl/match/pqtl/results/match.freeze12c.pgs", pattern="*.pc50.*match.all.tsv.gz", full.names = T)
res.all.files=res.all.files[-grep("mesa1-3", res.all.files)]
# need idmap file (since the *all.tsv columns are misaligned and the ID.self needs to be mapped)
idmap=read.table("/net/topmed11/working/jblamer/qc.xqtl/match/pqtl/idmaps/pqtl.combined.idmap.cardia.unmatched.nwds.rm.mesa.fails.added.tsv", sep="\t", header=F)
# read in each file, pull data
for (i in 1:length(res.all.files)) {
  res.all.list[[i]]=read.csv(res.all.files[i], sep="\t", header=T)
  #fix the shifted column issue with these "all" files
  colnames(res.all.list[[i]])=colnames(res.all.list[[i]])[c(1,2,4:21,3)]
  res.all.list[[i]]=res.all.list[[i]][,c(1,2,21,3:20)]
  res.all.list[[i]]$ID.self=idmap$V1[match(res.all.list[[i]]$ID.Pheno, idmap$V2)]
}
# convert to df for ease
res.all.df=do.call("rbind.data.frame", res.all.list)
res.all.df$file=rep(sapply(strsplit(res.all.files, "/"), '[', 11), unlist(lapply(res.all.list, nrow)))
res.all.df$ISMATCH=res.all.df$ID.self==res.all.df$ID.1st
res.all.df$studyname=sapply(strsplit(res.all.df$file, "[.]"), "[", 1)
res.all.df$mahalanobis=sapply(strsplit(res.all.df$file, "[.]"), "[", 7)




## create df to store results
df.out=data.frame(study=sapply(strsplit(sapply(strsplit(res.files, "/"), "[", 10), "[.]"), '[', 1),
                  pc.adj=sapply(strsplit(sapply(strsplit(res.files, "/"), "[", 10), "[.]"), '[', 3),
                  mahalanobis=sapply(strsplit(sapply(strsplit(res.files, "/"), "[", 10), "[.]"), '[', 7),
                  BEST_MATCH=NA, LENIENT_MATCH=NA, NO_MATCH=NA,
                  BEST_MATCH.perc=NA, LENIENT_MATCH.perc=NA, NO_MATCH.perc=NA)
# get match counts for each file
for (i in 1:length(res.files)) {
  df.out[i,4:6]=unclass(table(res.list[[i]]$MatchStatus))
  df.out[i,7:9]=as.numeric(unclass(df.out[i,4:6]/sum(df.out[i,4:6])))
}
# sum the best and lenient counts
df.out$BEST.LENIENT_MATCH=df.out$BEST_MATCH+df.out$LENIENT_MATCH
df.out$BEST.LENIENT_MATCH.perc=df.out$BEST_MATCH.perc+df.out$LENIENT_MATCH.perc
df.out$KEY=paste0(df.out$pc.adj,"_",df.out$mahalanobis)
df.out=df.out[order(df.out$study),]
# add a row to tally totals across studies
df.total=df.out[FALSE,]
# df.total[1,]=c("TOTAL","pc0","default",as.numeric(colSums(df.out[which(df.out$KEY=="pc0_default"),c(4:6)])),
#                as.numeric(colSums(df.out[which(df.out$KEY=="pc0_default"),c(4:6)]) / sum(colSums(df.out[which(df.out$KEY=="pc0_default"),c(4:6)]))),
#                sum(df.out[which(df.out$KEY=="pc0_default"),c(10)]),
#                sum(df.out[which(df.out$KEY=="pc0_default"),c(10)]) / sum(colSums(df.out[which(df.out$KEY=="pc0_default"),c(4:6)])),
#                "pc0_default")
df.total[1,]=c("TOTAL","pc50","default",as.numeric(colSums(df.out[which(df.out$KEY=="pc50_default"),c(4:6)])),
               as.numeric(colSums(df.out[which(df.out$KEY=="pc50_default"),c(4:6)]) / sum(colSums(df.out[which(df.out$KEY=="pc50_default"),c(4:6)]))),
               sum(df.out[which(df.out$KEY=="pc50_default"),c(10)]),
               sum(df.out[which(df.out$KEY=="pc50_default"),c(10)]) / sum(colSums(df.out[which(df.out$KEY=="pc50_default"),c(4:6)])),
               "pc50_default")
# df.total[3,]=c("TOTAL","pc0","mahalanobis",as.numeric(colSums(df.out[which(df.out$KEY=="pc0_mahalanobis"),c(4:6)])),
#                as.numeric(colSums(df.out[which(df.out$KEY=="pc0_mahalanobis"),c(4:6)]) / sum(colSums(df.out[which(df.out$KEY=="pc0_mahalanobis"),c(4:6)]))),
#                sum(df.out[which(df.out$KEY=="pc0_mahalanobis"),c(10)]),
#                sum(df.out[which(df.out$KEY=="pc0_mahalanobis"),c(10)]) / sum(colSums(df.out[which(df.out$KEY=="pc0_mahalanobis"),c(4:6)])),
#                "pc0_mahalanobis")
df.total[2,]=c("TOTAL","pc50","mahalanobis",as.numeric(colSums(df.out[which(df.out$KEY=="pc50_mahalanobis"),c(4:6)])),
               as.numeric(colSums(df.out[which(df.out$KEY=="pc50_mahalanobis"),c(4:6)]) / sum(colSums(df.out[which(df.out$KEY=="pc50_mahalanobis"),c(4:6)]))),
               sum(df.out[which(df.out$KEY=="pc50_mahalanobis"),c(10)]),
               sum(df.out[which(df.out$KEY=="pc50_mahalanobis"),c(10)]) / sum(colSums(df.out[which(df.out$KEY=="pc50_mahalanobis"),c(4:6)])),
               "pc50_mahalanobis")
df.out=rbind(df.out, df.total)



## plot % match between default and mahalanobis
# #v1 - list og ggplots
# gg.defvmah=data.frame(STUDY=unique(df.out$study), 
#                       BEST.DEFAULT=as.numeric(df.out$BEST_MATCH.perc[which(df.out$mahalanobis=="default")]),
#                       BEST.MAHALANOBIS=as.numeric(df.out$BEST_MATCH.perc[which(df.out$mahalanobis=="mahalanobis")]),
#                       LENIENT.DEFAULT=as.numeric(df.out$BEST.LENIENT_MATCH.perc[which(df.out$mahalanobis=="default")]),
#                       LENIENT.MAHALANOBIS=as.numeric(df.out$BEST.LENIENT_MATCH.perc[which(df.out$mahalanobis=="mahalanobis")]))
# gg.defvmah.list=list()
# gg.defvmah.list[[1]]=ggplot(gg.defvmah, aes(y=BEST.DEFAULT, x=BEST.MAHALANOBIS, label=STUDY))+
#   geom_text_repel(size=5)+
#   geom_point()+
#   ylab("Best Match % (Default)")+
#   xlab("Best Match % (Mahalanobis)")+
#   theme_bw()
# gg.defvmah.list[[2]]=ggplot(gg.defvmah, aes(y=LENIENT.DEFAULT, x=LENIENT.MAHALANOBIS, label=STUDY))+
#   geom_text_repel(size=5)+
#   geom_point()+
#   ylab("Lenient Match % (Default)")+
#   xlab("Lenient Match % (Mahalanobis)")+
#   theme_bw()
# gridExtra::grid.arrange(grobs=gg.defvmah.list, nrow=1)
#v2 - facet wrap match type
gg.defvmah2=data.frame(STUDY=rep(unique(df.out$study),2), 
                       TYPE=c(rep("BEST.MATCH", 11), rep("LENIENT.MATCH", 11)),
                       DEFAULT=c(as.numeric(df.out$BEST_MATCH.perc[which(df.out$mahalanobis=="default")]),
                                 as.numeric(df.out$BEST.LENIENT_MATCH.perc[which(df.out$mahalanobis=="default")])),
                       MAHALANOBIS=c(as.numeric(df.out$BEST_MATCH.perc[which(df.out$mahalanobis=="mahalanobis")]),
                                     as.numeric(df.out$BEST.LENIENT_MATCH.perc[which(df.out$mahalanobis=="mahalanobis")])))
ggplot(gg.defvmah2, aes(y=DEFAULT, x=MAHALANOBIS, label=STUDY))+
  geom_text_repel(size=5, max.overlaps = 10)+
  geom_point()+
  ylab("Match % (Default)")+
  xlab("Match % (Mahalanobis)")+
  theme_bw()+
  geom_abline(lty=2, col="grey75")+
  facet_wrap(~TYPE)







## plot % mismatch across conditions (one per study)
# #list of ggplots
# gg.df=data.frame(STUDY=rep(df.out$study,2), COND=rep(df.out$KEY), MATCH.TYPE=rep(c("strict","lenient"), each=nrow(df.out)), PERC=c(df.out$BEST_MATCH.perc, df.out$BEST.LENIENT_MATCH.perc))
# gg.list=list()
# for (i in 1:length(unique(gg.df$STUDY))) {
#   gg2list=gg.df[which(gg.df$STUDY==unique(gg.df$STUDY)[i]),]
#   #gg2list$COND=factor(gg2list$COND, levels=c("pc0_default","pc0_mahalanobis","pc50_default","pc50_mahalanobis"))
#   gg2list$COND=factor(gg2list$COND, levels=c("pc50_default","pc50_mahalanobis"))
#   gg.list[[i]]=ggplot(gg2list, aes(x=MATCH.TYPE, y=(100*(1-as.numeric(PERC))), fill=COND))+
#     geom_bar(position="dodge", stat="identity")+
#     ylab("% Mismatch")+xlab("")+
#     theme_bw()+
#     theme(legend.position="none")+
#     ggtitle(unique(gg.df$STUDY)[i])
# }
# gridExtra::grid.arrange(grobs=gg.list, nrow=3)
#or facet wrap
gg.df=data.frame(STUDY=rep(df.out$study,2), COND=rep(df.out$KEY), MATCH.TYPE=rep(c("strict","lenient"), each=nrow(df.out)), PERC=c(df.out$BEST_MATCH.perc, df.out$BEST.LENIENT_MATCH.perc))
ggplot(gg.df, aes(x=MATCH.TYPE, y=(100*(1-as.numeric(PERC))), fill=COND))+
  geom_bar(position="dodge", stat="identity")+
  ylab("% Mismatch")+xlab("")+
  theme_bw()+
  facet_wrap(~STUDY)



## plot best vs 2nd z-score
# best vs 2nd z-score for all NWDs -- facet wrap studies -- plot default and mahalanobis seperately
ggplot(res.all.df[which(res.all.df$mahalanobis=="default"),], aes(x=Z.1st, y=Z.2nd, color=ISMATCH))+
  geom_point(size=1)+
  xlab("Z-Score best match")+
  ylab("Z-Score 2nd best match")+
  geom_abline(lty=2)+
  ggtitle("DEFAULT")+
  theme_bw()+
  facet_wrap(~studyname)
ggplot(res.all.df[which(res.all.df$mahalanobis=="mahalanobis"),], aes(x=Z.1st, y=Z.2nd, color=ISMATCH))+
  geom_point(size=1)+
  xlab("Z-Score best match")+
  ylab("Z-Score 2nd best match")+
  geom_abline(lty=2)+
  ggtitle("MAHALANOBIS")+
  theme_bw()+
  facet_wrap(~studyname)
# best vs 2nd z-score for unassigned NWDs -- facet wrap studies -- plot default and mahalanobis seperately
ggplot(res.all.df[which(res.all.df$mahalanobis=="default"&is.na(res.all.df$ID.self)),], aes(x=Z.1st, y=Z.2nd, color=ISMATCH))+
  geom_point(size=1)+
  xlab("Z-Score best match")+
  ylab("Z-Score 2nd best match")+
  geom_abline(lty=2)+
  ggtitle("DEFAULT")+
  theme_bw()+
  facet_wrap(~studyname)
ggplot(res.all.df[which(res.all.df$mahalanobis=="mahalanobis"&is.na(res.all.df$ID.self)),], aes(x=Z.1st, y=Z.2nd, color=ISMATCH))+
  geom_point(size=1)+
  xlab("Z-Score best match")+
  ylab("Z-Score 2nd best match")+
  geom_abline(lty=2)+
  ggtitle("MAHALANOBIS")+
  theme_bw()+
  facet_wrap(~studyname)

# table of unmapped NWDs with clear (2x next best z-score) best matches
table(res.all.df$studyname[which((res.all.df$Z.1st>(2*res.all.df$Z.2nd))&is.na(res.all.df$ID.self))])
table(res.all.df$studyname[is.na(res.all.df$ID.self)])
# table of unmapped NWDs due to lack of metadata with clear (2x next best z-score) best matches
md6=read.table("/net/topmed3/working/exchange.area.mirror/freeze.12.sources/TOPMed_Combined_Omics_SampleAttributes_DS_20260212.txt", sep="\t", header=T)
res.all.df.missing.meta=res.all.df[which(res.all.df$ID.Pheno%in%setdiff(res.all.df$ID.Pheno, md6$SAMPLE_ID)),]
table(res.all.df.missing.meta$studyname[which((res.all.df.missing.meta$Z.1st>(2*res.all.df.missing.meta$Z.2nd))&is.na(res.all.df.missing.meta$ID.self))])
length(which((res.all.df.missing.meta$Z.1st>(2*res.all.df.missing.meta$Z.2nd))&is.na(res.all.df.missing.meta$ID.self)))
table(res.all.df.missing.meta$studyname[is.na(res.all.df.missing.meta$ID.self)])
length(which(is.na(res.all.df.missing.meta$ID.self)))
#*for protein extended cohorts: 
##   camp  cardia fhs-bcm    ltrc mesa1-4 
##    35      86      25     457      36 
 
#35+86+25+457+36 +53+10+51

## generate summary table for visualization
df.out.4vis2=data.frame(STUDY=sort(unique(df.out$study)), TotSamples=NA, WGS.Matched=NA)
df.out.4vis2$TotSamples[1:10]=unique(unlist(lapply(res.all.list, nrow)))[c(1,2,3,4,5,7,6,8,9,10)]
df.out.4vis2$TotSamples[11]=sum(unique(unlist(lapply(res.all.list, nrow)))[c(1,2,3,4,5,7,6,8,9,10)])
df.out.4vis2$WGS.Matched[1:10]=unique(unlist(lapply(res.list, nrow)))[c(1,2,3,4,5,7,6,8,9,10)]
df.out.4vis2$WGS.Matched[11]=sum(df.out.4vis2$WGS.Matched[1:10])
df.out.4vis2$WGS.Matched.perc=paste0(round(100*(df.out.4vis2$WGS.Matched/df.out.4vis2$TotSamples),2), "%")
#df.out.4vis2$PC0.DEFAULT.BEST=df.out$BEST_MATCH[which(df.out$pc.adj=="pc0"&df.out$mahalanobis=="default")]
#df.out.4vis2$PC0.DEFAULT.BEST.perc=paste0(round(100*(as.numeric(df.out$BEST_MATCH.perc[which(df.out$pc.adj=="pc0"&df.out$mahalanobis=="default")])),2), "%")
#df.out.4vis2$PC0.DEFAULT.LENIENT=df.out$BEST.LENIENT_MATCH[which(df.out$pc.adj=="pc0"&df.out$mahalanobis=="default")]
#df.out.4vis2$PC0.DEFAULT.LENIENT.perc=paste0(round(100*(as.numeric(df.out$BEST.LENIENT_MATCH.perc[which(df.out$pc.adj=="pc0"&df.out$mahalanobis=="default")])),2), "%")
#df.out.4vis2$PC0.MAHALANOBIS.BEST=df.out$BEST_MATCH[which(df.out$pc.adj=="pc0"&df.out$mahalanobis=="mahalanobis")]
#df.out.4vis2$PC0.MAHALANOBIS.BEST.perc=paste0(round(100*(as.numeric(df.out$BEST_MATCH.perc[which(df.out$pc.adj=="pc0"&df.out$mahalanobis=="mahalanobis")])),2), "%")
#df.out.4vis2$PC0.MAHALANOBIS.LENIENT=df.out$BEST.LENIENT_MATCH[which(df.out$pc.adj=="pc0"&df.out$mahalanobis=="mahalanobis")]
#df.out.4vis2$PC0.MAHALANOBIS.LENIENT.perc=paste0(round(100*(as.numeric(df.out$BEST.LENIENT_MATCH.perc[which(df.out$pc.adj=="pc0"&df.out$mahalanobis=="mahalanobis")])),2), "%")
df.out.4vis2$PC50.DEFAULT.BEST=df.out$BEST_MATCH[which(df.out$pc.adj=="pc50"&df.out$mahalanobis=="default")]
df.out.4vis2$PC50.DEFAULT.BEST.perc=paste0(round(100*(as.numeric(df.out$BEST_MATCH.perc[which(df.out$pc.adj=="pc50"&df.out$mahalanobis=="default")])),2), "%")
df.out.4vis2$PC50.MAHALANOBIS.BEST=df.out$BEST_MATCH[which(df.out$pc.adj=="pc50"&df.out$mahalanobis=="mahalanobis")]
df.out.4vis2$PC50.MAHALANOBIS.BEST.perc=paste0(round(100*(as.numeric(df.out$BEST_MATCH.perc[which(df.out$pc.adj=="pc50"&df.out$mahalanobis=="mahalanobis")])),2), "%")
df.out.4vis2$PC50.DEFAULT.LENIENT=df.out$BEST.LENIENT_MATCH[which(df.out$pc.adj=="pc50"&df.out$mahalanobis=="default")]
df.out.4vis2$PC50.DEFAULT.LENIENT.perc=paste0(round(100*(as.numeric(df.out$BEST.LENIENT_MATCH.perc[which(df.out$pc.adj=="pc50"&df.out$mahalanobis=="default")])),2), "%")
df.out.4vis2$PC50.MAHALANOBIS.LENIENT=df.out$BEST.LENIENT_MATCH[which(df.out$pc.adj=="pc50"&df.out$mahalanobis=="mahalanobis")]
df.out.4vis2$PC50.MAHALANOBIS.LENIENT.perc=paste0(round(100*(as.numeric(df.out$BEST.LENIENT_MATCH.perc[which(df.out$pc.adj=="pc50"&df.out$mahalanobis=="mahalanobis")])),2), "%")
openxlsx::write.xlsx(t(df.out.4vis2), "/net/topmed11/working/jblamer/qc.xqtl/match/analysis/pqtl.match.res.fz12.tab.xlsx", sep="\t", colNames=F, rowNames=T, quote=F)
#res.all.df[which(res.all.df$ID.Pheno=="TOP285217"),]


## what are the fail samples? Any systemic plate errors?
# import various *samples.tsv files
smp.list=list()
smp.files=list.files("/net/topmed11/working/jblamer/qc.xqtl/match/data", pattern="*samples.tsv", full.names = T)
smp.files=smp.files[-grep("mesa1-3", smp.files)]
# read in each file, pull data
for (i in 1:length(smp.files)) {
  smp.list[[i]]=read.table(smp.files[i], sep="\t", header=F)
  smp.list[[i]]$STUDY=strsplit(strsplit(smp.files[i], "/")[[1]][9], "[.]")[[1]][1]
}
#rm the extra col in MESA
smp.list[[10]]$V4=NULL
#bind 
lapply(smp.list, head)
smp.df=do.call("rbind.data.frame", smp.list)
#make a unique plateID for the smp files
smp.df$UniqPlateID=paste0(smp.df$STUDY,"-x-",smp.df$V3)
#get list of fail samples
qc.fails=c()
for (i in 1:length(res.list)) {
  if(i%in%grep("*.pc50.*mahalanobis*", res.files)) {
    #qc.fails=c(qc.fails, res.list[[i]]$ID.Pheno[which(res.list[[i]]$MatchStatus!="BEST_MATCH")]) #stringent fails
    qc.fails=c(qc.fails, res.list[[i]]$ID.Pheno[which(res.list[[i]]$MatchStatus=="NO_MATCH")])  #lenient fails
  }
}
plate.df=data.frame(study_plateID=unique(smp.df$UniqPlateID))
plate.count.tab=as.data.frame(unclass(table(smp.df$UniqPlateID)))
plate.df$no.samples=plate.count.tab[,1][match(plate.df$study_plateID, rownames(plate.count.tab))]
plate.fail.tab=as.data.frame(unclass(table(smp.df$UniqPlateID[match(qc.fails, smp.df$V1)])))
plate.df$no.fails=plate.fail.tab[,1][match(plate.df$study_plateID, rownames(plate.fail.tab))]
plate.df$no.fails[which(is.na(plate.df$no.fails))]=0
plate.df$perc.plate.fail=plate.df$no.fails/plate.df$no.samples
plate.df[order(plate.df$perc.plate.fail, decreasing=T)[1:20],]




## what about the weights?
# import various *weights files
wght.list=list()
wght.files=list.files("/net/topmed11/working/jblamer/qc.xqtl/match/results/match.freeze12c.pgs", pattern="*.pc50.*weights.tsv.gz", full.names = T)
wght.files=wght.files[-grep("mesa1-3", wght.files)]
wght.files=wght.files[-grep("default", wght.files)]
# read in each file, pull data
for (i in 1:length(wght.files)) {
  wght.list[[i]]=read.table(wght.files[i], sep="\t", header=T)
  wght.list[[i]]$STUDY=strsplit(strsplit(wght.files[i], "/")[[1]][10], "[.]")[[1]][1]
}
wght.df=do.call("rbind.data.frame", wght.list)
table(wght.df$STUDY)
table(wght.df$Trait, wght.df$STUDY)
cor.wght.df=data.frame(assay=sort(unique(wght.df$Trait)))
for (i in 1:length(wght.list)) {
  cor.wght.df[,(i+1)]=wght.list[[i]]$Weight[match(cor.wght.df$assay, wght.list[[i]]$Trait)]
  colnames(cor.wght.df)[(i+1)]=wght.list[[i]]$STUDY[1]
}
rownames(cor.wght.df)=cor.wght.df$assay
cor.wght.df$assay=NULL
pheatmap::pheatmap(cor(cor.wght.df, use="complete.obs", method = "spearman"), cluster_rows=T, cluster_cols=T, angle_col = 90)



## what are the z-scores with these lenient matches?
res.df.lenient=res.df[which(res.df$MatchStatus=="LENIENT_MATCH"),]
res.df.lenient=res.df.lenient[grep("mahalanobis", res.df.lenient$file),]
res.df.lenient$set=sapply(strsplit(res.df.lenient$file, "[.]"), '[', 1)
ggplot(res.df.lenient, aes(x=Rank.self, y=Z.self, color=MatchStatus))+
  geom_point()+
  theme_bw()+
  facet_wrap(~set)




##

length(setdiff(res.all.list[[1]]$ID.Pheno, res.list[[1]]$ID.Pheno))
f12$UNIQ_SUBJECT_ID2=paste0(f12$PHS,"_",f12$SUBJECT_ID)
f10$UNIQ_SUBJECT_ID2=paste0(f10$PHS,"_",f10$SUBJECT_ID)
md6$UNIQ_SUBJECT_ID[which(md6$SAMPLE_ID%in%setdiff(res.all.list[[1]]$ID.Pheno, res.list[[1]]$ID.Pheno))]
nrow(f12[which(f12$UNIQ_SUBJECT_ID2%in%md6$UNIQ_SUBJECT_ID[which(md6$SAMPLE_ID%in%setdiff(res.all.list[[1]]$ID.Pheno, res.list[[1]]$ID.Pheno))]),])
nrow(f10[which(f10$UNIQ_SUBJECT_ID2%in%md6$UNIQ_SUBJECT_ID[which(md6$SAMPLE_ID%in%setdiff(res.all.list[[1]]$ID.Pheno, res.list[[1]]$ID.Pheno))]),])
