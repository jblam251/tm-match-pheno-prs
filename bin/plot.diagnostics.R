#!/usr/bin/Rscript

library(ggplot2)
library(gridExtra)

## USAGE
## Rscript plot.diagnostics.R [results file]

## import results
args=commandArgs(trailingOnly = TRUE)
#args[1]="/net/topmed11/working/jblamer/id-checking/tm-match-pheno-prs/results/rnaseq.sagegalaii.results.match.all.tsv.gz"
res=read.csv(args[1], sep="\t", header=T)

## define color key
colkey=data.frame(v1=c("SELF_BEST","SELF_LENIENT","SINGLE_NEW_BEST","MULTI_NEW_BEST","UNCLEAR","NO_ASSIGNED_GT"), v2=c("#619CFF", "#F564E3", "#00BA38", "#F8766D", "#B79F00", "grey80"))

## generate plots
fout=gsub("results.match.all.tsv.gz","diagnostic.plot.png",strsplit(args[1], "/")[[1]][length(strsplit(args[1], "/")[[1]])])
#png(file=paste0("/net/topmed11/working/jblamer/id-checking/tm-match-pheno-prs/analysis/",fout), width=8, height=8, units="in", res=300)
png(file=fout, width=8, height=8, units="in", res=300)

res$MatchStatus[which(is.na(res$ID.self))]="NO_ASSIGNED_GT"
res$MatchStatus=factor(res$MatchStatus, levels=c("SELF_BEST","SELF_LENIENT","SINGLE_NEW_BEST","MULTI_NEW_BEST","UNCLEAR","NO_ASSIGNED_GT"))
status.tab=data.frame(as.data.frame(table(res$MatchStatus)), Group=c("PASS","PASS","WARN","WARN","WARN","NO.TEST"))
status.tab$Group=factor(status.tab$Group, levels=c("PASS","WARN","NO.TEST"))
if(any(status.tab$Freq==0)) {status.tab=status.tab[-which(status.tab$Freq==0),]}
cols=colkey$v2[match(status.tab$Var1, colkey$v1)]
gglist=list()
gglist[[4]]=ggplot(status.tab, aes(x=Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 0.75, color = "white") +
  scale_fill_manual(values=cols)+
  theme_minimal()+
  theme(legend.title = element_blank(),legend.text = element_text(size = 12),legend.key.size = unit(1, "cm"),legend.spacing.y = unit(0.4, "cm"))
gglist[[4]]=cowplot::get_legend(gglist[[4]])
gglist[[1]]=ggplot(status.tab[which(status.tab$Var1!="NO_ASSIGNED_GT"),], aes(x=Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 0.75, color = "white") +
  scale_fill_manual(values=cols)+
  xlab("")+ylab("Samples")+
  geom_text(aes(label = Freq),vjust = -0.5, size = 3)+
  theme_minimal()+
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  facet_wrap(~Group, scales="free", space = "free_x", strip.position = "bottom")
gglist[[2]]=ggplot(res, aes(x=Z.1st, y=Z.self, color=MatchStatus))+
  geom_point(size=3)+
  scale_color_manual(values=cols)+
  geom_abline(intercept=0, slope=1, color="grey80", lty=2)+
  ylab("Z-Score Self")+
  xlab("Z-Score Best Match")+
  theme_bw()+
  theme(legend.position ="none")
gglist[[3]]=ggplot(res[which(res$MatchStatus=="NO_ASSIGNED_GT"),], aes(x=Z.1st, y=Z.2nd))+
  geom_point(size=3, color="grey80")+
  geom_abline(intercept=0, slope=1, color="grey80", lty=2)+
  ylab("Z-Score 2nd Best Match")+
  xlab("Z-Score Best Match")+
  theme_bw()+
  theme(legend.position ="none")
do.call(grid.arrange, c(gglist, ncol = 2))  
dev.off()


