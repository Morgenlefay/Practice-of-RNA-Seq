rm(list=ls())
options(stringsAsFactors = F)
load(file = 'A315T_OE.Rdata')
library(pheatmap)
pheatmap(scale(cor(log2(CountData+1))))
colnames(CountData)
pheatmap::pheatmap(cor(CountData))
colData
tmp=data.frame(g=colData)
rownames(tmp)=colnames(CountData)
pheatmap::pheatmap(cor(CountData),annotation_col = tmp)

rm(list=ls())
options(stringsAsFactors = F)
dat<-read.csv("A315T_DEG.csv")
ENSEMBL <- dat$Row.names
row.names(dat) <- ENSEMBL
deg <- dat[ ,-1:-7]
x=dat$log2FoldChange
names(x)=rownames(deg)
cg=c(names(head(sort(x),100)),
     names(tail(sort(x),100)))
library(pheatmap)
pheatmap(deg[cg,],show_colnames =T,show_rownames = T)
n=t(scale(t(deg[cg,])))
n[n>2]=2
n[n< -2]= -2
pheatmap(n,show_colnames =T,show_rownames = T)
