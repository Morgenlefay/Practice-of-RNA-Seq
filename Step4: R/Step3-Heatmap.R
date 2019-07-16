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
dat <- read.csv("A315T_DEG.csv",header = TRUE, row.names = 1)
deg <- dat[ ,-1:-6]
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



rm(list=ls())
a<-read.csv(file="A315T.csv",header = T)
b<-read.csv(file="A315T_DEG.csv",header = T)
c<-merge(a,b,by=c("ID"))
write.csv(c, file="Venn.csv",row.names = F)
rm(list=ls())
options(stringsAsFactors = F)
dat <- read.csv("Venn.csv",header = TRUE,row.names = 1)
deg <- dat[ ,-1:-6]
x=dat$log2FoldChange
names(x)=rownames(deg)
cg=c(names(head(sort(x),200)),
     names(tail(sort(x),200)))
library(pheatmap)
pheatmap(deg[cg,],show_colnames =T,show_rownames = T)
n=t(scale(t(deg[cg,])))
n[n>2]=2
n[n< -2]= -2
pheatmap(n,show_colnames =T,show_rownames = F)
