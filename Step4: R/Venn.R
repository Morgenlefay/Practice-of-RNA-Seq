rm(list=ls())
install.packages("VennDiagram")
library(VennDiagram)
rt1<-read.csv(file="A315T.csv")
A315T <- as.vector(rt1[,1])
rt2<-read.csv(file="sh2.csv")
sh2 <- as.vector(rt2[,1])
rt3<-read.csv(file="sh8.csv")
sh8 <- as.vector(rt3[,1])
vn <- venn.diagram(list(A315T=A315T,sh2=sh2,sh8=sh8),filename='Venn.eps',
                   cat.col=c('red','green','blue'),
                   fill=c(colors()[148], colors()[589], colors()[116]))



A315T_sh2 <- intersect(A315T,sh2)
write.table(file="A315T_sh2.xls",A315T_sh2,sep="\t",quote=F,col.names=F,row.names=F)

