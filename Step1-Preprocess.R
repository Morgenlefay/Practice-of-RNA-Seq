###--------------------------
rm(list=ls())
options(stringsAsFactors = F)

###构建表达矩阵
Mydata <- read.table("A315T_OE_Matrix.txt", header=TRUE, sep = '\t')
ENSEMBL <- gsub("\\.\\d*", "", Mydata$Geneid)
row.names(Mydata) <- ENSEMBL
CountData <- Mydata[ ,-1]
CountData <- CountData[ ,-5]
condition <- factor(c("pCDH","pCDH","A315T","A315T"))
colData <- data.frame(row.names=colnames(CountData), condition)

save(colData,CountData,file = 'input.Rdata')
