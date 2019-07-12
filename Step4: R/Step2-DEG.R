###--------------------------
rm(list=ls())
options(stringsAsFactors = F)
load(file = 'A315T_OE.Rdata')
library(DESeq2)

###DEseq标准化dds
dds <- DESeqDataSetFromMatrix(countData=CountData, colData=colData, design=~condition)
dds <- DESeq(dds)
res <- results(dds)
mcols(res, use.names = TRUE)
summary(res)
plotMA(res, ylim = c(-6,6))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

###提取差异分析结果
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
diff_gene_deseq2 <- row.names(diff_gene_deseq2)
resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file= "A315T_DEG.csv",row.names = F)
subset(res,padj < 0.01) -> diff
subset(diff,log2FoldChange < -1) -> down
subset(diff,log2FoldChange > 1) -> up
as.data.frame(down) -> down_gene
as.data.frame(up) -> up_gene
write.csv(up_gene, file="A315T_Up.csv",row.names = T)
write.csv(down_gene, file="A315T_Down.csv",row.names = T)
