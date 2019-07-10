###--------------------------
rm(list=ls())
options(stringsAsFactors = F)
library(DESeq2)

###构建表达矩阵
Mydata <- read.table("A315T_OE_Matrix.txt", header=TRUE, sep = '\t')
ENSEMBL <- gsub("\\.\\d*", "", Mydata$Geneid)
row.names(Mydata) <- ENSEMBL
CountData <- Mydata[ ,-1]
CountData <- CountData[ ,-5]
condition <- factor(c("pCDH","pCDH","A315T","A315T"))
colData <- data.frame(row.names=colnames(CountData), condition)

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
write.csv(resdata,file= "DEG.csv",row.names = F)

subset(res,padj < 0.01) -> diff
subset(diff,log2FoldChange < -1) -> down
subset(diff,log2FoldChange > 1) -> up
as.data.frame(down) -> down_gene
as.data.frame(up) -> up_gene
write.csv(up_gene, file="Up_gene.csv",row.names = T)
write.csv(down_gene, file="Down_gene.csv",row.names = T)

###火山图-1
library(ggplot2)
volcano_data <-  read.csv("DEG.csv",header = TRUE)
loc_up <- intersect(which(volcano_data$padj<0.05),which(volcano_data$log2FoldChange>=1))
loc_down <- intersect(which(volcano_data$padj<0.05),which(volcano_data$log2FoldChange<=(-1)))
significant <- rep("Normal",times=nrow(volcano_data))
significant[loc_up] <- "Up"
significant[loc_down] <- "Down"
significant <- factor(significant,levels=c("Up","Down","Normal"))
p <- qplot(x=volcano_data$log2FoldChange,y=-log10(volcano_data$padj),xlab="log2(FC)",ylab="-log10(FDR)",colour=significant)
p <- p+ scale_color_manual(values=c("Up"="red","Normal"="black","Down"="green"))
xline=c(-log2(2),log2(2))
p <- p+geom_vline(xintercept=xline,lty=2,size=I(0.2),colour="grey11")
yline=-log(0.05,10)
p <- p+geom_hline(yintercept=yline,lty=2,size=I(0.2),colour="grey11")
p <- p+theme_bw()
p

###火山图-2
rm(list = ls())
library(ggplot2)
data <- read.csv("DEG.csv",header = TRUE)
data$color <- ifelse(data$padj<0.05 & abs(data$log2FoldChange)>= 1,ifelse(data$log2FoldChange > 1,'red','blue'),'black')
color <- c(red = "red",black = "black",blue = "blue")
p <- ggplot(data, aes(log2FoldChange, -log10(padj), col = color)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = color) +
  labs(x="log2 (fold change)",y="-log10 (q-value)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
p

###火山图-3
rm(list = ls())
library(EnhancedVolcano)
data <- read.csv("DEG.csv",header = TRUE, row.names = 1)
EnhancedVolcano(data,
                lab = rownames(data),
                x = "log2FoldChange",
                y = "padj",
                selectLab = c("ENSG00000182752","ENSG00000137558"),
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = 0.00000000000000001,
                FCcutoff = 2.0,
                xlim = c(-8,8),
                transcriptLabSize = 3.0,
                colAlpha = 1,
                legend=c("NS","Log2 FC","Adjusted p-value",
                         "Adjusted p-value & Log2 FC"),
                legendPosition = "bottom",
                legendLabSize = 10,
                legendIconSize = 3.0,
                DrawConnectors = FALSE,
                border = "full",
                borderWidth = 1.5,
                borderColour = "black",
                gridlines.major = FALSE,
                gridlines.minor = FALSE)

###功能富集分析
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
setwd("H:/RNA_seq")
sig.gene<-read.csv(file="DEG.csv")
head(sig.gene)
gene<-sig.gene[,1]
head(gene)
gene.df<-bitr(gene, fromType = "ENSEMBL", 
              toType = c("SYMBOL","ENTREZID"),
              OrgDb = org.Hs.eg.db)
head(gene.df)
ego_bp<-enrichGO(gene       = gene.df$ENSEMBL,
                 OrgDb      = org.Hs.eg.db,
                 keyType    = 'ENSEMBL',
                 ont        = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05)
barplot(ego_bp,showCategory = 25,title="The GO_BP enrichment analysis of all DEGs ")
dotplot(ego_bp,showCategory = 25,title="The GO_BP enrichment analysis of all DEGs ")
