rm(list=ls())
options(stringsAsFactors = F)
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(org.Hs.eg.db)
library(DOSE)
sig.gene<-read.csv(file="DEG.csv")
gene<-sig.gene[,1]
gene.df<-bitr(gene, fromType = "ENSEMBL", 
              toType = c("SYMBOL","ENTREZID"),
              OrgDb = org.Hs.eg.db)
head(gene.df)

### GO Analysis
ego_bp<-enrichGO(gene       = gene.df$ENSEMBL,
                 OrgDb      = org.Hs.eg.db,
                 keyType    = 'ENSEMBL',
                 ont        = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05)
barplot(ego_bp,showCategory = 25,title="The GO_BP enrichment analysis of all DEGs ")
dotplot(ego_bp,showCategory = 25,title="The GO_BP enrichment analysis of all DEGs ")
plotGOgraph(ego_bp)
ego_mf<-enrichGO(gene       = gene.df$ENSEMBL,
                 OrgDb      = org.Hs.eg.db,
                 keyType    = 'ENSEMBL',
                 ont        = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05)
barplot(ego_mf,showCategory = 25,title="The GO_MF enrichment analysis of all DEGs")
dotplot(ego_mf,showCategory = 25,title="The GO_MF enrichment analysis of all DEGs")
plotGOgraph(ego_mf)

### KEGG Analysis
kegg<-gene.df$ENTREZID
kk <- enrichKEGG(gene = kegg,
                 organism = 'hsa',
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05)
dotplot(kk)
barplot(kk)
browseKEGG(kk, 'hsa04974')

### DO Analysis
do<-gene.df$SYMBOL
do <- enrichDO(gene = kegg,
              ont = "DO",
              pvalueCutoff = 0.01,
              pAdjustMethod = "BH",
              minGSSize = 1,
              maxGSSize = 500,
              qvalueCutoff = 0.05,
              readable = TRUE)
dotplot(do)
barplot(do)
















