rm(list=ls())
options(stringsAsFactors = F)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
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
