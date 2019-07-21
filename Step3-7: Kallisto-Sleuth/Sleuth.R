--------#### DEG 
rm(list=ls())
options(stringsAsFactors = F)
library(rhdf5)
library(sleuth)
library(biomaRt)

#设置kallisto生成的路径
base_dir <- "KA"
#获取所有的simple_id
sample_id <- dir(file.path(base_dir))
#获取结果文件的路径
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
#读取实验设计表
s2c <- read.table("./design_matrix.txt", header = TRUE, sep='\t',stringsAsFactors=FALSE)
#与路径合并
s2c <- dplyr::mutate(s2c, path = kal_dirs)
#获取转录本ID信息
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

#读取Kallisto结果文件
so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g, extra_bootstrap_summary = TRUE)
#使用condition设计矩阵回归
so <- sleuth_fit(so)
#使用截距项回归
so <- sleuth_fit(so, ~1, 'reduced')
#使用LRT进行鉴定
so <- sleuth_lrt(so, 'reduced', 'full')
# LRT检验结果
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)

                   
                   
                   
--------#### GO-KEGG-DO Analysis                   
rm(list=ls())
options(stringsAsFactors = F)
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(org.Hs.eg.db)
library(DOSE)
sig.gene<-read.csv(file="DEG.csv")
gene<-sig.gene[,2]
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
