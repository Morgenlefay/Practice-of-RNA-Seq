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
      






