### Stringtie: DEG Analysis ###
rm(list=ls())
options(stringsAsFactors = F)
library(tximport)
library(readr)
library(biomaRt)
library(DESeq2)

#读取Stringtie结果文件
dir <- "ballgown"
list.files(dir)
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
files <- file.path(dir, samples$run, "t_data.ctab")
names(files) <- samples$run
tmp <- read_tsv(files[1])
tx2gene <- tmp[, c("t_name", "gene_name")]
txi <- tximport(files, type = "stringtie", tx2gene = tx2gene)


#构建DESeqDataSet
colnames(txi$counts)
condition=factor(c(rep("A315T",2), rep("pCDH", 2)))
sampleTable <- data.frame(condition, row.names = colnames(txi$counts))
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ condition)

#DEseq2





