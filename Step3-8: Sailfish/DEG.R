### Sailfish: DEG Analysis ###
rm(list=ls())
options(stringsAsFactors = F)
library(tximport)
library(readr)
library(biomaRt)
library(DESeq2)

#设置sailfish生成的路径
dir <- "Sailfish"
list.files(dir)
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
files <- file.path(dir, samples$run, "quant.sf")
names(files) <- samples$run

#获取转录本ID信息
tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))

#读取sailfish结果文件
txi <- tximport(files, type="sailfish", tx2gene=tx2gene)

#构建DESeqDataSet
colnames(txi$counts)
condition=factor(c(rep("A315T",2), rep("pCDH", 2)))
sampleTable <- data.frame(condition, row.names = colnames(txi$counts))
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ condition)

#DEseq2





