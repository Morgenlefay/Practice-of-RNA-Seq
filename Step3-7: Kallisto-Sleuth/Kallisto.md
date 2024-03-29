### Environment
```bash
conda create -n Kallisto
conda activate Kallisto
conda install -c bioconda kallisto
```
### Index
```bash
kallisto index -i Reference/Kallisto/hg38_ensembl_cds.idx Reference/Kallisto/Homo_sapiens.GRCh38.cds.all.fa
```
### Quant
```bash
kallisto quant -i Reference/Kallisto/hg38_ensembl_cds.idx -o KA/pCDH_1 -t 8 -b 100 RNA_seq/pCDH_1_1.fq.gz RNA_seq/pCDH_1_2.fq.gz
kallisto quant -i Reference/Kallisto/hg38_ensembl_cds.idx -o KA/pCDH_2 -t 8 -b 100 RNA_seq/pCDH_2_1.fq.gz RNA_seq/pCDH_2_2.fq.gz
kallisto quant -i Reference/Kallisto/hg38_ensembl_cds.idx -o KA/A315T_1 -t 8 -b 100 RNA_seq/A315T_1_1.fq.gz RNA_seq/A315T_1_2.fq.gz
kallisto quant -i Reference/Kallisto/hg38_ensembl_cds.idx -o KA/A315T_2 -t 8 -b 100 RNA_seq/A315T_2_1.fq.gz RNA_seq/A315T_2_2.fq.gz
kallisto quant -i Reference/Kallisto/hg38_ensembl_cds.idx -o KA/Scrb_1 -t 8 -b 100 RNA_seq/Scrb_1_1.fq.gz RNA_seq/Scrb_1_2.fq.gz
kallisto quant -i Reference/Kallisto/hg38_ensembl_cds.idx -o KA/Scrb_2 -t 8 -b 100 RNA_seq/Scrb_2_1.fq.gz RNA_seq/Scrb_2_2.fq.gz
kallisto quant -i Reference/Kallisto/hg38_ensembl_cds.idx -o KA/sh2_1 -t 8 -b 100 RNA_seq/sh2_1_1.fq.gz RNA_seq/sh2_1_2.fq.gz
kallisto quant -i Reference/Kallisto/hg38_ensembl_cds.idx -o KA/sh2_2 -t 8 -b 100 RNA_seq/sh2_2_1.fq.gz RNA_seq/sh2_2_2.fq.gz
kallisto quant -i Reference/Kallisto/hg38_ensembl_cds.idx -o KA/sh8_1 -t 8 -b 100 RNA_seq/sh8_1_1.fq.gz RNA_seq/sh8_1_2.fq.gz
kallisto quant -i Reference/Kallisto/hg38_ensembl_cds.idx -o KA/sh8_2 -t 8 -b 100 RNA_seq/sh8_2_1.fq.gz RNA_seq/sh8_2_2.fq.gz
```
