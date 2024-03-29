### Environment
```bash
conda create -n Sailfish
conda activate Sailfish
conda install -c bioconda sailfish
```
### Index
```bash
sailfish index -p 8 -t Reference/Sailfish/Homo_sapiens.GRCh38.cds.all.fa -o Reference/Sailfish/hg38_ensembl_cds.idx
```
### Quant
```bash
gunzip RNA_seq/pCDH_1_1.fq.gz
gunzip RNA_seq/pCDH_1_2.fq.gz
gunzip RNA_seq/pCDH_2_1.fq.gz
gunzip RNA_seq/pCDH_2_2.fq.gz
gunzip RNA_seq/A315T_1_1.fq.gz
gunzip RNA_seq/A315T_1_2.fq.gz
gunzip RNA_seq/A315T_2_1.fq.gz
gunzip RNA_seq/A315T_2_2.fq.gz
sailfish quant -i Reference/Sailfish/hg38_ensembl_cds.idx -l IU -1 RNA_seq/pCDH_1_1.fq -2 RNA_seq/pCDH_1_2.fq -o Sailfish/pCDH_1
sailfish quant -i Reference/Sailfish/hg38_ensembl_cds.idx -l IU -1 RNA_seq/pCDH_2_1.fq -2 RNA_seq/pCDH_2_2.fq -o Sailfish/pCDH_2
sailfish quant -i Reference/Sailfish/hg38_ensembl_cds.idx -l IU -1 RNA_seq/A315T_1_1.fq -2 RNA_seq/A315T_1_2.fq -o Sailfish/A315T_1
sailfish quant -i Reference/Sailfish/hg38_ensembl_cds.idx -l IU -1 RNA_seq/A315T_2_1.fq -2 RNA_seq/A315T_2_2.fq -o Sailfish/A315T_2
```

