### Environment
```bash
conda create -n SUPPA
conda activate SUPPA
conda install -c bioconda suppa
conda install -c bioconda salmon
```
### Index
```bash
salmon index -t Reference/Salmon/Homo_sapiens.GRCh38.cdna.all.fa -i Reference/Salmon/Ensembl_hg38_salmon_index
```

### Quantification
```bash
salmon quant -i Reference/Salmon/Gencode_hg19_salmon_index -l IU -1 RNA_seq/pCDH_1_1.fq.gz -2 RNA_seq/pCDH_1_2.fq.gz --validateMappings -o pCDH_1
```

