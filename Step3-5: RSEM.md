### Preparing Reference Sequences
```bash
rsem-prepare-reference --gtf Reference/GTF/gencode.v19.annotation.gtf \
		                   --bowtie2 \
                         Reference/Fasta/hg19.fa Ref/human_gencode
```
### Calculating Expression Values
```bash
rsem-calculate-expression --bowtie2 \
                          --paired-end \
			  RNA_seq/pCDH_1_1_val_1.fq.gz \
			  RNA_seq/pCDH_1_2_val_2.fq.gz \
			  -p 8 Ref/human_gencode pCDH_1
                        
rsem-calculate-expression --bowtie2 \
                          --paired-end \
			  RNA_seq/pCDH_2_1_val_1.fq.gz \
			  RNA_seq/pCDH_2_2_val_2.fq.gz \
			  -p 8 Ref/human_gencode pCDH_2

rsem-calculate-expression --bowtie2 \
                          --paired-end \
			  RNA_seq/A315T_1_1_val_1.fq.gz \
			  RNA_seq/A315T_1_2_val_2.fq.gz \
			  -p 8 Ref/human_gencode A315T_1

rsem-calculate-expression --bowtie2 \
                          --paired-end \
			  RNA_seq/A315T_2_1_val_1.fq.gz \
			  RNA_seq/A315T_2_2_val_2.fq.gz \
			  -p 8 Ref/human_gencode A315T_2
```
### Differential Expression Analysis
#### generate-data-matrix
```bash
rsem-generate-data-matrix \
pCDH_1.genes.results pCDH_2.genes.results \
A315T_1.genes.results A315T_2.genes.results \
> GeneMat.txt
```
#### run-ebseq
```bash
rsem-run-ebseq \
GeneMat.txt 2,2 GeneMat.results
```
#### control_fdr
```bash
rsem-control-fdr \
GeneMat.results 0.05 GeneMat.de.txt
```




