### Hisat2
```bash
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 RNA_seq/pCDH_1_1_val_1.fq.gz -2 RNA_seq/pCDH_1_2_val_2.fq.gz | samtools sort -@4 -O bam -o RNA_seq/pCDH_1.bam	
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 RNA_seq/pCDH_2_1_val_1.fq.gz -2 RNA_seq/pCDH_2_2_val_2.fq.gz | samtools sort -@4 -O bam -o RNA_seq/pCDH_2.bam	
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 RNA_seq/A315T_1_1_val_1.fq.gz -2 RNA_seq/A315T_1_2_val_2.fq.gz | samtools sort -@4 -O bam -o RNA_seq/A315T_1.bam	
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 RNA_seq/A315T_2_1_val_1.fq.gz -2 RNA_seq/A315T_2_2_val_2.fq.gz | samtools sort -@4 -O bam -o RNA_seq/A315T_2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 RNA_seq/Scrb_1_1_val_1.fq.gz -2 RNA_seq/Scrb_1_2_val_2.fq.gz | samtools sort -@4 -O bam -o RNA_seq/Scrb_1.bam	
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 RNA_seq/Scrb_2_1_val_1.fq.gz -2 RNA_seq/Scrb_2_2_val_2.fq.gz | samtools sort -@4 -O bam -o RNA_seq/Scrb_2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 RNA_seq/sh2_1_1_val_1.fq.gz -2 RNA_seq/sh2_1_2_val_2.fq.gz | samtools sort -@4 -O bam -o RNA_seq/sh2_1.bam	
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 RNA_seq/sh2_2_1_val_1.fq.gz -2 RNA_seq/sh2_2_2_val_2.fq.gz | samtools sort -@4 -O bam -o RNA_seq/sh2_2.bam
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 RNA_seq/sh8_1_1_val_1.fq.gz -2 RNA_seq/sh8_1_2_val_2.fq.gz | samtools sort -@4 -O bam -o RNA_seq/sh8_1.bam	
hisat2 -t -p 8 -x Reference/index/hg19/genome -1 RNA_seq/sh8_2_1_val_1.fq.gz -2 RNA_seq/sh8_2_2_val_2.fq.gz | samtools sort -@4 -O bam -o RNA_seq/sh8_2.bam
```
### FeatureCounts
```bash
featureCounts -T 8 -p -t exon -g gene_id -a Reference/GTF/gencode.v19.annotation.gtf -o pCDH_1.txt RNA_seq/pCDH_1.bam	
featureCounts -T 8 -p -t exon -g gene_id -a Reference/GTF/gencode.v19.annotation.gtf -o pCDH_2.txt RNA_seq/pCDH_2.bam	
featureCounts -T 8 -p -t exon -g gene_id -a Reference/GTF/gencode.v19.annotation.gtf -o A315T_1.txt RNA_seq/A315T_1.bam	
featureCounts -T 8 -p -t exon -g gene_id -a Reference/GTF/gencode.v19.annotation.gtf -o A315T_2.txt RNA_seq/A315T_2.bam
featureCounts -T 8 -p -t exon -g gene_id -a Reference/GTF/gencode.v19.annotation.gtf -o Scrb_1.txt RNA_seq/Scrb_1.bam	
featureCounts -T 8 -p -t exon -g gene_id -a Reference/GTF/gencode.v19.annotation.gtf -o Scrb_2.txt RNA_seq/Scrb_2.bam	
featureCounts -T 8 -p -t exon -g gene_id -a Reference/GTF/gencode.v19.annotation.gtf -o sh2_1.txt RNA_seq/sh2_1.bam	
featureCounts -T 8 -p -t exon -g gene_id -a Reference/GTF/gencode.v19.annotation.gtf -o sh2_2.txt RNA_seq/sh2_2.bam
featureCounts -T 8 -p -t exon -g gene_id -a Reference/GTF/gencode.v19.annotation.gtf -o sh8_1.txt RNA_seq/sh8_1.bam	
featureCounts -T 8 -p -t exon -g gene_id -a Reference/GTF/gencode.v19.annotation.gtf -o sh8_2.txt RNA_seq/sh8_2.bam

cut -f 1,7 pCDH_1.txt |grep -v '^#' >RNA_seq/Count/pCDH_1.txt
cut -f 1,7 pCDH_2.txt |grep -v '^#' >RNA_seq/Count/pCDH_2.txt
cut -f 1,7 A315T_1.txt |grep -v '^#' >RNA_seq/Count/A315T_1.txt
cut -f 1,7 A315T_2.txt |grep -v '^#' >RNA_seq/Count/A315T_2.txt
cut -f 1,7 Scrb_1.txt |grep -v '^#' >RNA_seq/Count/Scrb_1.txt
cut -f 1,7 Scrb_2.txt |grep -v '^#' >RNA_seq/Count/Scrb_2.txt
cut -f 1,7 sh2_1.txt |grep -v '^#' >RNA_seq/Count/sh2_1.txt
cut -f 1,7 sh2_2.txt |grep -v '^#' >RNA_seq/Count/sh2_2.txt
cut -f 1,7 sh8_1.txt |grep -v '^#' >RNA_seq/Count/sh8_1.txt
cut -f 1,7 sh8_2.txt |grep -v '^#' >RNA_seq/Count/sh8_2.txt
paste RNA_seq/Count/pCDH_1.txt RNA_seq/Count/pCDH_2.txt RNA_seq/Count/A315T_1.txt RNA_seq/Count/A315T_2.txt | awk '{printf $1"\t";for(i=2;i<=NF;i=i+2)printf $i"\t";print $i}' > A315T_OE_Matrix.txt
paste RNA_seq/Count/Scrb_1.txt RNA_seq/Count/Scrb_2.txt RNA_seq/Count/sh2_1.txt RNA_seq/Count/sh2_2.txt | awk '{printf $1"\t";for(i=2;i<=NF;i=i+2)printf $i"\t";print $i}' > sh2_KD_Matrix.txt
paste RNA_seq/Count/Scrb_1.txt RNA_seq/Count/Scrb_2.txt RNA_seq/Count/sh8_1.txt RNA_seq/Count/sh8_2.txt | awk '{printf $1"\t";for(i=2;i<=NF;i=i+2)printf $i"\t";print $i}' > sh8_KD_Matrix.txt
```
### Bam2Bed
```bash
bedtools bamtobed -i RNA_seq/pCDH_1.bam  > RNA_seq/pCDH_1.bed
bedtools bamtobed -i RNA_seq/pCDH_1.bam  > RNA_seq/pCDH_2.bed
bedtools bamtobed -i RNA_seq/A315T_1.bam  > RNA_seq/A315T_1.bed
bedtools bamtobed -i RNA_seq/A315T_2.bam  > RNA_seq/A315T_2.bed
bedtools bamtobed -i RNA_seq/Scrb_1.bam  > RNA_seq/Scrb_1.bed
bedtools bamtobed -i RNA_seq/Scrb_2.bam  > RNA_seq/Scrb_2.bed
bedtools bamtobed -i RNA_seq/sh2_1.bam  > RNA_seq/sh2_1.bed
bedtools bamtobed -i RNA_seq/sh2_2.bam  > RNA_seq/sh2_2.bed
bedtools bamtobed -i RNA_seq/sh8_1.bam  > RNA_seq/sh8_1.bed
bedtools bamtobed -i RNA_seq/sh8_2.bam  > RNA_seq/sh8_2.bed
```
### Bam2Bigwig
```bash
chmod 775 bam2bigwig.sh
./bam2bigwig.sh RNA_seq/pCDH_1.bam
./bam2bigwig.sh RNA_seq/pCDH_2.bam
./bam2bigwig.sh RNA_seq/A315T_1.bam
./bam2bigwig.sh RNA_seq/A315T_2.bam
./bam2bigwig.sh RNA_seq/Scrb_1.bam
./bam2bigwig.sh RNA_seq/Scrb_2.bam
./bam2bigwig.sh RNA_seq/sh2_1.bam
./bam2bigwig.sh RNA_seq/sh2_2.bam
./bam2bigwig.sh RNA_seq/sh8_1.bam
./bam2bigwig.sh RNA_seq/sh8_2.bam
```

