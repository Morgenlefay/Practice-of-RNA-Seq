### Build index
```bash
subread-buildindex -o hg19_subread.index Reference/Fasta/hg19.fa
```
### Alignment
```bash
subread-align -T 8 -i Reference/index/hg19_subread/hg19_subread.index -r RNA_seq/pCDH_1_1_val_1.fq.gz -R RNA_seq/pCDH_1_2_val_2.fq.gz -t 0 -o RNA_seq/pCDH_1.bam
subread-align -T 8 -i Reference/index/hg19_subread/hg19_subread.index -r RNA_seq/pCDH_2_1_val_1.fq.gz -R RNA_seq/pCDH_2_2_val_2.fq.gz -t 0 -o RNA_seq/pCDH_2.bam
subread-align -T 8 -i Reference/index/hg19_subread/hg19_subread.index -r RNA_seq/A315T_1_1_val_1.fq.gz -R RNA_seq/A315T_1_2_val_2.fq.gz -t 0 -o RNA_seq/A315T_1.bam
subread-align -T 8 -i Reference/index/hg19_subread/hg19_subread.index -r RNA_seq/A315T_2_1_val_1.fq.gz -R RNA_seq/A315T_2_2_val_2.fq.gz -t 0 -o RNA_seq/A315T_2.bam
subread-align -T 8 -i Reference/index/hg19_subread/hg19_subread.index -r RNA_seq/Scrb_1_1_val_1.fq.gz -R RNA_seq/Scrb_1_2_val_2.fq.gz -t 0 -o RNA_seq/Scrb_1.bam
subread-align -T 8 -i Reference/index/hg19_subread/hg19_subread.index -r RNA_seq/Scrb_2_1_val_1.fq.gz -R RNA_seq/Scrb_2_2_val_2.fq.gz -t 0 -o RNA_seq/Scrb_2.bam
subread-align -T 8 -i Reference/index/hg19_subread/hg19_subread.index -r RNA_seq/sh2_1_1_val_1.fq.gz -R RNA_seq/sh2_1_2_val_2.fq.gz -t 0 -o RNA_seq/sh2_1.bam
subread-align -T 8 -i Reference/index/hg19_subread/hg19_subread.index -r RNA_seq/sh2_2_1_val_1.fq.gz -R RNA_seq/sh2_2_2_val_2.fq.gz -t 0 -o RNA_seq/sh2_2.bam
subread-align -T 8 -i Reference/index/hg19_subread/hg19_subread.index -r RNA_seq/sh8_1_1_val_1.fq.gz -R RNA_seq/sh8_1_2_val_2.fq.gz -t 0 -o RNA_seq/sh8_1.bam
subread-align -T 8 -i Reference/index/hg19_subread/hg19_subread.index -r RNA_seq/sh8_2_1_val_1.fq.gz -R RNA_seq/sh8_2_2_val_2.fq.gz -t 0 -o RNA_seq/sh8_2.bam
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
