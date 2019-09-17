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
### HTseq
```bash
htseq-count -s no -r name -f bam RNA_seq/pCDH_1.bam Reference/GTF/gencode.v19.annotation.gtf > pCDH_1.count	
htseq-count -s no -r name -f bam RNA_seq/pCDH_2.bam Reference/GTF/gencode.v19.annotation.gtf > pCDH_2.count	
htseq-count -s no -r name -f bam RNA_seq/A315T_1.bam Reference/GTF/gencode.v19.annotation.gtf > A315T_1.count	
htseq-count -s no -r name -f bam RNA_seq/A315T_2.bam Reference/GTF/gencode.v19.annotation.gtf > A315T_2.count	
paste pCDH_1.count pCDH_2.count A315T_1.count A315T_2.count | awk '{printf $1"\t";for(i=2;i<=NF;i=i+2)printf $i"\t";print $i}' > A315T_OE_Matrix.count

htseq-count -s no -r name -f bam RNA_seq/Scrb_1.bam Reference/GTF/gencode.v19.annotation.gtf > Scrb_1.count	
htseq-count -s no -r name -f bam RNA_seq/Scrb_2.bam Reference/GTF/gencode.v19.annotation.gtf > Scrb_2.count	
htseq-count -s no -r name -f bam RNA_seq/sh2_1.bam Reference/GTF/gencode.v19.annotation.gtf > sh2_1.count	
htseq-count -s no -r name -f bam RNA_seq/sh2_2.bam Reference/GTF/gencode.v19.annotation.gtf > sh2_2.count	
paste Scrb_1.count Scrb_2.count sh2_1.count sh2_2.count | awk '{printf $1"\t";for(i=2;i<=NF;i=i+2)printf $i"\t";print $i}' > sh2_KD_Matrix.count

htseq-count -s no -r name -f bam RNA_seq/Scrb_1.bam Reference/GTF/gencode.v19.annotation.gtf > Scrb_1.count	
htseq-count -s no -r name -f bam RNA_seq/Scrb_2.bam Reference/GTF/gencode.v19.annotation.gtf > Scrb_2.count	
htseq-count -s no -r name -f bam RNA_seq/sh8_1.bam Reference/GTF/gencode.v19.annotation.gtf > sh8_1.count	
htseq-count -s no -r name -f bam RNA_seq/sh8_2.bam Reference/GTF/gencode.v19.annotation.gtf > sh8_2.count	
paste Scrb_1.count Scrb_2.count sh8_1.count sh8_2.count | awk '{printf $1"\t";for(i=2;i<=NF;i=i+2)printf $i"\t";print $i}' > sh8_KD_Matrix.count
```
