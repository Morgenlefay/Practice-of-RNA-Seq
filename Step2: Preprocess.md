### FastQC and MultiQC of Raw Data
```bash
fastqc -f fastq RNA_seq/pCDH_1_1.fq.gz -o RNA_seq/FastQC_1	
fastqc -f fastq RNA_seq/pCDH_1_2.fq.gz -o RNA_seq/FastQC_1	
fastqc -f fastq RNA_seq/pCDH_2_1.fq.gz -o RNA_seq/FastQC_1	
fastqc -f fastq RNA_seq/pCDH_2_2.fq.gz -o RNA_seq/FastQC_1
fastqc -f fastq RNA_seq/A315T_1_1.fq.gz -o RNA_seq/FastQC_1	
fastqc -f fastq RNA_seq/A315T_1_2.fq.gz -o RNA_seq/FastQC_1	
fastqc -f fastq RNA_seq/A315T_2_1.fq.gz -o RNA_seq/FastQC_1	
fastqc -f fastq RNA_seq/A315T_2_2.fq.gz -o RNA_seq/FastQC_1
fastqc -f fastq RNA_seq/Scrb_1_1.fq.gz -o RNA_seq/FastQC_1	
fastqc -f fastq RNA_seq/Scrb_1_2.fq.gz -o RNA_seq/FastQC_1	
fastqc -f fastq RNA_seq/Scrb_2_1.fq.gz -o RNA_seq/FastQC_1	
fastqc -f fastq RNA_seq/Scrb_2_2.fq.gz -o RNA_seq/FastQC_1
fastqc -f fastq RNA_seq/sh2_1_1.fq.gz -o RNA_seq/FastQC_1	
fastqc -f fastq RNA_seq/sh2_1_2.fq.gz -o RNA_seq/FastQC_1	
fastqc -f fastq RNA_seq/sh2_2_1.fq.gz -o RNA_seq/FastQC_1	
fastqc -f fastq RNA_seq/sh2_2_2.fq.gz -o RNA_seq/FastQC_1
fastqc -f fastq RNA_seq/sh8_1_1.fq.gz -o RNA_seq/FastQC_1	
fastqc -f fastq RNA_seq/sh8_1_2.fq.gz -o RNA_seq/FastQC_1	
fastqc -f fastq RNA_seq/sh8_2_1.fq.gz -o RNA_seq/FastQC_1	
fastqc -f fastq RNA_seq/sh8_2_2.fq.gz -o RNA_seq/FastQC_1
multiqc RNA_seq/FastQC_1 -o RNA_seq/FastQC_1	
```

### Trim_Galore
```bash
trim_galore -q 25 --phred33 --length 20 -e 0.1 --stringency 3 --paired RNA_seq/pCDH_1_1.fq.gz RNA_seq/pCDH_1_2.fq.gz --gzip -o RNA_seq	
trim_galore -q 25 --phred33 --length 20 -e 0.1 --stringency 3 --paired RNA_seq/pCDH_2_1.fq.gz RNA_seq/pCDH_2_2.fq.gz --gzip -o RNA_seq
trim_galore -q 25 --phred33 --length 20 -e 0.1 --stringency 3 --paired RNA_seq/A315T_1_1.fq.gz RNA_seq/A315T_1_2.fq.gz --gzip -o RNA_seq	
trim_galore -q 25 --phred33 --length 20 -e 0.1 --stringency 3 --paired RNA_seq/A315T_2_1.fq.gz RNA_seq/A315T_2_2.fq.gz --gzip -o RNA_seq
trim_galore -q 25 --phred33 --length 20 -e 0.1 --stringency 3 --paired RNA_seq/Scrb_1_1.fq.gz RNA_seq/Scrb_1_2.fq.gz --gzip -o RNA_seq	
trim_galore -q 25 --phred33 --length 20 -e 0.1 --stringency 3 --paired RNA_seq/Scrb_2_1.fq.gz RNA_seq/Scrb_2_2.fq.gz --gzip -o RNA_seq	
trim_galore -q 25 --phred33 --length 20 -e 0.1 --stringency 3 --paired RNA_seq/sh2_1_1.fq.gz RNA_seq/sh2_1_2.fq.gz --gzip -o RNA_seq	
trim_galore -q 25 --phred33 --length 20 -e 0.1 --stringency 3 --paired RNA_seq/sh2_2_1.fq.gz RNA_seq/sh2_2_2.fq.gz --gzip -o RNA_seq
trim_galore -q 25 --phred33 --length 20 -e 0.1 --stringency 3 --paired RNA_seq/sh8_1_1.fq.gz RNA_seq/sh8_1_2.fq.gz --gzip -o RNA_seq	
trim_galore -q 25 --phred33 --length 20 -e 0.1 --stringency 3 --paired RNA_seq/sh8_2_1.fq.gz RNA_seq/sh8_2_2.fq.gz --gzip -o RNA_seq
```
### FastQC and MultiQC of Clean Data
```bash
fastqc -f fastq RNA_seq/pCDH_1_1_val_1.fq.gz -o RNA_seq/FastQC_2	
fastqc -f fastq RNA_seq/pCDH_1_2_val_2.fq.gz -o RNA_seq/FastQC_2	
fastqc -f fastq RNA_seq/pCDH_2_1_val_1.fq.gz -o RNA_seq/FastQC_2	
fastqc -f fastq RNA_seq/pCDH_2_2_val_2.fq.gz -o RNA_seq/FastQC_2
fastqc -f fastq RNA_seq/A315T_1_1_val_1.fq.gz -o RNA_seq/FastQC_2	
fastqc -f fastq RNA_seq/A315T_1_2_val_2.fq.gz -o RNA_seq/FastQC_2	
fastqc -f fastq RNA_seq/A315T_2_1_val_1.fq.gz -o RNA_seq/FastQC_2	
fastqc -f fastq RNA_seq/A315T_2_2_val_2.fq.gz -o RNA_seq/FastQC_2
fastqc -f fastq RNA_seq/Scrb_1_1_val_1.fq.gz -o RNA_seq/FastQC_2	
fastqc -f fastq RNA_seq/Scrb_1_2_val_2.fq.gz -o RNA_seq/FastQC_2	
fastqc -f fastq RNA_seq/Scrb_2_1_val_1.fq.gz -o RNA_seq/FastQC_2	
fastqc -f fastq RNA_seq/Scrb_2_2_val_2.fq.gz -o RNA_seq/FastQC_2
fastqc -f fastq RNA_seq/sh2_1_1_val_1.fq.gz -o RNA_seq/FastQC_2	
fastqc -f fastq RNA_seq/sh2_1_2_val_2.fq.gz -o RNA_seq/FastQC_2	
fastqc -f fastq RNA_seq/sh2_2_1_val_1.fq.gz -o RNA_seq/FastQC_2	
fastqc -f fastq RNA_seq/sh2_2_2_val_2.fq.gz -o RNA_seq/FastQC_2
fastqc -f fastq RNA_seq/sh8_1_1_val_1.fq.gz -o RNA_seq/FastQC_2	
fastqc -f fastq RNA_seq/sh8_1_2_val_2.fq.gz -o RNA_seq/FastQC_2	
fastqc -f fastq RNA_seq/sh8_2_1_val_1.fq.gz -o RNA_seq/FastQC_2	
fastqc -f fastq RNA_seq/sh8_2_2_val_2.fq.gz -o RNA_seq/FastQC_2
multiqc RNA_seq/FastQC_2 -o RNA_seq/FastQC_2	
```
