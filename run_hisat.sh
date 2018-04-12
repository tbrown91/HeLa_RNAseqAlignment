#!/bin/bash
# Alignment of HeLa DMSO RNAseq experiments from:
# Molecular basis of differential 3´splice site sensitivity to anti-tumor drugs targeting U2 snRNP. Luisa Vigevani, André Gohr, Thomas Webb, Manuel Irimia and Juan Valcárcel. 

#Align paired-end fastq files using Hisat - good fitting for splced reads
#hisat 2 version 2.0.0-beta

hisat2 -q -k 2 --no-mixed --no-discordant -p 8 -x /home/trin2108/hg38/genome -1 DMSO_tot_1_14444_ACAGTG_read1.fastq -2 DMSO_tot_1_14444_ACAGTG_read2.fastq -S HeLa_RNAseq_DMSOrep1.sam
hisat2 -q -k 2 --no-mixed --no-discordant -p 8 -x /home/trin2108/hg38/genome -1 DMSO_tot_2_13460_ACAGTG_read1.fastq -2 DMSO_tot_2_13460_ACAGTG_read2.fastq -S HeLa_RNAseq_DMSOrep2.sam

#Extract reads above mapping quality (unique and aligns well)
#Split reads by strand
cd ../sam_files/
samtools view -h -q 40 -f 99 -F 3852 -bS HeLa_RNAseq_DMSOrep1.sam > HeLa_RNAseq_DMSOrep1Minus_1.bam &
samtools view -h -q 40 -f 147 -F 3852 -bS HeLa_RNAseq_DMSOrep1.sam > HeLa_RNAseq_DMSOrep1Minus_2.bam &
samtools view -h -q 40 -f 163 -F 3852 -bS HeLa_RNAseq_DMSOrep1.sam > HeLa_RNAseq_DMSOrep1Plus_1.bam &
samtools view -h -q 40 -f 83 -F 3852 -bS HeLa_RNAseq_DMSOrep1.sam > HeLa_RNAseq_DMSOrep1Plus_2.bam &
samtools view -h -q 40 -f 99 -F 3852 -bS HeLa_RNAseq_DMSOrep2.sam > HeLa_RNAseq_DMSOrep2Minus_1.bam &
samtools view -h -q 40 -f 147 -F 3852 -bS HeLa_RNAseq_DMSOrep2.sam > HeLa_RNAseq_DMSOrep2Minus_2.bam &
samtools view -h -q 40 -f 163 -F 3852 -bS HeLa_RNAseq_DMSOrep2.sam > HeLa_RNAseq_DMSOrep2Plus_1.bam &
samtools view -h -q 40 -f 83 -F 3852 -bS HeLa_RNAseq_DMSOrep2.sam > HeLa_RNAseq_DMSOrep2Plus_2.bam 

#Combine pairs of reads from same strand
samtools merge HeLa_RNAseq_DMSOrep1Plus_merge.bam HeLa_RNAseq_DMSOrep1Plus_1.bam HeLa_RNAseq_DMSOrep1Plus_1.bam &
samtools merge HeLa_RNAseq_DMSOrep1Minus_merge.bam HeLa_RNAseq_DMSOrep1Minus_1.bam HeLa_RNAseq_DMSOrep1Minus_1.bam &
samtools merge HeLa_RNAseq_DMSOrep2Plus_merge.bam HeLa_RNAseq_DMSOrep2Plus_1.bam HeLa_RNAseq_DMSOrep2Plus_1.bam &
samtools merge HeLa_RNAseq_DMSOrep2Minus_merge.bam HeLa_RNAseq_DMSOrep2Minus_1.bam HeLa_RNAseq_DMSOrep2Minus_1.bam

#Sort bam files and index
samtools sort HeLa_RNAseq_DMSOrep1Plus_merge.bam HeLa_RNAseq_DMSOrep1Plus_sorted &
samtools sort HeLa_RNAseq_DMSOrep1Minus_merge.bam HeLa_RNAseq_DMSOrep1Minus_sorted &
samtools sort HeLa_RNAseq_DMSOrep2Plus_merge.bam HeLa_RNAseq_DMSOrep2Plus_sorted &
samtools sort HeLa_RNAseq_DMSOrep2Minus_merge.bam HeLa_RNAseq_DMSOrep2Minus_sorted

samtools index HeLa_RNAseq_DMSOrep1Plus_sorted.bam &
samtools index HeLa_RNAseq_DMSOrep1Minus_sorted.bam &
samtools index HeLa_RNAseq_DMSOrep2Plus_sorted.bam &
samtools index HeLa_RNAseq_DMSOrep2Minus_sorted.bam &

#Split files by chromosome to be read by R script

bamtools split -in rep1_plus/HeLa_RNAseq_DMSOrep1Plus_sorted.bam -reference &
bamtools split -in rep1_minus/HeLa_RNAseq_DMSOrep1Minus_sorted.bam -reference &
bamtools split -in rep2_plus/HeLa_RNAseq_DMSOrep2Plus_sorted.bam -reference &
bamtools split -in rep2_minus/HeLa_RNAseq_DMSOrep2Minus_sorted.bam -reference &

