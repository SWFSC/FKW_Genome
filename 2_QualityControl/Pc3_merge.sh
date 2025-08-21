#!/bin/bash

#SBATCH --job-name=Pc3_merge
#SBATCH -e Pc3_merge_%j.e.txt
#SBATCH -o Pc3_merge_%j.log
#SBATCH --mail-user=keith.hernandez@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -D /home/khernandez

#load samtools
module load bio/samtools

#run merge to combine the deduplicated and repeat-masked bams into one file
#Change into the individual subdirectory first
cd bamtools/z0045928

samtools merge Pc_n0025803_sorted_dedup.bam Pc_n0025805_sorted_dedup.bam \
Pc_n0025807_sorted_dedup.bam Pc_n0025809_sorted_dedup.bam Pc_n0025811_sorted_dedup.bam \
Pc_n0025813_sorted_dedup.bam Pc_n0025815_sorted_dedup.bam Pc_n0025817_sorted_dedup.bam \
-o Pcra_z0045928_dedup_noRepeats.bam

#get coverage for individual. This will be the COV variable for heterozygosity and PSMC input
samtools coverage Pcra_z0045928_dedup_noRepeats.bam -o z0045928_coverage.txt
