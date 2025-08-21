#!/bin/bash

#SBATCH --job-name=Pc_nuc_repair3
#SBATCH -e Pc_nuc_repair3_%j.e.txt
#SBATCH -o Pc_nuc_repair3_%j.log
#SBATCH --mail-user=keith.hernandez@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 10
#SBATCH --mem=40G
#SBATCH -t 15:00:00
#SBATCH -D /home/khernandez

#load bbtools
module load bio/bbmap

#Directory path
DIR=/home/khernandez/Pcra_Springer/nuc_reads/z0045928

repair.sh in=${DIR}/n0025803_R1.fastq.gz in2=${DIR}/n0025803_R2.fastq.gz \
out=${DIR}/n0025803_paired_R1.fastq.gz out2=${DIR}/n0025803_paired_R2.fastq.gz outs=${DIR}/n0025803_s.fastq.gz

repair.sh in=${DIR}/n0025805_R1.fastq.gz in2=${DIR}/n0025805_R2.fastq.gz \
out=${DIR}/n0025805_paired_R1.fastq.gz out2=${DIR}/n0025805_paired_R2.fastq.gz outs=${DIR}/n0025805_s.fastq.gz

repair.sh in=${DIR}/n0025807_R1.fastq.gz in2=${DIR}/n0025807_R2.fastq.gz \
out=${DIR}/n0025807_paired_R1.fastq.gz out2=${DIR}/n0025807_R2.fastq.gz outs=${DIR}/n0025807_s.fastq.gz

repair.sh in=${DIR}/n0025809_R1.fastq.gz in2=${DIR}/n0025809_R2.fastq.gz \
out=${DIR}/n0025809_paired_R1.fastq.gz out2=${DIR}/n0025809_paired_R2.fastq.gz outs=${DIR}/n0025809_s.fastq.gz

repair.sh in=${DIR}/n0025811_R1.fastq.gz in2=${DIR}/n0025811_R2.fastq.gz \
out=${DIR}/n0025811_paired_R1.fastq.gz out2=${DIR}/n0025811_paired_R2.fastq.gz outs=${DIR}/n0025811_s.fastq.gz

repair.sh in=${DIR}/n0025813_R1.fastq.gz in2=${DIR}/n0025813_R2.fastq.gz \
out=${DIR}/n0025813_paired_R1.fastq.gz out2=${DIR}/n0025813_paired_R2.fastq.gz outs=${DIR}/n0025813_s.fastq.gz

repair.sh in=${DIR}/n0025815_R1.fastq.gz in2=${DIR}/n0025815_R2.fastq.gz \
out=${DIR}/n0025815_paired_R1.fastq.gz out2=${DIR}/n0025815_paired_R2.fastq.gz outs=${DIR}/n0025815_s.fastq.gz

repair.sh in=${DIR}/n0025817_R1.fastq.gz in2=${DIR}/n0025817_R2.fastq.gz \
out=${DIR}/n0025817_paired_R1.fastq.gz out2=${DIR}/n0025817_paired_R2.fastq.gz outs=${DIR}/n0025817_s.fastq.gz

