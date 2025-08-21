#!/bin/bash

#SBATCH --job-name=Pc_repair3
#SBATCH -e Pc_repair3_%j.e.txt
#SBATCH -o Pc_repair3_%j.log
#SBATCH --mail-user=keith.hernandez@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 10
#SBATCH --mem=40G
#SBATCH -t 15:00:00
#SBATCH -D /home/khernandez

#load bbmap
module load bio/bbmap

#Directory reference
DIR=/home/khernandez/Pcra_Springer/bbduk_out/z0045928

#run repair on paired reads. This can be memory intensive
repair.sh in=${DIR}/n0025803_bbduk_R1.fastq.gz in2=${DIR}/n0025804_bbduk_R2.fastq.gz \
out=${DIR}/n0025803_bbduk_paired_R1.fastq.gz out2=${DIR}/n0025804_bbduk_paired_R2.fastq.gz \
outs=${DIR}/n0025803_s.fastq.gz

repair.sh in=${DIR}/n0025805_bbduk_R1.fastq.gz in2=${DIR}/n0025806_bbduk_R2.fastq.gz \
out=${DIR}/n0025805_bbduk_paired_R1.fastq.gz out2=${DIR}/n0025806_bbduk_paired_R2.fastq.gz \
outs=${DIR}/n0025805_s.fastq.gz

repair.sh in=${DIR}/n0025807_bbduk_R1.fastq.gz in2=${DIR}/n0025808_bbduk_R2.fastq.gz \
out=${DIR}/n0025807_bbduk_paired_R1.fastq.gz out2=${DIR}/n0025808_bbduk_paired_R2.fastq.gz \
outs=${DIR}/n0025807_s.fastq.gz

repair.sh in=${DIR}/n0025809_bbduk_R1.fastq.gz in2=${DIR}/n0025810_bbduk_R2.fastq.gz \
out=${DIR}/n0025809_bbduk_paired_R1.fastq.gz out2=${DIR}/n0025810_bbduk_paired_R2.fastq.gz \
outs=${DIR}/n0025809_s.fastq.gz

repair.sh in=${DIR}/n0025811_bbduk_R1.fastq.gz in2=${DIR}/n0025812_bbduk_R2.fastq.gz \
out=${DIR}/n0025811_bbduk_paired_R1.fastq.gz out2=${DIR}/n0025812_bbduk_paired_R2.fastq.gz \
outs=${DIR}/n0025811_s.fastq.gz

repair.sh in=${DIR}/n0025813_bbduk_R1.fastq.gz in2=${DIR}/n0025814_bbduk_R2.fastq.gz \
out=${DIR}/n0025813_bbduk_paired_R1.fastq.gz out2=${DIR}/n0025814_bbduk_paired_R2.fastq.gz \
outs=${DIR}/n0025813_s.fastq.gz

repair.sh in=${DIR}/n0025815_bbduk_R1.fastq.gz in2=${DIR}/n0025816_bbduk_R2.fastq.gz \
out=${DIR}/n0025815_bbduk_paired_R1.fastq.gz out2=${DIR}/n0025816_bbduk_paired_R2.fastq.gz \
outs=${DIR}/n0025815_s.fastq.gz

repair.sh in=${DIR}/n0025817_bbduk_R1.fastq.gz in2=${DIR}/n0025818_bbduk_R2.fastq.gz \
out=${DIR}/n0025817_bbduk_paired_R1.fastq.gz out2=${DIR}/n0025818_bbduk_paired_R2.fastq.gz \
outs=${DIR}/n0025817_s.fastq.gz

repair.sh in=/home/khernandez/Pcra_Springer/bbduk_out/z0018462/n0026019_bbduk_R1.fastq.gz \
in2=/home/khernandez/Pcra_Springer/bbduk_out/z0018462/n0026020_bbduk_R2.fastq.gz \
out=/home/khernandez/Pcra_Springer/bbduk_out/z0018462/n0026019_bbduk_paired_R1.fastq.gz \
out2=/home/khernandez/Pcra_Springer/bbduk_out/z0018462/n0026020_bbduk_paired_R2.fastq.gz \
outs=/home/khernandez/Pcra_Springer/bbduk_out/z0018462/n0026019_s.fastq.gz

repair.sh in=/home/khernandez/Pcra_Springer/bbduk_out/z0027510/n0025987_bbduk_R1.fastq.gz \
in2=/home/khernandez/Pcra_Springer/bbduk_out/z0027510/n0025988_bbduk_R2.fastq.gz \
out=/home/khernandez/Pcra_Springer/bbduk_out/z0027510/n0025987_bbduk_paired_R1.fastq.gz \
out2=/home/khernandez/Pcra_Springer/bbduk_out/z0027510/n0025988_bbduk_paired_R2.fastq.gz \
outs=/home/khernandez/Pcra_Springer/bbduk_out/z0027510/n0025987_s.fastq.gz

repair.sh in=/home/khernandez/Pcra_Springer/bbduk_out/z0027510/n0025989_bbduk_R1.fastq.gz \
in2=/home/khernandez/Pcra_Springer/bbduk_out/z0027510/n0025990_bbduk_R2.fastq.gz \
out=/home/khernandez/Pcra_Springer/bbduk_out/z0027510/n0025989_bbduk_paired_R1.fastq.gz \
out2=/home/khernandez/Pcra_Springer/bbduk_out/z0027510/n0025990_bbduk_paired_R2.fastq.gz \
outs=/home/khernandez/Pcra_Springer/bbduk_out/z0027510/n0025989_s.fastq.gz
