#!/bin/bash

#SBATCH --job-name=bbduk_Pc3
#SBATCH -e bbduk_Pc3_%j.e.txt
#SBATCH -o bbduk_Pc3_%j.log
#SBATCH --mail-user=keith.hernandez@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task 10
#SBATCH -t 10:00:00
#SBATCH -D /home/khernandez

#load bbtools module
module load bio/bbmap

#Directory references
INDIR=/home/khernandez/Pcra_Springer/Raw_reads/z0045928
OUTDIR=/home/khernandez/Pcra_Springer/bbduk_out/z0045928
REF=/opt/bioinformatics/bio/bbmap/bbmap-39.06/resources

#run bbduk.sh on paired files. This removes adapters, trims ends to quality >15, min. avg. quality > 20 and length>40 bp
bbduk.sh in1=${INDIR}/z0045928_Pcra_z0045928.Pcra.Springer2017_n0025803_R1.fastq.gz \
in2=${INDIR}/z0045928_Pcra_z0045928.Pcra.Springer2017_n0025804_R2.fastq.gz out1=${OUTDIR}/n0025803_bbduk_R1.fastq.gz \
out2=${OUTDIR}/n0025804_bbduk_R2.fastq.gz ref=${REF}/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo \
ftl=0 qtrim=rl trimq=15 maq=20 minlen=40

bbduk.sh in1=${INDIR}/z0045928_Pcra_z0045928.Pcra.Springer2017_n0025805_R1.fastq.gz \
in2=${INDIR}/z0045928_Pcra_z0045928.Pcra.Springer2017_n0025806_R2.fastq.gz out1=${OUTDIR}/n0025805_bbduk_R1.fastq.gz \
out2=${OUTDIR}/n0025806_bbduk_R2.fastq.gz ref=${REF}/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo \
ftl=0 qtrim=rl trimq=15 maq=20 minlen=40

bbduk.sh in1=${INDIR}/z0045928_Pcra_z0045928.Pcra.Springer2017_n0025807_R1.fastq.gz \
in2=${INDIR}/z0045928_Pcra_z0045928.Pcra.Springer2017_n0025808_R2.fastq.gz out1=${OUTDIR}/n0025807_bbduk_R1.fastq.gz \
out2=${OUTDIR}/n0025808_bbduk_R2.fastq.gz ref=${REF}/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo \
ftl=0 qtrim=rl trimq=15 maq=20 minlen=40

bbduk.sh in1=${INDIR}/z0045928_Pcra_z0045928.Pcra.Springer2017_n0025809_R1.fastq.gz \
in2=${INDIR}/z0045928_Pcra_z0045928.Pcra.Springer2017_n0025810_R2.fastq.gz out1=${OUTDIR}/n0025809_bbduk_R1.fastq.gz \
out2=${OUTDIR}/n0025810_bbduk_R2.fastq.gz ref=${REF}/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo \
ftl=0 qtrim=rl trimq=15 maq=20 minlen=40

bbduk.sh in1=${INDIR}/z0045928_Pcra_z0045928.Pcra.Springer2017_n0025811_R1.fastq.gz \
in2=${INDIR}/z0045928_Pcra_z0045928.Pcra.Springer2017_n0025812_R2.fastq.gz out1=${OUTDIR}/n0025811_bbduk_R1.fastq.gz \
out2=${OUTDIR}/n0025812_bbduk_R2.fastq.gz ref=${REF}/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo \
ftl=0 qtrim=rl trimq=15 maq=20 minlen=40

bbduk.sh in1=${INDIR}/z0045928_Pcra_z0045928.Pcra.Springer2017_n0025813_R1.fastq.gz \
in2=${INDIR}/z0045928_Pcra_z0045928.Pcra.Springer2017_n0025814_R2.fastq.gz out1=${OUTDIR}/n0025813_bbduk_R1.fastq.gz \
out2=${OUTDIR}/n0025814_bbduk_R2.fastq.gz ref=${REF}/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo \
ftl=0 qtrim=rl trimq=15 maq=20 minlen=40

bbduk.sh in1=${INDIR}/z0045928_Pcra_z0045928.Pcra.Springer2017_n0025815_R1.fastq.gz \
in2=${INDIR}/z0045928_Pcra_z0045928.Pcra.Springer2017_n0025816_R2.fastq.gz out1=${OUTDIR}/n0025815_bbduk_R1.fastq.gz \
out2=${OUTDIR}/n0025816_bbduk_R2.fastq.gz ref=${REF}/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo \
ftl=0 qtrim=rl trimq=15 maq=20 minlen=40

bbduk.sh in1=${INDIR}/z0045928_Pcra_z0045928.Pcra.Springer2017_n0025817_R1.fastq.gz \
in2=${INDIR}/z0045928_Pcra_z0045928.Pcra.Springer2017_n0025818_R2.fastq.gz out1=${OUTDIR}/n0025817_bbduk_R1.fastq.gz \
out2=${OUTDIR}/n0025818_bbduk_R2.fastq.gz ref=${REF}/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo \
ftl=0 qtrim=rl trimq=15 maq=20 minlen=40

bbduk.sh in1=/home/khernandez/Pcra_Springer/Raw_reads/z0018462/z0018462_Pcra_z0018462.Pcra.Springer2018_n0026019_R1.fastq.gz \
in2=/home/khernandez/Pcra_Springer/Raw_reads/z0018462/z0018462_Pcra_z0018462.Pcra.Springer2018_n0026020_R2.fastq.gz \
out1=/home/khernandez/Pcra_Springer/bbduk_out/z0018462/n0026019_bbduk_R1.fastq.gz \
out2=/home/khernandez/Pcra_Springer/bbduk_out/z0018462/n0026020_bbduk_R2.fastq.gz ref=${REF}/adapters.fa ktrim=r k=23 \
mink=11 hdist=1 tbo ftl=0 qtrim=rl trimq=15 maq=20 minlen=40
