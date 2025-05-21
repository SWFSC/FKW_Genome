#!/bin/bash

#SBATCH --job-name=RM_index
#SBATCH -e RM_index_%j.e.txt
#SBATCH -o RM_index_%j.log
#SBATCH --mail-user=keith.hernandez@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 1
#SBATCH -D /home/khernandez

#load modules
module load aligners/bwa
module load bio/samtools

bwa index -p /home/khernandez/bwa_masked/GCF_039906515.1_mPseCra1.hap1_genomic /home/khernandez/Genomes/GCF_039906515.1_mPseCra1.hap1_genomic.fasta.masked

samtools faidx /home/khernandez/Genomes/GCF_039906515.1_mPseCra1.hap1_genomic.fasta.masked
