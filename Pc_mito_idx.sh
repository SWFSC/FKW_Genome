#!/bin/bash

#SBATCH --job-name=Pc_mito_index
#SBATCH -e Pc_mito_index_%j.e.txt
#SBATCH -o Pc_mito_index_%j.log
#SBATCH --mail-user=keith.hernandez@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 1
#SBATCH -D /home/khernandez

#load modules
module load aligners/bwa
module load bio/samtools

bwa index -p /home/khernandez/bwa/Pcra_mitochondrion /home/khernandez/Genomes/Pcra_mitochondrion.fasta

samtools faidx /home/khernandez/Genomes/Pcra_mitochondrion.fasta
