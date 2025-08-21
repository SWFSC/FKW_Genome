#!/bin/bash

#SBATCH --job-name=Pc_mask
#SBATCH -e Pc_mask_%j.e.txt
#SBATCH -o Pc_mask_%j.log
#SBATCH -c 20
#SBATCH -p medmem
#SBATCH -t 20-00:00:00
#SBATCH --mem=180G
#SBATCH --mail-user=keith.hernandez@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -D /scratch/khernandez

#load module
module load bio/repeatmasker/4.1.6
source /opt/bioinformatics/venv/repeatmasker-4.1.6/bin/activate

#define variables
REFDIR=/home/khernandez/Genome_masked
now=$(date +'%d.%m.%Y_%H.%M')

#save copy of repeat for the species indicated
famdb.py -i /opt/bioinformatics/bio/repeatmasker/repeatmasker-4.1.6/Libraries/famdb lineage cetacean -d > ${REFDIR}/cetacea_Dfam3.8_list_${now}.txt

#Run RepeatMasker
RepeatMasker -species cetacea -pa 10 -s GCF_039906515.1_mPseCra1.hap1_genomic.fasta
