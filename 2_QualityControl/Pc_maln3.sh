#!/bin/bash

#SBATCH --job-name=Pc_mito_align3
#SBATCH -e Pc_ma3_%j.e.txt
#SBATCH -o Pc_ma3_%j.log
#SBATCH --mail-user=keith.hernandez@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task 10
#SBATCH -t 1-00:00:00
#SBATCH -D /home/khernandez
#SBATCH --array=1-8%2

#load modules
module load aligners/bwa
module load bio/samtools

#Define variables
FILE=/home/khernandez/z0045928_mito_align.txt
IDS=$(cat ${FILE})
DIR=/home/khernandez/Pcra_Springer/bbduk_out/z0045928

#create job index
for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	fq_r1=$(echo ${sample_line} | awk -F ":" '{print $2}')
	fq_r2=$(echo ${sample_line} | awk -F ":" '{print $3}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
	break
fi
done

sample_id=$(echo $fq_r1 | sed 's!^.*/!!')
sample_id=${sample_id%%_*}

#align paired reads aginst the mitochondrial genome
bwa mem -t 10 /home/khernandez/bwa/Pcra_mitochondrion ${DIR}/${fq_r1} ${DIR}/${fq_r2} > /home/khernandez/Pcra_Springer/mito_align/z0045928/Malign_${sample_id}.sam

#convert sam to bam and remove the original sam file
samtools view -bS /home/khernandez/Pcra_Springer/mito_align/z0045928/Malign_${sample_id}.sam > /home/khernandez/Pcra_Springer/mito_align/z0045928/Malign_${sample_id}.bam
rm /home/khernandez/Pcra_Springer/mito_align/z0045928/Malign_${sample_id}.sam

#retain unmapped reads (-f 4) from the mitochondrial alignment as nuclear reads
samtools view -f 4 /home/khernandez/Pcra_Springer/mito_align/z0045928/Malign_${sample_id}.bam -o /home/khernandez/Pcra_Springer/mito_align/z0045928/${sample_id}_nuc.bam
rm /home/khernandez/Pcra_Springer/mito_align/z0045928/Malign_${sample_id}.bam
