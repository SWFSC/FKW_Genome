#!/bin/bash

#SBATCH --job-name=bam2fq3
#SBATCH -e bam2fq3_%j.e.txt
#SBATCH -o bam2fq3_%j.log
#SBATCH --mail-user=keith.hernandez@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task 10
#SBATCH -t 10:00:00
#SBATCH -D /home/khernandez
#SBATCH --array=1-8%2

#load samtools
module load bio/samtools

#Define variables
FILE=/home/khernandez/z0045928_bams.txt
IDS=$(cat ${FILE})
INDIR=/home/khernandez/Pcra_Springer/mito_align/z0045928
OUTDIR=/home/khernandez/Pcra_Springer/nuc_reads/z0045928

#create job index
for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	sample=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
	break
fi
done

sample_id=$(echo $sample | sed 's!^.*/!!')
sample_id=${sample%%_*}

#files have to be sorted first before R1 and R2 can be extracted. Pipe sort into the fastq command
#-n in sort command sorts by name. -n in fastq does not append /1 or 2/ to respective read names
samtools sort -n ${INDIR}/${sample_id}_nuc.bam | \
samtools fastq -1 ${OUTDIR}/${sample_id}_R1.fastq.gz -2 ${OUTDIR}/${sample_id}_R2.fastq.gz -n
