#!/bin/bash

#SBATCH --job-name=RM_align3
#SBATCH -e RM_align3_%j.e.txt
#SBATCH -o RM_align3_%j.log
#SBATCH --mail-user=keith.hernandez@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task 10
#SBATCH -t 1-00:00:00
#SBATCH -D /home/khernandez
#SBATCH --array=1-8%2

#load modules
module load aligners/bwa
module load bio/samtools
module load bio/bamtools
module load bio/picard

#Define variables
FILE=/home/khernandez/z0045928_align.txt
IDS=$(cat ${FILE})
DIR=/home/khernandez/Pcra_Springer/nuc_reads/z0045928

#create job index to run through array and paired files
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

#align paired samples to the repeat-masked genome and create sam file
bwa mem -M -t 10 /home/khernandez/bwa_masked/GCF_039906515.1_mPseCra1.hap1_genomic ${DIR}/${fq_r1} ${DIR}/${fq_r2} > /home/khernandez/bamtools/Pc_${sample_id}.sam

#convert sam to bam
samtools view -bS -F 4 /home/khernandez/bamtools/Pc_${sample_id}.sam > /home/khernandez/bamtools/Pc_${sample_id}.bam
rm /home/khernandez/bamtools/Pc_${sample_id}.sam

#sort bam
samtools view -h /home/khernandez/bamtools/Pc_${sample_id}.bam | samtools view -buS | \
samtools sort -o /home/khernandez/bamtools/Pc_${sample_id}_sorted.bam
rm /home/khernandez/bamtools/Pc_${sample_id}.bam

#Use Picard to remove duplicates
java -jar ${PICARD} MarkDuplicates I=/home/khernandez/bamtools/Pc_${sample_id}_sorted.bam \
O=/home/khernandez/bamtools/Pc_${sample_id}_sorted_dedup.bam \
M=/home/khernandez/bamtools/Pc_${sample_id}_dups.log VALIDATION_STRINGENCY=SILENT \
REMOVE_DUPLICATES=true

#calculate depth of coverage
samtools depth -aa /home/khernandez/bamtools/Pc_${sample_id}_sorted_dedup.bam | cut -f 3 | \
gzip > /home/khernandez/bamtools/Pc_${sample_id}.depth.gz

#index the aligned bam files
samtools index /home/khernandez/bamtools/Pc_${sample_id}_sorted_dedup.bam
