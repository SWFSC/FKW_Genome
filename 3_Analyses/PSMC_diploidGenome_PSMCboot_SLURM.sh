#!/bin/bash

#SBATCH --job-name=PSMC1
#SBATCH -e PSMC1_%j.e.txt
#SBATCH -o PSMC1_%j.log
#SBATCH --mail-user=keith.hernandez@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH --mem=80G
#SBATCH -t 3-0
#SBATCH -D /scratch/khernandez/temp

#programs to load
module load bio/samtools/1.11
module load bio/bcftools/1.11
module load tools/gnuplot/5.4.1
module load bio/bbmap
PSMC=/opt/bioinformatics/bio/psmc/psmc-0.6.5
set -eux

########################################
#variables for bam and reference files
BAMDIR=/home/khernandez/bamtools/z0018462
BAMFILE=Pcra_z0018462_dedup_noRepeats.bam
REFDIR=/home/khernandez/Genomes
REF=GCF_039906515.1_mPseCra1.hap1_genomic.fna
OUTDIR=${BAMDIR}/PSMC
XYFILE=/home/khernandez/Pcra_XY.txt

#Define variables for PSMC 
#generation time for Pcra, average depth of coverage and number of bootstraps to run in parallel
GEN=24.0
COV=34
THREADS=50

#Extract portions of the BAM filename for output files
ID=`echo $BAMFILE | cut -f1,2 -d "_"`

DIPLOID=${ID}_diploid.fq.gz

mkdir -p ${OUTDIR}

#######################################
##PSMC Parameters
TEMP_DIR=/scratch/khernandez/temp
MUT=4.90E-10 #mutation rate (mu/site/yr)
MUTRATE=$(echo $GEN $MUT | awk '{ printf "%.12f", $1*$2 }') #mutation rate by generation time
t=15 #default from Li and Durbin
PSMC_INT="4+25*2+4+6" #default from Li and Durbin

#diploid genome parameters
MEANDEPTH=${COV} #has to be an interger for min and max calculations below
let "MIN=MEANDEPTH/3" #1/3 average coverage
let "MAX=MEANDEPTH*2" #2x average coverage

########################################
#Create a diploid consensus fastq file from the bam file. May take >24hr
samtools mpileup -uf ${REFDIR}/${REF} ${BAMDIR}/${BAMFILE} | bcftools call -c --threads ${THREADS} | vcfutils.pl vcf2fq -d ${MIN} -D ${MAX} > ${OUTDIR}/${DIPLOID}

########################################
#PSMC
#generate psmcfa file from diploid genome file
${PSMC}/utils/fq2psmcfa -q20 ${OUTDIR}/${DIPLOID} > ${OUTDIR}/${ID}_diploid.psmcfa

#Remove X and Y chromosomes before running PSMC, following Carroll et al. for M. mirus/M. eueu
filterbyname.sh in=${OUTDIR}/${ID}_diploid.psmcfa out=${OUTDIR}/${ID}_noXY_diploid.psmcfa names=${XYFILE}

#generate the psmc file using the default settings for humans (-N25, -t20, -r5)
${PSMC}/psmc -N25 -t${t} -r5 -p ${PSMC_INT} -o ${OUTDIR}/${ID}_${MUT}_t${t}.psmc ${OUTDIR}/${ID}_noXY_diploid.psmcfa

#Make the PSMC plots and adapt the scaling using this PSMC file
nice ${PSMC}/utils/psmc_plot.pl -u ${MUTRATE} -g ${GEN} -RM ${ID}"_"${MUT} ${OUTDIR}/${ID}_${MUT}_t${t}_psmc.out ${OUTDIR}/${ID}_${MUT}_t${t}.psmc

#######################################
### bootstrap PSMC
#Split the PSMC
${PSMC}/utils/splitfa ${OUTDIR}/${ID}_noXY_diploid.psmcfa > ${OUTDIR}/split_${ID}_diploid.psmcfa

#PSMC bootstrap, multithread (written to temp directory)
seq 100 | xargs -P ${THREADS} -i ${PSMC}/psmc -N25 -t${t} -r5 -b -p ${PSMC_INT} -o ${TEMP_DIR}/${ID}_round-{}.psmc ${OUTDIR}/split_${ID}_diploid.psmcfa | sh

#merge original.psmc round-*.psmc > merged.psmc
#within folder containing original sample psmc and bootstrap files
cat ${OUTDIR}/${ID}_${MUT}_t${t}.psmc ${TEMP_DIR}/${ID}_round-*.psmc > ${OUTDIR}/merged_${ID}_t${t}_boot.psmc

#and then plot it: (- -> ID placed on plot)
${PSMC}/utils/psmc_plot.pl -u ${MUTRATE} -g ${GEN} -RM "" ${OUTDIR}/merged_${ID}_t${t}_boot.out ${OUTDIR}/merged_${ID}_t${t}_boot.psmc

#cp ${BAMDIR}/*PSMCboot.SEDNA.sh ${OUTDIR}
