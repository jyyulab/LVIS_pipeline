#!/bin/bash
#
#BSUB -P insertionSite         # project code
#BSUB -J bwa      # job name
#BSUB -W 10:00                # wall-clock time (hrs:mins)
#BSUB -n 4	# number of cpu
#BSUB -R "rusage[mem=8000]"     # memory to reserve, in MB
#BSUB -e errors.%J     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J     # output file name in which %J is replaced by the job ID

module load bwa
module load samtools


if [ ! -d bwa ]; then mkdir bwa; fi 
for fq_R1 in `ls newFastq/ | grep -P "fastq$" | grep R1`
do
	fq_R2=`echo $fq_R1| sed 's/R1/R2/'`
	bwa_out=`echo $fq_R1| sed 's/fastq/bam/' `
	if [ ! -e bwa/$bwa_out ]; then 
	bwa mem -t 4 bwa_index/hg19_noHap.fa.gz newFastq/$fq_R1 newFastq/$fq_R2 | samtools view -hbS - > bwa/$bwa_out
	fi
done
