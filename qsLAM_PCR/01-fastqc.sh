#!/bin/bash
#
#BSUB -P insertionSite         # project code
#BSUB -J fastqc      # job name
#BSUB -W 10:00                # wall-clock time (hrs:mins)
#BSUB -n 1      # number of cpu
#BSUB -R "rusage[mem=8000]"     # memory to reserve, in MB
#BSUB -e errors.%J     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J     # output file name in which %J is replaced by the job ID



module load fastqc/0.11.2

if [ ! -d fastqc/beforeCutAdapt ]; then mkdir -p fastqc/beforeCutAdapt; fi 
for i in `ls -1 rawdata/ | grep fastq.gz`
do
	if [ ! -e fastqc/beforeCutAdapt/${i}_fastqc.html ]; then 
	fastqc -o fastqc/beforeCutAdapt rawdata/$i
	fi
done
