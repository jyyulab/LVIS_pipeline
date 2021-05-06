#!/bin/bash
#
#BSUB -P insertionSite         # project code
#BSUB -J cutPrimer      # job name
#BSUB -W 40:00                # wall-clock time (hrs:mins)
#BSUB -R "rusage[mem=8000]"     # memory to reserve, in MB
#BSUB -e errors.%J.hybrid     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J.hybrid     # output file name in which %J is replaced by the job ID

module load python/3.6.1

cutadapt=/hpcf/apps/python/install/3.6.1/bin/cutadapt


if [ ! -d cutPrimer ]; then mkdir cutPrimer; fi
for f in `ls -1 rawdata | grep -P "fastq.gz$" | grep R1`
do
	if [ ! -e cutPrimer/${f%.gz} ]; then
	zcat rawdata/$f | $cutadapt -e 0.1 -g ATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTC --info-file cutPrimer/$f.R1 - > cutPrimer/${f%.gz} 
	#zcat rawdata/$f | cutadapt -e 0.1 -g ATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTC --trimmed-only -O 10 -m 20 - > cutPrimer/${f%.gz} 
	fi
done
for f in `ls -1 rawdata | grep -P "fastq.gz$" | grep R2`
do
	if [ ! -e cutPrimer/${f%.gz} ]; then
	zcat rawdata/$f | $cutadapt -e 0.1 -g GACTGCGTATCAGT --info-file cutPrimer/$f.R2 - > cutPrimer/${f%.gz} 
	#zcat rawdata/$f | cutadapt -e 0.1 -g GACTGCGTATCAGT --trimmed-only -O 10  -m 20  - > cutPrimer/${f%.gz} 
	fi
done

cd cutPrimer
rename fastq.gz fastq *fastq.gz*


# adaptor sequences from: http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences_1000000002694-01.pdf
