#!/bin/bash
#
#BSUB -P insertionSite         # project code
#BSUB -J makeNewFastQ      # job name
#BSUB -W 40:00                # wall-clock time (hrs:mins)
#BSUB -R "rusage[mem=8000]"     # memory to reserve, in MB
#BSUB -e errors.%J.hybrid     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J.hybrid     # output file name in which %J is replaced by the job ID

if [ ! -d newFastq ]; then mkdir newFastq ; fi
for fq_R1 in `ls cutPrimer/ | grep -P "fastq$" | grep R1`
do
	if [ ! -e newFastq/$fq_R1 ]; then
	fq_R2=`echo $fq_R1| sed 's/R1/R2/'`
	paste <(cat cutPrimer/$fq_R1 | paste - - - -| sed 's/^.//') <(cat cutPrimer/$fq_R2 | paste - - - -|sed 's/^.//') > temp.fastq
	./myjoin temp.fastq cutPrimer/$fq_R1.R1 | grep -P "^=" | cut -f 2,3,4,5,6,7,8,9,12,13 | awk -F "\t"  '{if($10-$9>30){print}}' | cut -f 1-4 | sed 's/^/@/'| sed 's/\t/\n/g' > newFastq/$fq_R1 
	./myjoin temp.fastq cutPrimer/$fq_R1.R1 | grep -P "^=" | cut -f 2,3,4,5,6,7,8,9,12,13 | awk -F "\t"  '{if($10-$9>30){print}}' | cut -f 5-8 | sed 's/^/@/'| sed 's/\t/\n/g' > newFastq/$fq_R2 
	rm temp.fastq
	fi
done

#cd cutPrimer
#rename fastq.gz fastq *fastq.gz*
