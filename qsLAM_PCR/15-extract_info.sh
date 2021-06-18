#!/bin/bash
#
#BSUB -P insertionSite         # project code
#BSUB -J extractInfo      # job name
#BSUB -W 40:00                # wall-clock time (hrs:mins)
#BSUB -R "rusage[mem=8000]"     # memory to reserve, in MB
#BSUB -e errors.%J.hybrid     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J.hybrid     # output file name in which %J is replaced by the job ID

if [ ! -d qcInfo ]; then mkdir qcInfo ; fi

for fq_R1 in `ls newFastq/ | grep -P "fastq$" | grep R1`
do
        wc -l newFastq/$fq_R1  > qcInfo/$fq_R1\.qc
	fq_R1_short=`echo $fq_R1 | sed 's/\.fastq//'`
	
	for i in `ls -1 fastqc/beforeCutAdapt | grep fastqc.zip`
	do
        	unzip fastqc/beforeCutAdapt/$i -d fastqc/beforeCutAdapt/
		i2=`echo $i| sed 's/\.zip//'`
		grep "Total Sequences" fastqc/beforeCutAdapt/$i2/fastqc_data.txt >> qcInfo/$fq_R1\.qc
	done
	wc -l bam2bed/$fq_R1_short\.bed >> qcInfo/$fq_R1\.qc
	wc -l bam2bed/$fq_R1_short\.rmdup.bed >> qcInfo/$fq_R1\.qc
	wc -l bed2peak_0/$fq_R1_short\.rmdup.all_peaks.txt >> qcInfo/$fq_R1\.qc
	for d in `seq 0 10`
	do
		wc -l bed2peak_$d/$fq_R1_short\.rmdup.peak.merge.txt >> qcInfo/$fq_R1\.qc
	done
done
