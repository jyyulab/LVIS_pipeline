#!/bin/bash
#
#BSUB -P insertionSite         # project code
#BSUB -J bed2peak      # job name
#BSUB -W 10:00                # wall-clock time (hrs:mins)
#BSUB -n 1      # number of cpu
#BSUB -R "rusage[mem=8000]"     # memory to reserve, in MB
#BSUB -e errors.%J     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J     # output file name in which %J is replaced by the job ID


module load bedtools
#module load R/3.5.1

for d in `seq 0 10`
do
if [ ! -d bed2peak_${d} ]; then mkdir bed2peak_${d}; fi
for i in `ls -1 bam2bed | grep rmdup.bed`
do
        out_prefix=`echo $i | sed 's/\.bed//'`
echo $out_prefix


#out_prefix=1555624_A10W027_11Jul18_PB_CD3_S10_L001_R1_001.rmdup
		
#	awk -F "\t" '{if($6=="+"){print $1"\t"$2"\t"$6} else {print $1"\t"$3"\t"$6}}' bam2bed/${out_prefix}.bed | sort | uniq -c | sort -k2,2 -k3,3n > bed2peak/${out_prefix}.peak.xls
#	awk -F "\t" '{if($6=="+"){print $1"\t"$2"\t"$6} else {print $1"\t"$3"\t"$6}}' bam2bed/${out_prefix}.bed | grep -v vector | sort | uniq -c | awk '{print $2"\t"$3"\t"$3+1"\t.\t"$1"\t"$4}' | sort -k1,1 -k2,2n > bed2peak_${d}/${out_prefix}.peak.txt
	
	awk -F "\t" '{if($6=="+"){print $1"\t"$2"\t"$6} else {print $1"\t"$3"\t"$6}}' bam2bed/${out_prefix}.bed | grep -v vector | sort | uniq -c | awk '{print $2"\t"$3"\t"$3+1"\t.\t"$1"\t"$4}' | sort -k1,1 -k2,2n | bedtools merge -c 5 -s -o sum -d $d -i - | awk '{print $1"\t"$2"\t"$3"\t.\t"$5"\t"$4}' | grep -v vector | bedtools intersect -b <(grep -v vector bam2bed/${out_prefix%.rmdup}.bed| awk '{if($6=="+"){print $1"\t"$2"\t"$2+1"\t.\t.\t"$6} else {print $1"\t"$3"\t"$3+1"\t.\t.\t"$6}}') -c -s -a - | awk '{print $1"\t"$2"\t"$3"\t""'$out_prefix'""\t"$5"\t"$6"\t"$7}' > bed2peak_${d}/${out_prefix}.peak.merge.txt

	awk -F "\t" '{if($6=="+"){print $1"\t"$2"\t"$6} else {print $1"\t"$3"\t"$6}}' bam2bed/${out_prefix}.bed | grep -v vector | sort | uniq -c | awk '{print $2"\t"$3"\t"$3+1"\t.\t"$1"\t"$4}' | sort -k1,1 -k2,2n | grep -v vector | bedtools intersect -b <(grep -v vector bam2bed/${out_prefix%.rmdup}.bed| awk '{if($6=="+"){print $1"\t"$2"\t"$2+1"\t.\t.\t"$6} else {print $1"\t"$3"\t"$3+1"\t.\t.\t"$6}}') -c -s -a - | awk '{print $1"\t"$2"\t"$3"\t""'$out_prefix'""\t"$5"\t"$6"\t"$7}' > bed2peak_${d}/${out_prefix}.all_peaks.txt

	bedtools window -w $d -a bed2peak_${d}/${out_prefix}.peak.merge.txt -b known20sites.bed -sm -c > bed2peak_${d}/${out_prefix}.peak.merge.2.txt

done
done
