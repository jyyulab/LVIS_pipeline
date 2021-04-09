#!/bin/bash
#
#BSUB -P insertionSite         # project code
#BSUB -J bam2bed      # job name
#BSUB -W 10:00                # wall-clock time (hrs:mins)
#BSUB -n 1      # number of cpu
#BSUB -R "rusage[mem=8000]"     # memory to reserve, in MB
#BSUB -e errors.%J     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J     # output file name in which %J is replaced by the job ID


module load bedtools/2.25.0

if [ ! -d bam2bed ]; then mkdir bam2bed; fi
for i in `ls bwa | grep -P "bam$" `
do
	out_prefix=`echo $i | sed 's/\.bam//'`
	if [ ! -e bam2bed/${out_prefix}.rmdup.bed ]; then
	samtools view -H bwa/$i > bam2bed/${out_prefix}.sam
	#samtools view  -F 4 -q 1 bwa/$i | perl -ne '$line=$_; @rec=split("\t",$line); $rec[5]=~ s/[SH].*//; if( !($rec[5] =~ /^[0-9]+$/ && $rec[5] > 6)){print $line;}' >> bam2bed/${out_prefix}.sam
	samtools view  -F 4  bwa/$i | perl -ne '$line=$_; @rec=split("\t",$line); $rec[5]=~ s/[SH].*//; if( !($rec[5] =~ /^[0-9]+$/ && $rec[5] > 6)){print $line;}' >> bam2bed/${out_prefix}.sam
	samtools sort -n <(samtools view -hbS bam2bed/${out_prefix}.sam ) bam2bed/${out_prefix}
	bedtools bamtobed -i bam2bed/${out_prefix}.bam -bedpe -mate1 > bam2bed/${out_prefix}.temp
	cat bam2bed/${out_prefix}.temp | perl -e 'while(<STDIN>){$line=$_; @rec=split("\t", $line); if($rec[0] eq $rec[3]){ $start = ($rec[1], $rec[4])[$rec[1] > $rec[4]]; $end = ($rec[2], $rec[5])[$rec[2] < $rec[5]]; print "$rec[0]\t$start\t$end\t$rec[6]\t$rec[7]\t$rec[8]\n";}}' | awk -F "\t" '{if(($3-$2)<1000){print}}' > bam2bed/${out_prefix}.bed
	cut -f 1,2,3,6 bam2bed/${out_prefix}.bed | sort -u | awk '{print $1"\t"$2"\t"$3"\t.\t.\t"$4}' > bam2bed/${out_prefix}.rmdup.bed
	rm bam2bed/${out_prefix}.temp bam2bed/${out_prefix}.bam
	fi
done
