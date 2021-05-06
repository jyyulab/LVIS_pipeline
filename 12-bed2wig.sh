#!/bin/bash
#
#BSUB -P insertionSite         # project code
#BSUB -J bed2wig      # job name
#BSUB -W 10:00                # wall-clock time (hrs:mins)
#BSUB -n 1	# number of cpu
#BSUB -R "rusage[mem=8000]"     # memory to reserve, in MB
#BSUB -e errors.%J     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J     # output file name in which %J is replaced by the job ID


module load bedtools

if [ ! -d bed2wig ]; then mkdir bed2wig; fi
for i in `ls -1 bam2bed | grep -P ".rmdup.bed$"`
do
	out_prefix=`echo $i | sed 's/.bed//'`

	if [ ! -e bed2wig/${out_prefix}.pos.bw ]; then
	genomeCoverageBed -i <(sort -k 1,1 -k 2,2n bam2bed/$i) -bg -g hg19.chrom.sizes -strand + | sort -k1,1 -k2,2n > bed2wig/${out_prefix}.pos.bedgraph
	 ./bedGraphToBigWig bed2wig/${out_prefix}.pos.bedgraph hg19.chrom.sizes bed2wig/${out_prefix}.pos.bw
	rm bed2wig/${out_prefix}.pos.bedgraph

	genomeCoverageBed -i <(sort -k 1,1 -k 2,2n bam2bed/$i) -bg -g hg19.chrom.sizes -strand - | sort -k1,1 -k2,2n > bed2wig/${out_prefix}.neg.bedgraph
	 ./bedGraphToBigWig bed2wig/${out_prefix}.neg.bedgraph hg19.chrom.sizes bed2wig/${out_prefix}.neg.bw
	rm bed2wig/${out_prefix}.neg.bedgraph
	fi
done
