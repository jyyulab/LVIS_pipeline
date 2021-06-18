#!/bin/bash
#
#BSUB -P insertionSite         # project code
#BSUB -J peakrefine      # job name
#BSUB -W 10:00                # wall-clock time (hrs:mins)
#BSUB -n 1      # number of cpu
#BSUB -R "rusage[mem=8000]"     # memory to reserve, in MB
#BSUB -e errors.%J     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J     # output file name in which %J is replaced by the job ID


module load bedtools/2.25.0
module load conda3
source activate r_env

for d in `seq 8 8`
do
for i in `ls -1 bam2bed | grep rmdup.bed`
do
        out_prefix=`echo $i | sed 's/\.bed//'`
echo $out_prefix

R --slave <<EOF

library(openxlsx)
library(bedr);

vis_res <- read.xlsx("bed2peak_${d}/${out_prefix}.peak.merge.xlsx" );
res<-read.table("bed2peak_${d}/${out_prefix}.all_peaks.txt");
outfile<-"bed2peak_${d}/${out_prefix}.peak.merge_refine.xlsx";
vis_res<-cbind(vis_res,refine_st=matrix(NA,dim(vis_res)[1],1),refine_ed=matrix(NA,dim(vis_res)[1],1));

aux<-paste0(res[,1],":",res[,2],"-",res[,3]);
res<-res[is.valid.region(aux),];

for (i in 1:dim(vis_res)[1]){
	print(i);
	vis<-paste0(vis_res[i,]$seqnames,":",vis_res[i,]$start,"-",vis_res[i,]$end);
	i_pick<-which(in.region(paste0(res[,1],":",res[,2],"-",res[,3]),vis)&(res[,6]==vis_res[i,]$strand));
	if (length(i_pick)>0){
		im<-which.max(res[i_pick,]$V5);
		vis_res[i,]$refine_st<-res[i_pick[im],]$V2;
		vis_res[i,]$refine_ed<-res[i_pick[im],]$V3;
	}
}
write.xlsx(vis_res,file=outfile);

EOF

done
done
