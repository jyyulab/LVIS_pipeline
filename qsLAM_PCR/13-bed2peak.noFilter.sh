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
module load R/3.5.1

for d in 8
do
if [ ! -d bed2peak_${d} ]; then mkdir bed2peak_${d}; fi
if [ ! -d bed2peak_filtered_${d} ]; then mkdir bed2peak_filtered_${d}; fi
for i in `ls -1 bam2bed | grep rmdup.bed`
do
        out_prefix=`echo $i | sed 's/\.bed//'`
echo $out_prefix


#out_prefix=1555624_A10W027_11Jul18_PB_CD3_S10_L001_R1_001.rmdup
		
#	awk -F "\t" '{if($6=="+"){print $1"\t"$2"\t"$6} else {print $1"\t"$3"\t"$6}}' bam2bed/${out_prefix}.bed | sort | uniq -c | sort -k2,2 -k3,3n > bed2peak/${out_prefix}.peak.xls
	awk -F "\t" '{if($6=="+"){print $1"\t"$2"\t"$6} else {print $1"\t"$3"\t"$6}}' bam2bed/${out_prefix}.bed | grep -v vector | sort | uniq -c | awk '{print $2"\t"$3"\t"$3+1"\t.\t"$1"\t"$4}' | sort -k1,1 -k2,2n | bedtools merge -c 5 -s -o sum -d $d -i - | awk '{print $1"\t"$2"\t"$3"\t.\t"$5"\t"$4}' | grep -v vector | bedtools intersect -b <(grep -v vector bam2bed/${out_prefix%.rmdup}.bed| awk '{if($6=="+"){print $1"\t"$2"\t"$2+1"\t.\t.\t"$6} else {print $1"\t"$3"\t"$3+1"\t.\t.\t"$6}}') -c -s -a - | awk '{print $1"\t"$2"\t"$3"\t""'$out_prefix'""\t"$5"\t"$6"\t"$7}' > bed2peak_${d}/${out_prefix}.peak.merge.xls

	bedtools window -w $d -a bed2peak_${d}/${out_prefix}.peak.merge.xls -b known20sites.bed -sm -c > bed2peak_${d}/${out_prefix}.2.xls 
R --slave <<EOF

source("target_gene_prediction.R")
inputBed <- read.table("bed2peak_${d}/${out_prefix}.2.xls", sep="\t", check.names =FALSE )
colnames(inputBed) <- c("seqnames", "start", "end", "name", "nUniqueReads", "strand", "nReads", "nOverlapWithSpike")
inputBed[["name"]] <- gsub("-", ".", paste("X", inputBed[["name"]], 1:nrow(inputBed), sep="_"))
d <- target_gene_prediction(inputBed)

inputBed=inputBed[order(inputBed$name),]


rownames(d) <- d[["peak_id"]]
inputBed[["gene"]] <- d[inputBed[["name"]], "nearest_gene_symbol"]
inputBed[["tss_distances"]] <- d[inputBed[["name"]], "tss_distances"]
inputBed[["gene_region"]] <- d[inputBed[["name"]], "gene_region"]
inputBed[["name"]] <- sub("^X_", "", inputBed[["name"]])
inputBed <- inputBed[order(inputBed[["nUniqueReads"]], decreasing=TRUE),]
write.table(inputBed, file="bed2peak_${d}/${out_prefix}.peak.merge.xls", sep="\t", quote=FALSE, row.names=FALSE)
EOF


R --slave <<EOF
library("xlsx")
options(stringsAsFactors = FALSE)
inputBed <- read.table("bed2peak_${d}/${out_prefix}.peak.merge.xls", sep="\t", header=TRUE, check.names =FALSE )
inputBed <- subset(inputBed, nUniqueReads > 1 | nReads > 5 )
inputBed[["percent"]] <- inputBed[,"nUniqueReads"] / sum(inputBed[,"nUniqueReads"])
inputBed[["percent"]] <- round(100*inputBed[["percent"]], 2)

inputBed[["percentAllReads"]] <- inputBed[,"nReads"] / sum(inputBed[,"nReads"])
inputBed[["percentAllReads"]] <- round(100*inputBed[["percentAllReads"]], 2)

write.xlsx(inputBed, file="bed2peak_${d}/${out_prefix}.peak.merge.xlsx")
#inputBed <- subset(inputBed, nOverlapWithSpike == 0)
inputBed[["percent"]] <- inputBed[,"nUniqueReads"] / sum(inputBed[,"nUniqueReads"])
inputBed[["percent"]] <- round(100*inputBed[["percent"]], 2)

inputBed[["percentAllReads"]] <- inputBed[,"nReads"] / sum(inputBed[,"nReads"])
inputBed[["percentAllReads"]] <- round(100*inputBed[["percentAllReads"]], 2)

write.xlsx(inputBed, file="bed2peak_filtered_${d}/${out_prefix}.peak.merge.xlsx")
write.xlsx(head(inputBed, n=20), file="bed2peak_filtered_${d}/${out_prefix}.top20.peak.merge.xlsx")
EOF

#rm bed2peak_${d}/${out_prefix}.2.xls 
#rm bed2peak_${d}/${out_prefix}.peak.merge.xls
done
done
