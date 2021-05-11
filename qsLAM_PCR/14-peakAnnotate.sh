#!/bin/bash
#
#BSUB -P insertionSite         # project code
#BSUB -J bed2peak      # job name
#BSUB -W 10:00                # wall-clock time (hrs:mins)
#BSUB -n 1      # number of cpu
#BSUB -R "rusage[mem=8000]"     # memory to reserve, in MB
#BSUB -e errors.%J     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J     # output file name in which %J is replaced by the job ID


module load bedtools/2.25.0
module load R/3.5.1

for d in `seq 0 10`
do
for i in `ls -1 bam2bed | grep rmdup.bed`
do
        out_prefix=`echo $i | sed 's/\.bed//'`
echo $out_prefix

R --slave <<EOF

source("target_gene_prediction.R")
inputBed <- read.table("bed2peak_${d}/${out_prefix}.peak.merge.2.txt", sep="\t", check.names =FALSE )
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
write.table(inputBed, file="bed2peak_${d}/${out_prefix}.peak.merge.3.txt", sep="\t", quote=FALSE, row.names=FALSE)
EOF


R --slave <<EOF
library("xlsx")
options(stringsAsFactors = FALSE)
inputBed <- read.table("bed2peak_${d}/${out_prefix}.peak.merge.3.txt", sep="\t", header=TRUE, check.names =FALSE )
inputBed <- subset(inputBed, nUniqueReads > 1 | nReads > 5 )
inputBed[["percent"]] <- inputBed[,"nUniqueReads"] / sum(inputBed[,"nUniqueReads"])
inputBed[["percent"]] <- round(100*inputBed[["percent"]], 2)
 
inputBed[["percentAllReads"]] <- inputBed[,"nReads"] / sum(inputBed[,"nReads"])
inputBed[["percentAllReads"]] <- round(100*inputBed[["percentAllReads"]], 2)

write.xlsx(inputBed, file="bed2peak_output_${d}/${out_prefix}.peak.merge.xlsx")

#we don't do spike
#inputBed <- subset(inputBed, nOverlapWithSpike == 0)
#inputBed[["percent"]] <- inputBed[,"nUniqueReads"] / sum(inputBed[,"nUniqueReads"])
#inputBed[["percent"]] <- round(100*inputBed[["percent"]], 2)

#inputBed[["percentAllReads"]] <- inputBed[,"nReads"] / sum(inputBed[,"nReads"])
#inputBed[["percentAllReads"]] <- round(100*inputBed[["percentAllReads"]], 2)

#write.xlsx(inputBed, file="bed2peak_filtered_${d}/${out_prefix}.peak.merge.xlsx")
#write.xlsx(head(inputBed, n=20), file="bed2peak_filtered_${d}/${out_prefix}.top20.peak.merge.xlsx")
EOF

#rm bed2peak_${d}/${out_prefix}.2.xls 
#rm bed2peak_${d}/${out_prefix}.peak.merge.xls
done
done
