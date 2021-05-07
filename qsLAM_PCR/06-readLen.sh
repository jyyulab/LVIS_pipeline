#!/bin/bash
#
#BSUB -P insertionSite         # project code
#BSUB -J readLen     # job name
#BSUB -W 40:00                # wall-clock time (hrs:mins)
#BSUB -R "rusage[mem=8000]"     # memory to reserve, in MB
#BSUB -e errors.%J.hybrid     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J.hybrid     # output file name in which %J is replaced by the job ID
module load R/3.4.0

for i in `ls newFastq | grep -P "fastq$"`
do
	cat newFastq/$i | paste - - - - | cut -f 2 | perl -e 'while(<STDIN>){chomp;print length($_)."\n"; }' | sort | uniq -c | sort -nk 2,2 | sed 's/^ *//' > newFastq/$i.readLen
done

cd newFastq
R --slave <<EOF
options(stringsAsFactors=FALSE)
library("ggplot2")
pdf("readLen.pdf")
files <- dir("./", ".readLen")
for(f in files){
	temp <- read.table(f, sep=" ", header=FALSE)
	for(i in (nrow(temp)-1):1){
		temp[i,1] <- temp[i + 1, 1 ] + temp[i,1]
	}
	temp[,1] <- temp[,1] / temp[1,1]
	temp[,2] <- as.numeric(temp[,2])
	print(ggplot(temp, aes(V2, V1)) + geom_line(colour="blue") + labs(title=f) + xlim(1, 151))
}
dev.off()
EOF
