# Calling vector integration sites from qsLAM PCR assay

Download everything in the folder qsLAM_PCR into the working directory. The pipeline consists of a number shell scripts usable in Linux. Some of the scripts contain module load commands to load environmental variables.  This may need to be changed for a particular system. Create a folder named rawdata in the working directory for the input files. The input files are paired-end reads sequenced from the qsLAM-PCR assay. A optional step using [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for quality control can be performed using [01-fastqc.sh](https://github.com/jyyulab/LVIS_pipeline/blob/master/qsLAM_PCR/01-fastqc.sh).

Step 1 (Reads preprocessing) Primer sequences are trimmed using [cutadapt](https://cutadapt.readthedocs.io/en/stable/). For forward reads this is ATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTC and for the reverse the sequence is GACTGCGTATCAGT. A maximum error rate of 0.1 is allowed. This is performed using [02-cutadapt.sh](https://github.com/jyyulab/LVIS_pipeline/blob/master/qsLAM_PCR/02-cutadapt.sh). Reads whose length are less than 30 bp after trimming are further filtered [05-makeNewFastq.sh](https://github.com/jyyulab/LVIS_pipeline/blob/master/qsLAM_PCR/05-makeNewFastq.sh), resulting at new fastq files for mapping. 


### 10.bwa.sh

Bwa (http://bio-bwa.sourceforge.net/) is used to align reads in the newFastq directory to the hg19 human genome.  Results are placed in the bwa directory.

### 11.bam2bed.sh

Samtools (http://www.htslib.org/) is used to view and sort bam files produced in the bwa directory into sam files.  (Soft and hard clipped reads are...?)  .  The bedtools bamtobed (https://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html) program is used with the  -bedpe and -mate1 options to convert the  sam files into bed files.
(????? more steps I dont understand)

    samtools view -H bwa/$i > bam2bed/${out_prefix}.sam
	#samtools view  -F 4 -q 1 bwa/$i | perl -ne '$line=$_; @rec=split("\t",$line); $rec[5]=~ s/[SH].*//; if( !($rec[5] =~ /^[0-9]+$/ && $rec[5] > 6)){print $line;}' >> bam2bed/${out_prefix}.sam
	samtools view  -F 4  bwa/$i | perl -ne '$line=$_; @rec=split("\t",$line); $rec[5]=~ s/[SH].*//; if( !($rec[5] =~ /^[0-9]+$/ && $rec[5] > 6)){print $line;}' >> bam2bed/${out_prefix}.sam
	samtools sort -n <(samtools view -hbS bam2bed/${out_prefix}.sam ) bam2bed/${out_prefix}
	bedtools bamtobed -i bam2bed/${out_prefix}.bam -bedpe -mate1 > bam2bed/${out_prefix}.temp
	cat bam2bed/${out_prefix}.temp | perl -e 'while(<STDIN>){$line=$_; @rec=split("\t", $line); if($rec[0] eq $rec[3]){ $start = ($rec[1], $rec[4])[$rec[1] > $rec[4]]; $end = ($rec[2], $rec[5])[$rec[2] < $rec[5]]; print "$rec[0]\t$start\t$end\t$rec[6]\t$rec[7]\t$rec[8]\n";}}' | awk -F "\t" '{if(($3-$2)<1000){print}}' > bam2bed/${out_prefix}.bed
	cut -f 1,2,3,6 bam2bed/${out_prefix}.bed | sort -u | awk '{print $1"\t"$2"\t"$3"\t.\t.\t"$4}' > bam2bed/${out_prefix}.rmdup.bed
	rm bam2bed/${out_prefix}.temp bam2bed/${out_prefix}.bam


Final bed files are placed in the bam2bed directory.

### 12-bed2wig.sh

The bed files from the prvious step are converted to bigwig (bw) files in this step.  The bedtools genomecov (https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html) program is used to create a bedgraph file.  The -bg (output in bedgraph fortmat) and -strand (calculate coverage from specific strands)  options are used.  
bedGraphToBigWig (https://www.encodeproject.org/software/bedgraphtobigwig/) is then used to convert the bedgraph files into bigwig files.



### 13-bed2peak.noFilter.sh

The final step uses the bed files to find the location of viral integration site and generates a  report.  An window size parameter d (used as $d) is a window size to merge ends of reads.



    awk -F "\t" '{if($6=="+"){print $1"\t"$2"\t"$6} else {print $1"\t"$3"\t"$6}}' bam2bed/${out_prefix}.bed | 
    grep -v vector |
    sort |
    uniq -c |
    awk '{print $2"\t"$3"\t"$3+1"\t.\t"$1"\t"$4}' |
    sort -k1,1 -k2,2n |
    bedtools merge -c 5 -s -o sum -d $d -i - |
    awk '{print $1"\t"$2"\t"$3"\t.\t"$5"\t"$4}' | 
    grep -v vector |
    bedtools intersect -b <(grep -v vector bam2bed/${out_prefix%.rmdup}.bed|
    awk '{if($6=="+"){print $1"\t"$2"\t"$2+1"\t.\t.\t"$6} else {print $1"\t"$3"\t"$3+1"\t.\t.\t"$6}}') -c -s -a - | \
    awk '{print $1"\t"$2"\t"$3"\t""'$out_prefix'""\t"$5"\t"$6"\t"$7}' > bed2peak_${d}/${out_prefix}.peak.merge.xls

Breaking this up:



The first step of the line, prints the chromosome, and fist position (0 based) of the match if forward strand and end position (1 based) if not on forward strand.  (why 0 or 1 based?)
    awk -F "\t" '{if($6=="+"){print $1"\t"$2"\t"$6} else {print $1"\t"$3"\t"$6}}' bam2bed/${out_prefix}.bed 

Any sequence match vector is removed.   (But this is not in bwa index?)
    grep -v vector



Sequences are sorted and duplicate line are removed (with their counts retained as first column).  
    sort | uniq -c



This is then converted back to a bed file (containing the first position of reads).  

    awk '{print $2"\t"$3"\t"$3+1"\t.\t"$1"\t"$4}'

Sort by chromosome and position
     sort -k1,1 -k2,2n 


Bedtools merge is usede to merge positions within N bases ($d) in line above.  And the number of counts is retained) 
    bedtools merge -c 5 -s -o sum -d $d -i -
    



Revert this back to a bed format.
    awk '{print $1"\t"$2"\t"$3"\t.\t"$5"\t"$4}' 
  

Remove lines with vector (why again)

    grep -v vector


 Intersect these positions with the original bed file??
      bedtools intersect -b <(grep -v vector bam2bed/${out_prefix%.rmdup}.bed|



---
In the next step the bedtools window program is used to identify any control sites.  the window size is identified as $d .  the -sm parameter signals to use reads from same strand.

    bedtools window -w $d -a bed2peak_${d}/${out_prefix}.peak.merge.xls -b known20sites.bed -sm -c > bed2peak_${d}/${out_prefix}.2.xls 

---
A series of R commands used with an R script annotate hits

    R --slave <<EOF
    source("target_gene_prediction.R")
    inputBed <- read.table("bed2peak_${d}/${out_prefix}.2.xls", sep="\t", check.names =FALSE )
    colnames(inputBed) <- c("seqnames", "start", "end", "name", "nUniqueReads",   "strand", "nReads", "nOverlapWithSpike")
    inputBed[["name"]] <- gsub("-", ".", paste("X", inputBed[["name"]], 1:nrow(inputBed), sep="_"))
    d <- target_gene_prediction(inputBed)

    inputBed=inputBed[order(inputBed$name),]
    rownames(d) <- d[["peak_id"]]
    inputBed[["gene"]] <- d[inputBed[["name"]], "nearest_gene_symbol"]
    inputBed[["tss_distances"]] <- d[inputBed[["name"]], "tss_distances"]
     inputBed[["gene_region"]] <- d[inputBed[["name"]], "gene_region"]
    inputBed[["name"]] <- sub("^X_", "", inputBed[["name"]])
    inputBed <- inputBed[order(inputBed[["nUniqueReads"]], decreasing=TRUE),]
    write.table(inputBed, file="bed2peak_${d}/${out_prefix}.peak.merge.xls",     sep="\t", quote=FALSE, row.names=FALSE)
    EOF


---
A set of R commands then process data to an exce

    R --slave <<EOF
    options(java.parameters = "-Xmx8000m")
    library("xlsx")
    options(stringsAsFactors = FALSE)
    inputBed <- read.table("bed2peak_${d}/${out_prefix}.peak.merge.xls", sep="\t",     header=TRUE, check.names =FALSE )
    inputBed <- subset(inputBed, nUniqueReads > 1 | nReads > 5 )
    inputBed[["percent"]] <- inputBed[,"nUniqueReads"] /   sum(inputBed[,"nUniqueReads"])
     inputBed[["percent"]] <- round(100*inputBed[["percent"]], 2)

    inputBed[["percentAllReads"]] <- inputBed[,"nReads"] / sum(inputBed[,"nReads"])
    [["percentAllReads"]] <- round(100*inputBed[["percentAllReads"]], 2)

    write.xlsx(inputBed, file="bed2peak_${d}/${out_prefix}.peak.merge.xlsx")
    #inputBed <- subset(inputBed, nOverlapWithSpike == 0)
    inputBed[["percent"]] <- inputBed[,"nUniqueReads"] / sum(inputBed[,"nUniqueReads"])
    inputBed[["percent"]] <- round(100*inputBed[["percent"]], 2)

    inputBed[["percentAllReads"]] <- inputBed[,"nReads"] / sum(inputBed[,"nReads"])
    inputBed[["percentAllReads"]] <- round(100*inputBed[["percentAllReads"]], 2)

    write.xlsx(inputBed,  file="bed2peak_filtered_${d}/${out_prefix}.peak.merge.xlsx")
    write.xlsx(head(inputBed, n=20),    file="bed2peak_filtered_${d}/${out_prefix}.top20.peak.merge.xlsx")
    EOF

##  Output

Output is placed in the bed2peak_filtered_N directory,  where N is the distance used for combining reads.  Are two excel files for each set of FASTQ sequences.  One is the complete list of results, while the other is the top 20.


![image](https://user-images.githubusercontent.com/20668533/117071120-43894780-acf4-11eb-9b22-bc7000587fba.png)

