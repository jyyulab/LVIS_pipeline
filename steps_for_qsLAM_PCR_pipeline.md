# Calling vector integration sites from qsLAM PCR assay

Download everything in the folder qsLAM_PCR into the working directory. The pipeline consists of a number shell scripts usable in Linux. Some of the scripts contain module load commands to load environmental variables.  This may need to be changed for a particular system. Create a folder named rawdata in the working directory for the input files. The input files are paired-end reads sequenced from the qsLAM-PCR assay. A optional step using [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for quality control can be performed using [01-fastqc.sh](https://github.com/jyyulab/LVIS_pipeline/blob/master/qsLAM_PCR/01-fastqc.sh).

### Step 1 (Reads preprocessing) 
Primer sequences are trimmed using [cutadapt](https://cutadapt.readthedocs.io/en/stable/). For forward reads this is ATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTC (3' end of LTR) and for the reverse the sequence is GACTGCGTATCAGT (PCR adaptor). A maximum error rate of 0.1 is allowed. This is performed using [02-cutadapt.sh](https://github.com/jyyulab/LVIS_pipeline/blob/master/qsLAM_PCR/02-cutadapt.sh). Reads whose length are less than 30 bp after trimming are further filtered [05-makeNewFastq.sh](https://github.com/jyyulab/LVIS_pipeline/blob/master/qsLAM_PCR/05-makeNewFastq.sh), resulting at new fastq files for mapping. An optional script [06-readLen.sh](https://github.com/jyyulab/LVIS_pipeline/blob/master/qsLAM_PCR/06-readLen.sh) can be used to examine the read length.

### Step 2 (Reads alignment) 
Trimmed reads are aligned to the reference genome using the script [10-bwa.sh](https://github.com/jyyulab/LVIS_pipeline/blob/master/qsLAM_PCR/10-bwa.sh) by [BWA](http://bio-bwa.sourceforge.net/).  

### Step 3 (Post-mapping processing)
The script [11-bam2bed.sh](https://github.com/jyyulab/LVIS_pipeline/blob/master/qsLAM_PCR/11-bam2bed.sh) is used for filtering singleton reads, pairs mapped on different chromosomes, pairs with insertion size >1000bp using bedtools [bamtobed](https://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html). The reads (and genomic coordinates) are stored in bed files, including one for only non-duplicated reads. An optional script [12-bed2wig.sh](https://github.com/jyyulab/LVIS_pipeline/blob/master/qsLAM_PCR/12-bed2wig.sh) can be used to convert the bed file to bedgraph and bigwig files.

### Step 4 (Integration sites calling)
13-bed2peak.noFilter.sh
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
### Step 5 (Annotate the vector integration sites)
The script [14-peakAnnotate.sh](https://github.com/jyyulab/LVIS_pipeline/blob/master/qsLAM_PCR/14-peakAnnotate.sh) is used to annotate the integration sites like the nearest genes, distance to TSS. Two excel files are generated in the output folder bed2peak_output_d as the final output for each set of FASTQ sequences.  One is the complete list of results, while the other is the top 20.

  
![image](https://user-images.githubusercontent.com/20668533/117071120-43894780-acf4-11eb-9b22-bc7000587fba.png)

