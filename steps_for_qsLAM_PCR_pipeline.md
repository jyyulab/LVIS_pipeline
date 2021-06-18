# Calling vector integration sites from qsLAM PCR assay

Download everything in the folder qsLAM_PCR into the working directory. The pipeline consists of a number shell scripts usable in Linux. Some of the scripts contain module load commands to load environmental variables.  This may need to be changed for a particular system. Create a folder named rawdata in the working directory for the input files. The input files are paired-end reads sequenced from the qsLAM-PCR assay. A optional step using [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for quality control can be performed using [01-fastqc.sh](https://github.com/jyyulab/LVIS_pipeline/blob/master/qsLAM_PCR/01-fastqc.sh).

### Step 1 (Reads preprocessing) 
Primer sequences are trimmed using [cutadapt](https://cutadapt.readthedocs.io/en/stable/). For forward reads this is ATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTC (3' end of LTR) and for the reverse the sequence is GACTGCGTATCAGT (PCR adaptor). A maximum error rate of 0.1 is allowed. This is performed using [02-cutadapt.sh](https://github.com/jyyulab/LVIS_pipeline/blob/master/qsLAM_PCR/02-cutadapt.sh). Reads whose length are less than 30 bp after trimming are further filtered [05-makeNewFastq.sh](https://github.com/jyyulab/LVIS_pipeline/blob/master/qsLAM_PCR/05-makeNewFastq.sh), resulting at new fastq files for mapping. An optional script [06-readLen.sh](https://github.com/jyyulab/LVIS_pipeline/blob/master/qsLAM_PCR/06-readLen.sh) can be used to examine the read length.

### Step 2 (Reads alignment) 
Trimmed reads are aligned to the reference genome using the script [10-bwa.sh](https://github.com/jyyulab/LVIS_pipeline/blob/master/qsLAM_PCR/10-bwa.sh) by [BWA](http://bio-bwa.sourceforge.net/).  

### Step 3 (Post-mapping processing)
The script [11-bam2bed.sh](https://github.com/jyyulab/LVIS_pipeline/blob/master/qsLAM_PCR/11-bam2bed.sh) is used for filtering singleton reads, pairs mapped on different chromosomes, pairs with insertion size >1000bp using bedtools [bamtobed](https://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html). The reads (and genomic coordinates) are stored in bed files, including one for all reads and one for reads without duplicates (same start and end locations). An optional script [12-bed2wig.sh](https://github.com/jyyulab/LVIS_pipeline/blob/master/qsLAM_PCR/12-bed2wig.sh) can be used to convert the bed file to bedgraph and bigwig files.

### Step 4 (Integration sites calling)
The genomic starting coordinate of a read corresponds to the location of a vector integration site. After filtering of duplicated reads, the number of reads sharing the same integration sites are essentially originated from different cells. Counting the number of such reads measures the abundance of a clone. Due to a resolution limit of the assay, two vector integration sites (in the same strand) very close together are considered to be identical. This procedure is carried out by bedtools with an window size parameter d, using the script [13-bed2peak.noFilter.sh](https://github.com/jyyulab/LVIS_pipeline/blob/master/qsLAM_PCR/13-bed2peak.noFilter.sh). After arriving at a set of integration sites, apart from the number of unique reads, the total number of reads (including duplicates) starting from a site is also counted. The two counts nUniqueReads and nReads serve as two measures of clonal abundance.  

### Step 5 (Annotate the vector integration sites)
The script [14-peakAnnotate.sh](https://github.com/jyyulab/LVIS_pipeline/blob/master/qsLAM_PCR/14-peakAnnotate.sh) is used to annotate the integration sites like the nearest genes, distance to TSS. Two excel files are generated in the output folder bed2peak_output_d as the final output for each set of FASTQ sequences.  One is the complete list of results, while the other is the top 20. An example of output file is shown below.

  
![image](https://user-images.githubusercontent.com/20668533/117071120-43894780-acf4-11eb-9b22-bc7000587fba.png)

