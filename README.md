# Identification of lentiviral integration sites and analysis
<img src="./img/track.png" width="600" height="600">


This repository includes a few analysis pipelines for lentiviral integrome analysis.
  <ol>
    <li><a href="#sequencing-pipeline-for-our-qsLAM-PCR-assay">Sequencing pipeline for our qsLAM PCR assay</a></li>
    <li><a href="#steps-for-profiling-integration-sites-from-scATAC-seq-and-scMultiome-data">Steps for profiling integration sites from scATAC-seq and scMultiome data</a></li>
    <li><a href="#Downstream-analysis-of-vector-integration-sites">Downstream analysis of vector integration sites</a></li>
    <li><a href="#A-classifier-of-integration-sites-by-integrome-signatures">A classifier of integration sites by integrome signatures</a></li>
  </ol>



## Sequencing pipeline for our qsLAM PCR assay

## Steps for profiling integration sites from scATAC-seq and scMultiome data

### Prerequisites
  * [samtools]
  * [bwa]
  * [R]

See the steps-by-steps [instructions](https://github.com/jyyulab/LVIS_pipeline/blob/master/steps_profile_VIS_scMultiome.md) for identifying integration sites from ATAC-seq data.


## Downstream analysis of vector integration sites

### Prerequisites
  * [R >3.6.2]
  * [bedr 1.0.7]
  * [bedtools 2.29.0]

### Usage
The function pipeline.R shows step-by-step the compilation of VISs, and the functions in LVIS_functions.R. can further be used for determining hotspots and quantifying clonal diversity.

```R
$ source("LVIS_functions.R")
```

### Inputs
The main input is a tab-separated file generated by an in-house pipeline. The pipeline took the reads from a quantitative shearing linear amplification (qsLAM) PCR method, and mapped the reads with insertion junctions the reference genome hg19. The mapped reads correspond to potential VISs specified by the chromosome coordinates. Two potential integration sites that differed by 8 base pairs would be merged. 
An example of input file:
![picture2](./img/sample_input.png)

### Outputs
The pipeline.R script shows the the compilation of a final list of VISs. Based on the compiled list, a frequency matrix storing the total number of reads mapped to each of the integration sites across samples from different cell types and time points was constructed for quantifying clonal populations.

To reduce false-positive results, only integration sites with at least two unique reads or at least five total reads were kept. The total number of reads mapped to a VIS was used to quantify the relative population of the corresponding clone. To compile a complete list of integration sites for a patient, integration sites across samples were matched. In short, overlapping integration sites (defined from the start position to the end position in the VIS calling output) from different samples were considered to be the same site, even though the procedure might potentially merge a few integration sites that were close to one another. 

## A classifier of integration sites by integrome signatures



