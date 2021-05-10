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

### Dependencies
  * [fastqc >0.11.5]
  * [cutadapt]
  * [samtools >1.10]
  * [bwa >0.7.17]
  * [bedtools >2.25.0]
  * [R >3.6.2]

### Usage
Download everything inside the folder qsLAM into the working directory. Create a folder called rawdata and put the paired-end reads files inside. Follow the steps-by-steps [instructions](https://github.com/jyyulab/LVIS_pipeline/blob/master/steps_for_qsLAM_PCR_pipeline.md).

## Steps for profiling integration sites from scATAC-seq and scMultiome data

### Prerequisites
  * [fastqc >0.11.5]
  * [samtools >1.10]
  * [bwa >0.7.17]

### Usage
See the steps-by-steps [instructions](https://github.com/jyyulab/LVIS_pipeline/blob/master/steps_profile_VIS_scMultiome.md) for identifying integration sites from ATAC-seq data.


## Downstream analysis of vector integration sites

### Prerequisites
  * [R >3.6.2]
  * [bedr 1.0.7]
  * [bedtools 2.29.0]

### Usage
Useful functions are implemented in LVIS_functions.R. See examples and steps-by-steps [instructions](https://github.com/jyyulab/LVIS_pipeline/edit/master/more_vis_analysis.md) for the compilation of VISs across samples and other downstream analysis.

```R
$ source("LVIS_functions.R")
```

## A classifier of integration sites by integrome signatures



