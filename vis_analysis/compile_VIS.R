##This pipeline starts from the outputs of the qsLAM pipeline.
##Outputs of different samples are concatenated vertically, and an extra column called sample_ori is added. 
##The format of the column sample.ori is
##patient_celltype_time, e.g. P1_CD3_15mo
##Ususally, patients are analyzed separately.

load("~/Desktop/LVXSCID_manuscript/parsed_VIS/all_P1_collected_IS_rename.RData");
dim(all_P1_collected_IS_rename)
setwd("~/Github/LVIS_pipeline/vis_analysis/")

source("~/Github/LVIS_pipeline/vis_analysis/LVIS_functions.R");
#load("hg19.basic.annotation.update.Rdata");

###In some of the early experiments, spikeIn were used. We remove all sites overlap with the spike.

iz1<-which(all_P1_collected_IS_rename$nOverlapWithSpike==0)
all_P1_collected_IS_rename<-all_P1_collected_IS_rename[iz1,];

all_expts<-unique(all_P1_collected_IS_rename$sample_ori);

###To analyze the sites in different samples of the same patient , we need to match sites across samples. 
###In short, two sites from two samples which overlap with each other are identifed as the same site (strand specific).
###Bedtools and the R package bedr are required.

res<-get_merged_ID_from_all_VIS_collection(all_P1_collected_IS_rename);

u_VIS_merge<-res[['u_VIS_merge']];
dim(u_VIS_merge)
##The data frame u_VIS_merge is essentially the dictionary of all the unique insertion sites compiled across all samples.
##Each site is specified by its genome coordinates, and the last column will be used as the index of the site.

X<-res[['X']];
##The dataframe X matches will be the input collection of all insertion sites. 
##Using the last column "u_merge_id", one could tell the unique insertion site index of any site called in a sample.

###The list of insertion sites are re-annotated using HOMER. First, an input file is prepared.
u_VIS_merge_out<-cbind(PeakID=rownames(u_VIS_merge),u_VIS_merge[,1:3],Strand=u_VIS_merge[,6]);
colnames(u_VIS_merge_out)<-c("PeakID","Chr","Start","End","Strand");
write.table(u_VIS_merge_out,file="all_P1_VIS_merge.txt",sep="\t",quote=FALSE,row.names = FALSE);

###Then HOMER is run separately. The annotation is stored in the out.txt file.
VIS_list="all_P1_VIS_merge";
module load homer/4.9.1
annotatePeaks.pl $VIS_list.txt hg19 -annStats $VIS_list.homer.stat.txt  > $VIS_list.homer.out.txt

library(openxlsx);
tmp<-read.xlsx('all_P1_VIS_merge.homer.out.xlsx');
rownames(tmp)<-tmp[,1];
u_VIS_merge_homer_annot<-tmp[rownames(u_VIS_merge),];

####Finally, we build two count matrices, one for the number of reads mapped to the sites in all samples
####one for the number of unique shear sites reads mapped to the sites in all samples
iu<-which(!is.na(X$u_merge_ind))
X<-X[iu,];
counts_nR_IS_expts<-generate_u_VIS_merge_samples_matrix(X);
dim(counts_nR_IS_expts)#10270    58

counts_nSS_IS_expts<-generate_u_VIS_merge_samples_matrix_unique(X);
dim(counts_nSS_IS_expts)#

all_P1_collected_IS_clean<-all_P1_collected_IS_rename;

save(all_P1_collected_IS_clean,counts_nR_IS_expts,counts_nSS_IS_expts,u_VIS_merge_homer_annot,u_VIS_merge,file="all_P1_source.Rdata");











