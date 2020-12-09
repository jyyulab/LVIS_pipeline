##This pipeline starts from the outputs of the VIS calling pipeline.
##Outputs of different samples are concatenated vertically, and an extra column called sample_ori is added. 
##The format of the column sample.ori is
##patient_celltype_time, e.g. P1_CD3_15mo
##Ususally, patients are analyzed separately.

source("LVIS_functions.R");
#load("hg19.basic.annotation.update.Rdata");

load("/Volumes/yu3grp/IO_JY/yu3grp/LVXSCID/VIS_20190708/parsed_VIS/all_P1_collected_IS_rename.RData");
dim(all_P1_collected_IS_rename)

###In some of the early experiments, spikeIn were used. We remove all sites overlap with the spike.

iz1<-which(all_P1_collected_IS_rename$nOverlapWithSpike==0)
all_P1_collected_IS_rename<-all_P1_collected_IS_rename[iz1,];

all_expts<-unique(all_P1_collected_IS_rename$sample_ori);

###To make pie charts of clone diversity. 
res<-get_top20_composition_and_nVIS_all_samples(all_P1_collected_IS_rename);
num_VIS<-res[['num_VIS']];
top20_IS_freq_vs_expts<-res[['top20_IS_freq']];

out<-prepare_pie_charts(top20_IS_freq_vs_expts,num_VIS);
df<-out$df;
dt<-out$dt;

##focus on sites from P1 only
df1<-df[which(df$expt.patients=='P1'),];
dt1<-dt[which(dt$expt.patients=='P1'),];

df1$expt.Celltype<-factor(df1$expt.Celltype,levels=c('CD3','CD19','CD56','CD14CD15','TNC'));
df1$expt.time<-factor(df1$expt.time,levels=c('12wks','16wks','6mo','9mo','12mo','15mo','16mo','18mo','21mo','24mo','27mo','30mo'));

dt1$expt.Celltype<-factor(dt$expt.Celltype,levels=c('CD3','CD19','CD56','CD14CD15','TNC'));
dt1$expt.time<-factor(dt$expt.time,levels=c('12wks','16wks','6mo','9mo','12mo','18mo','24mo'));


p1<-draw_pie_charts(df1,dt1);p1;
ggsave("all_P1_samples_VIS_composition.pdf",plot = last_plot(), device = NULL, path = NULL,scale = 1, width = 10, height = 14, dpi = 300);

###We implemented a few metrics to quantify sample diveristy, including entropy, Chao estimator, OCI, UC50.
samples_diversity<-get_diversity_measures_all_samples(all_P1_collected_IS_rename);

###To analyze the sites in different samples of the same patient , we need to match sites across samples. In short, two sites from two samples which overlap with each other are identifed as the same site.
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

####Finally, we build a count matrix for each of the unique insertion sites across different experiments.
iu<-which(!is.na(X$u_merge_ind))
X<-X[iu,];
counts_IS_expts<-generate_u_VIS_merge_samples_matrix(X);
dim(counts_IS_expts)#10270    58

all_P1_collected_IS_clean<-all_P1_collected_IS_rename;

save(all_P1_collected_IS_clean,counts_IS_expts,u_VIS_merge_homer_annot,u_VIS_merge,file="all_P1_source.Rdata");











