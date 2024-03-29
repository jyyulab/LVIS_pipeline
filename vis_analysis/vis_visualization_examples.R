source("LVIS_functions.R");
load("hg19.basic.annotation.update.Rdata");

load("all_P4_collected_IS.RData");
dim(all_P4_collected_IS)

###To make pie charts of clone diversity. 
res<-get_top20_composition_and_nVIS_all_samples(all_P4_collected_IS);
num_VIS<-res[['num_VIS']];
top20_IS_freq_vs_expts<-res[['top20_IS_freq']];

out<-prepare_pie_charts(top20_IS_freq_vs_expts,num_VIS);
df<-out$df;
dt<-out$dt;

##focus on sites from P1 only
df1<-df[which(df$expt.patients=='P4'),];
dt1<-dt[which(dt$expt.patients=='P4'),];

df1$expt.Celltype<-factor(df1$expt.Celltype,levels=c('CD3','CD19','CD56','CD14CD15','TNC'));
df1$expt.time<-factor(df1$expt.time,levels=c('12wks','16wks','6mo','9mo','12mo','15mo','16mo','18mo','21mo','24mo','27mo','30mo'));

dt1$expt.Celltype<-factor(dt$expt.Celltype,levels=c('CD3','CD19','CD56','CD14CD15','TNC'));
dt1$expt.time<-factor(dt$expt.time,levels=c('12wks','16wks','6mo','9mo','12mo','18mo','24mo'));

p1<-draw_pie_charts(df1,dt1);p1;
ggsave("all_P4_samples_VIS_composition.pdf",plot = last_plot(), device = NULL, path = NULL,scale = 1, width = 10, height = 14, dpi = 300);

###We implemented a few metrics to quantify sample diveristy, including entropy, Chao estimator, OCI, UC50.
samples_diversity<-get_diversity_measures_all_samples(all_P4_collected_IS);


####To generate GenesCloud
load(paste0(fcts_loc,"hg19.basic.annotation.update.Rdata"));
u_VIS_merge_homer_annot_modify<-simplify_homer_annotation(u_VIS_merge_homer_annot,M);

geneFreq<-get_geneFreq(u_VIS_merge_homer_annot_modify);
generate_GeneCloud(geneFreq,'P4_geneCloud.pdf');


