###the starting point is all_P1_collected_IS
###based on the sample_ori column, we subset all VIS from a certain expt.
###we standardize sample_ori in the form Px_celltype_time, make sure boost is transformed.
###
library(ggplot2);

count_VIS_freq <- function(tmp){
	U<-tmp$nReads;#we decided to use nReads to count abundance
	if (length(U)<21){
		U<-c(U,rep(0,1,21-length(U)));
	}
	ord<-order(U,decreasing=TRUE);
  	slices<-c(U[ord][1:20],sum(U[ord][21:length(U)]));
  	x<-slices/sum(slices);
  	return(x);
}

#use dataframe of all_collected_VIS as input
get_top20_composition_and_nVIS_all_samples <- function(all_P1_collected_IS_rename){
	print('remember to get rid of the VISs with spike');
	all_expts<-unique(all_P1_collected_IS_rename$sample_ori);
	top20_IS_freq_vs_expts<-matrix(0,21,length(all_expts));
	colnames(top20_IS_freq_vs_expts)<-all_expts;
	rownames(top20_IS_freq_vs_expts)<-c(1:20,'rest');
	nVIS<-matrix(0,length(all_expts),1);
	colnames(nVIS)<-'nVIS';
	rownames(nVIS)<-all_expts;
	for (e in all_expts){
		tmp<-all_P1_collected_IS_rename[which(all_P1_collected_IS_rename$sample_ori==e),];
		nVIS[e,1]<-dim(tmp)[1];
		x<-count_VIS_freq(tmp);
		top20_IS_freq_vs_expts[,e]<-x;
	}
	res<-list();
	res[['num_VIS']]<-nVIS;
	res[['top20_IS_freq']]<-top20_IS_freq_vs_expts;
	return(res);
}

get_samples_diversity_measures <- function(tmp){

	rr<-tmp$nReads;
	y<-rr[rr>0];
	y<-as.numeric(cumsum(sort(y)));
	y<-c(0,y/y[length(y)]);
	x<-seq(0,1,1/(length(y)-1));
	library(pracma);
	sample_OCI<-1-2*trapz(x,y);
	y<-rr[rr>0];
	y<-as.numeric(y);
	y<-y/sum(y);
	sample_entropy<-sum(y*log2(y))*-1;
	y<-rr[rr>0];
	D<-length(y);
	ny<-table(y);
	f1<-sum(ny==1);
	f2<-sum(ny==2);
	sample_Chao<-D+f1*(f1-1)/(2*(f2+1));
	y<-rr[rr>0];
	y<-as.numeric(sort(y,decreasing = TRUE));
	p<-y/sum(y);
	cp<-cumsum(p);
	sample_UC50<-which(cp>0.5)[1];
	res<-list();
	res[['sample_OCI']]<-sample_OCI;
	res[['sample_entropy']]<-sample_entropy;
	res[['sample_Chao']]<-sample_Chao;
	res[['sample_UC50']]<-sample_UC50;
	res[['sample_nVIS']]<-length(rr);

	return(res);
}

get_diversity_measures_all_samples <- function(all_P1_collected_IS_rename){
	print('remember to get rid of the VISs with spike');
	all_expts<-unique(all_P1_collected_IS_rename$sample_ori);
	diversity_vs_expts<-matrix(0,4,length(all_expts));
	colnames(diversity_vs_expts)<-all_expts;
	rownames(diversity_vs_expts)<-c('sample_OCI','sample_entropy','sample_Chao','sample_UC50');
	for (e in all_expts){
		tmp<-all_P1_collected_IS_rename[which(all_P1_collected_IS_rename$sample_ori==e),];
		out<-get_samples_diversity_measures(tmp);
		diversity_vs_expts[1,e]<-out[['sample_OCI']];
		diversity_vs_expts[2,e]<-out[['sample_entropy']];
		diversity_vs_expts[3,e]<-out[['sample_Chao']];
		diversity_vs_expts[4,e]<-out[['sample_UC50']];
	}
	return(diversity_vs_expts);
}

prepare_pie_charts <- function(top20_IS_freq,num_VIS){
	library(reshape2);
	data<-melt(top20_IS_freq);
	colnames(data)<-c("VIS_label","sample_ori","composition");
	df = transform(data, expt=colsplit(sample_ori, "_", names = c('patients','Celltype','time')));
	#df$expt.Celltype<-as.character(df$expt$Celltype);
	#colnames(df)<-c("VIS","expt","composition","patients","time","Celltype");
	iz<-which(df$expt.Celltype=='bulk');
	df$expt.Celltype[iz]<-'TNC';
	#df$expt.Celltype<-factor(df$expt.Celltype,levels=c('CD3','CD14CD15','CD19','CD56','TNC'));
	#df<-cbind(df,patient_time=paste0(df$patients,":",df$time))
	#df$patient_time<-factor(df$patient_time,levels=c('P1:30mo','P2:24mo','P3:24mo','P4:18mo','P5:18mo','P6:18mo','P7:12mo','P8:12mo','P10:6mo','P11:16wks'));
	all_expts<-colnames(top20_IS_freq);
	dat_text <- data.frame(
  	sample_ori=all_expts,
  	nVIS = num_VIS#
	)
	dt = transform(dat_text, expt=colsplit(sample_ori, "_", names = c('patients','Celltype','time')));
	#dt$expt.Celltype<-as.character(dt$expt.Celltype);
	#colnames(dt)<-c("sample_ori","nVIS",'expt.patients', 'expt.time','expt.Celltype');
	iz<-which(dt$expt.Celltype=='bulk');
	dt$expt.Celltype[iz]<-'TNC';
	#dt$expt.Celltype<-factor(dt$expt.Celltype,levels=c('CD3','CD14CD15','CD19','CD56','TNC'));

	res<-list();
	res[['df']]<-df;
	res[['dt']]<-dt;
	print('modify the columns of df, and generate the required factors before using the function draw_pie_charts()');
	return(res);

}

draw_pie_charts <- function(df,dt){

	library(ggplot2);
	library(scales);
	color_p<-hue_pal()(21);
	color_p[21]<-"grey";

	p<- ggplot(df, aes(x="", y=composition, color=VIS_label, fill=VIS_label)) + 
	geom_bar(width = 1, stat = "identity") +
	#coord_polar("y", start=0) + facet_grid(expt.time ~ expt.Celltype ) +
	coord_polar("y", start=0) + facet_grid(expt.Celltype ~ expt.time ) +
	theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12)) + scale_fill_manual(values=color_p) + scale_color_manual(values=color_p);

	p <- p + geom_text(aes(x=0, y=0,label = as.character(nVIS)),data=dt,inherit.aes=FALSE, parse=FALSE,size=5);

	return(p);
	#ggsave("all_patients_last_samples_composition.pdf",plot = last_plot(), device = NULL, path = NULL,scale = 1, width = 10, height = 14, dpi = 300);
}


preprocess_all_VIS_collection<- function(all_P1_collected_IS){
	#all_P1_collected_IS is the collection of output from the same patient
	all_IS_P1<-paste0(all_P1_collected_IS$seqnames,'.',all_P1_collected_IS$start,'.',all_P1_collected_IS$end,'.',all_P1_collected_IS$strand);
	all_IS_P1_aux<-paste0(all_P1_collected_IS$seqnames,':',all_P1_collected_IS$start,'-',all_P1_collected_IS$end);
	
	#20152
	u_all_IS_P1<-unique(all_IS_P1);
	u_all_IS_P1<-read.table(text = u_all_IS_P1, sep = ".");
	library(bedr);
	x<-paste0(u_all_IS_P1[,1],":",u_all_IS_P1[,2],'-',u_all_IS_P1[,3]);
	is.x.valid  <- is.valid.region(x);
	#get rid of other chr
	x <- x[is.x.valid];
	u_all_IS_P1<-u_all_IS_P1[is.x.valid,];
	dim(u_all_IS_P1)

	u_all_IS_P1<-cbind(u_all_IS_P1[,c(1:3)],tmp1='.',tmp2='.',strand=u_all_IS_P1[,4]);
	colnames(u_all_IS_P1)<-c("chr","st","ed","tmp1","tmp2","strand");
	rownames(u_all_IS_P1)<-paste0(u_all_IS_P1$chr,':',u_all_IS_P1$st,'-',u_all_IS_P1$ed);

	x<-rownames(u_all_IS_P1);
	x.sort <- bedr.sort.region(x);
	u_all_IS_P1<-u_all_IS_P1[x.sort,];
	u_all_IS_P1[,1]<-as.character(u_all_IS_P1[,1]);
	
	map_back<-match(all_IS_P1,paste0(u_all_IS_P1[,1],'.',u_all_IS_P1[,2],'.',u_all_IS_P1[,3],'.',u_all_IS_P1[,6]));
	#max(map_back[!is.na(map_back)]);
	#how original collected IS go to the u_all_IS_P1

	u_all_IS_P1<-cbind(u_all_IS_P1[,c(1:3)],ind=c(1:dim(u_all_IS_P1)[1]),aux=1,strand=u_all_IS_P1[,6]);
	

	all_P1_collected_IS_mapped<-u_all_IS_P1[map_back,];
	all_P1_collected_IS_mapped<-cbind(all_P1_collected_IS,u_ind=all_P1_collected_IS_mapped$ind);

	res<-list();
	res[['u_all_IS_sort']]<-u_all_IS_P1;
	res[['all_collected_IS_mapped']]<-all_P1_collected_IS_mapped;
	#back to the input, use the ind col to back to the unique 
	return(res);

}

merge_all_preprocessed_VIS_collection <-function(u_all_IS_P1){
	###we can merge adjacinet VIS..with d=8? d=4? d=0?
	u_all_IS_P1.merge<-bedr(
		engine = "bedtools", 
     	input = list(i = u_all_IS_P1), 
        method = "merge", 
        params = "-s -d 0 -c 5,4,6 -o sum,collapse,distinct"
        );

	u_all_IS_P1.merge<-cbind(u_all_IS_P1.merge,eff_ind=c(1:dim(u_all_IS_P1.merge)[1]));

	U<-u_all_IS_P1.merge$aux;
	V<-u_all_IS_P1.merge$eff_ind;
	U_final<-c();
	V_final<-c();
	for (i in 1:length(V)){
		tmp<-as.numeric(strsplit(U[i],",")[[1]]);
		n<-length(tmp);
		U_final<-c(U_final,tmp);
		V_final<-c(V_final,rep(V[i],n));
	}
	map<-cbind(ori_index=U_final,new_index=V_final);
	rownames(map)<-U_final;
	res<-list();
	res[['u_VIS_merge']]<-u_all_IS_P1.merge;
	res[['map']]<-data.frame(map);

	return(res);
}

get_merged_ID_from_all_VIS_collection <- function(all_P1_collected_IS){
	
	res1<-preprocess_all_VIS_collection(all_P1_collected_IS);
	u_all_IS_P1<-res1[['u_all_IS_sort']];
	all_P1_collected_IS_mapped<-res1[['all_collected_IS_mapped']];

	#x<-all_IS_P1_collected_mapped$u_ind;
	#x<-x[!is.na(x)];
	#max(x)

	res2<-merge_all_preprocessed_VIS_collection(u_all_IS_P1);

	u_VIS_merge<-res2$u_VIS_merge;
	u_VIS_to_merged_VIS<-res2$map;

	#head(all_P1_collected_IS)
	#head(all_IS_P1_collected_mapped)
	#iz<-head(all_IS_P1_collected_mapped)$u_ind;

	#the unique VIS index to the merged index..
	#u_VIS_to_merged_VIS[iz,]

	X<-cbind(all_P1_collected_IS_mapped,u_merge_ind=u_VIS_to_merged_VIS[all_P1_collected_IS_mapped$u_ind,]$new_index);
	#X[order(X$nReads,decreasing = TRUE),];

	res<-list();
	res[['X']]<-X;
	res[['u_VIS_merge']]<-u_VIS_merge;

	return(res);
}

generate_u_VIS_merge_samples_matrix <- function(X){
	all_expts<-unique(X$sample_ori);
	n<-max(X$u_merge_ind)
	uID<-matrix(0,n,length(all_expts));
	rownames(uID)<-c(1:n);
	colnames(uID)<-all_expts;

	for (e in all_expts){
		iz<-which(X$sample_ori==e);
		X2<-X[iz,];
		utmp<-unique(X2$u_merge_ind);
		vtmp<-c();
		for (uu in utmp){
			iuu<-which(X2$u_merge_ind==uu);
			vtmp<-c(vtmp,sum(X2$nReads[iuu]));
		}
		uID[utmp,e]<-vtmp;
	}
	return(uID);

}

prepare_venn_diagram_samples <- function(counts_IS_expts,pick_samples){
	aux<-list();
	for (i in 1:length(pick_samples)){
		aux[[pick_samples[i]]]<-names(which(counts_IS_expts[,pick_samples[i]]>0));
	}
	return(aux);
}






get_genomicDensity<-function(u_VIS_merge,win_size){
	
	tmp<-u_VIS_merge[,c(1,2,3,6)];
	colnames(tmp)<-c('chr','start','end','strand');
	d<-(tmp$end-tmp$start)/2;
	tmp$start<-floor(tmp$start+d);
	tmp$end<-floor(tmp$start+d)+1;
	VIS_list_clean_overall<-tmp;
	
	VIS_list_clean_pos<-VIS_list_clean_overall[which(VIS_list_clean_overall$strand=='+'),];
	VIS_list_clean_neg<-VIS_list_clean_overall[which(VIS_list_clean_overall$strand=='-'),];

	library(circlize);
	VIS1_overall<-genomicDensity(VIS_list_clean_overall, window.size = win_size);
	VIS1_pos<-genomicDensity(VIS_list_clean_pos, window.size = win_size);
	VIS1_neg<-genomicDensity(VIS_list_clean_neg, window.size = win_size);

	VIS_density<-list();
	VIS_density[['overall']]<-VIS1_overall;
	VIS_density[['pos']]<-VIS1_pos;
	VIS_density[['neg']]<-VIS1_neg;

	return(VIS_density);
}

#input a few bed files as a list. for each file, we have chr, start, end storing the VIS
#row names not matter..
draw_circos<-function(VIS_list){
	library(circlize);
	gap<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);gap[24]<-5;
	circos.par(gap.after=gap);
	circos.initializeWithIdeogram();
	n<-length(VIS_list);
	for (i in 1:n){
		circos.genomicDensity(VIS_list[[i]], col = c("#0000FF80"), track.height = 0.1,window.size=1e6,ylim.force = TRUE);
	}
	circos.clear();
}

draw_composition_stackbars<-function(u_VIS_merge_homer_annot,counts_IS_expts,pick_samples,pick_VIS){

	if (missing(pick_VIS)){
		aux<-prepare_venn_diagram_samples(counts_IS_expts,pick_samples);
		common_VIS<-aux[[1]];
		for (i in 1:length(aux)){
  			common_VIS<-intersect(common_VIS,aux[[i]]);
		}
		pick_VIS<-common_VIS;
	}
	
	X<-counts_IS_expts[,pick_samples];
	Z<-sweep(X,2,colSums(X),`/`);
	write.table(Z[pick_VIS,],file="tmp.txt");

	data<-read.table('tmp.txt',header=TRUE);
	rest<-t(data.frame(1-colSums(data)));
	rownames(rest)<-"rest";
	tmp<-u_VIS_merge_homer_annot[as.numeric(rownames(data)),];
	rownames(data)<-paste0(rownames(tmp),':',tmp$Gene.Name);
	data<-rbind(data,rest);

	df<-melt(cbind(VIS=rownames(data),data));
	df$VIS<-factor(df$VIS,levels=rownames(data));

	library(scales)
	n<-dim(data)[1];
	color_p<-hue_pal()(n);
	color_p[n]<-"grey";

	p <- ggplot(df, aes(x = variable, y = value, fill = VIS)) + 
  	geom_bar(stat = "identity") +
  	ylab("Fraction") +
  	theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.text=element_text(size=12)) +
  	scale_fill_manual(values=color_p) + scale_color_manual(values=color_p);

  	res<-list();
  	res[['ggplot_obj']]<-p;
  	res[['pick_VIS_composition']]<-data;
  	return(res);

}


# get_VIS_chr<-function(VIS1,chr){
# 	iz<-which(VIS1$chr==chr);
# 	vec<-VIS1[iz,]$pct;
# 	return(vec);
# }

# get_hotspots_thres<-function(P2_hot,nVIS){
# 	window=1e5;
# 	genome_size=3e9;
# 	n<-genome_size/window;
# 	lambda<-nVIS/n;
# 	expect_count=P2_hot$pct*window;
# 	all_Pvalue<-1-ppois(expect_count,lambda);
# 	P2_hot$P<-all_Pvalue;
# 	return(P2_hot);
# }



