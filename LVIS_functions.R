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
	nVIS<-matrix(0,length(all_expts),1);
	top20_IS_freq_vs_expts<-matrix(0,21,length(all_expts));
	colnames(top20_IS_freq_vs_expts)<-all_expts;
	rownames(nVIS)<-all_expts;
	colnames(nVIS)<-'nVIS';
	rownames(top20_IS_freq_vs_expts)<-c(c(1:20),'rest');
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
	res[['sample_entropy']]<-sample_entropy;#disad is just on prob. not exactly the number of species
	res[['sample_Chao']]<-sample_Chao;#give a estimate of total number
	res[['sample_UC50']]<-sample_UC50;#how many species to cover half of abuadance
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

#this function is unnecessary. the use of scientific=FALSE in the format() function fix the problem. 
validate_coordinates <- function(all_P5_collected_IS){
	z1<-as.character(all_P5_collected_IS$start);
	z2<-as.character(all_P5_collected_IS$end);
	if (isempty(grep("e",z1))){
		print("start ok");
	} else{
		print(paste0("check start:",as.character(grep("e",z1))));
	}
	if (isempty(grep("e",z2))){
		print("end ok");
	} else{
		print(paste0("check end:",as.character(grep("e",z2))));
	}
}


preprocess_all_VIS_collection<- function(all_P1_collected_IS){
	#all_P1_collected_IS is the collection of output from the same patient
	all_IS_P1<-paste0(all_P1_collected_IS$seqnames,'.',format(all_P1_collected_IS$start,trim=TRUE,scientific=FALSE),'.',format(all_P1_collected_IS$end,trim=TRUE,scientific=FALSE),'.',all_P1_collected_IS$strand);
	all_IS_P1_aux<-paste0(all_P1_collected_IS$seqnames,':',format(all_P1_collected_IS$start,trim=TRUE,scientific=FALSE),'-',format(all_P1_collected_IS$end,trim=TRUE,scientific=FALSE));
	
	#20152
	u_all_IS_P1<-unique(all_IS_P1);
	u_all_IS_P1<-read.table(text = u_all_IS_P1, sep = ".");
	library(bedr);
	x<-paste0(u_all_IS_P1[,1],":",u_all_IS_P1[,2],'-',u_all_IS_P1[,3]);
	is.x.valid  <- is.valid.region(x); #get rid of other chr
	x <- x[is.x.valid];
	u_all_IS_P1<-u_all_IS_P1[is.x.valid,];
	dim(u_all_IS_P1)

	u_all_IS_P1<-cbind(u_all_IS_P1[,c(1:3)],tmp1='.',tmp2='.',strand=u_all_IS_P1[,4]);
	colnames(u_all_IS_P1)<-c("chr","st","ed","tmp1","tmp2","strand");
	#rownames(u_all_IS_P1)<-paste0(u_all_IS_P1$chr,':',u_all_IS_P1$st,'-',u_all_IS_P1$ed);
	#x<-rownames(u_all_IS_P1);
	x<-paste0(u_all_IS_P1$chr,':',u_all_IS_P1$st,'-',u_all_IS_P1$ed);
	x.sort <- bedr.sort.region(x);

	i_map<-match(x.sort,x);
	u_all_IS_P1<-u_all_IS_P1[i_map,];
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
	#currently, we use d=0 in avoid VIS sites getting too big
	#the Lance pipeline tested for d=5,6,7,8 to merge in in the raw read level
	#it is different from here..
	u_all_IS_P1.merge<-bedr(
		engine = "bedtools", 
     	input = list(i = u_all_IS_P1), 
        method = "merge", 
        params = "-s -d 0 -c 5,4,6 -o sum,collapse,distinct"
        );

	u_all_IS_P1.merge<-cbind(u_all_IS_P1.merge,eff_ind=c(1:dim(u_all_IS_P1.merge)[1]));
	#note that column ind refer to how many sites in the input list that merged list covers
	#aux stores the label those sites (indexed by row # in u_all_IS_P1. the aux information is then used to generate th
	#output, u_VIS_merg and map

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
	
	#in the preprocess functions, a simple unique is performed and sites mapped to the non-canonical chr are removed
	#the result is u_all_IS_P1, an additional column is added to the input frame, it maps every site in the input to the u_all_IS. NA will be used
	#if the site cannot be mapped, because it is found in one of the non-canonical chr,
	res1<-preprocess_all_VIS_collection(all_P1_collected_IS);
	u_all_IS_P1<-res1[['u_all_IS_sort']];
	all_P1_collected_IS_mapped<-res1[['all_collected_IS_mapped']];


	#NB, essentially, all sites in the  overlapping unique sites were merged using bedtools. 
	#u_VIS_merge stores the merged list, the mapping between the merged list and the input list is 
	#defined in the array map
	res2<-merge_all_preprocessed_VIS_collection(u_all_IS_P1);

	u_VIS_merge<-res2$u_VIS_merge;
	u_VIS_to_merged_VIS<-res2$map;

	#head(all_P1_collected_IS)
	#head(all_IS_P1_collected_mapped)
	#iz<-head(all_IS_P1_collected_mapped)$u_ind;

	#the unique VIS index to the merged index..
	#u_VIS_to_merged_VIS[iz,]

	#after 2 steps mapping, the original list is mapped all_P1_collected_IS got a new merged id...listed in u_VIS_merge
	#(shown in data.frame X);
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
		aux[[pick_samples[i]]]<-rownames(counts_IS_expts)[which(counts_IS_expts[,pick_samples[i]]>0)];
	}
	return(aux);
}


collect_IS_per_celltypes <- function(counts_IS_expts){
	library(dplyr)
	library(tidyr)
	X<-colnames(counts_IS_expts)
	X<-data.frame(X);
	X<-X %>% separate(X,c("patients","type","time"),"_");
	rownames(X)<-colnames(counts_IS_expts);
	uT<-c('bulk','CD3','CD19','CD14CD15','CD56','CD34');
	#unique(X$type);
	i<-1;
	tt<-uT[i];
	iz<-which(X$type==tt);
	r<-rowSums(as.matrix(counts_IS_expts[,rownames(X)[iz]]));
	all_r<-data.frame(r);
	for (tt in uT[2:length(uT)]){
		iz<-which(X$type==tt);
		r<-rowSums(as.matrix(counts_IS_expts[,rownames(X)[iz]]));
		all_r<-cbind(all_r,r);
	}
	colnames(all_r)<-uT;
	return(all_r);
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

#VIS1, 2 are bed files like VIS1_before$overall
cor_genomicDensity <- function(VIS1,VIS2){

	all_chr<-c(paste0('chr',as.character(1:22)),'chrX','chrY');
	all_cc<-matrix(0,length(all_chr),1);
	rownames(all_cc)<-all_chr;
	colnames(all_cc)<-'correlation';
	for (cc in all_chr){
		iz1<-which(VIS1$chr==cc);
		iz2<-which(VIS2$chr==cc);
		iz<-intersect(iz1,iz2);
		x<-VIS1$pct[iz];
		y<-VIS2$pct[iz];
		all_cc[cc,1]<-cor(x,y);
	}
	
	return(all_cc);
}

#input a few bed files as a list. for each file, we have chr, start, end storing the VIS
#row names not matter..
#trace("circos.genomicDensity",edit=TRUE)
#max:4e-5
#https://stackoverflow.com/questions/53600926/how-do-you-add-track-label-in-r-circlize
draw_circos<-function(VIS_list,win_size,bed){

	library(circlize);
	#bed = generateRandomBed(nr = 50, fun = function(k) sample(letters, k, replace = TRUE))
	gap<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);gap[24]<-5;

	circos.par(gap.after=gap);
	circos.initializeWithIdeogram(plotType = NULL);
	#circos.genomicLabels(bed, labels.column = 4, side = "outside", cex=0.5,
    #	col = as.numeric(factor(bed[[1]])), line_col = as.numeric(factor(bed[[1]])))
	circos.genomicIdeogram()

	n<-length(VIS_list);
	track_names<-names(VIS_list);
	for (i in 1:n){
		circos.genomicDensity(VIS_list[[i]], col = c("#600000"), track.height = 0.1,window.size=win_size,ylim.force = FALSE);
		#circos.text(sector.index="chr1",track.index = 2*i,get.cell.meta.data("cell.xlim")-mean(get.cell.meta.data("cell.xlim"))/2,
        #    get.cell.meta.data("cell.ylim"), labels = track_names[i],facing = "clockwise", niceFacing = TRUE, adj = c(0,0),cex=0.5)
	}
	circos.clear();
}

#use for plotting the hotspots..
draw_circos_rainfall<-function(VIS_list){
	library(circlize);
	gap<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);gap[24]<-5;
	circos.par(gap.after=gap);
	circos.initializeWithIdeogram();
	circos.genoimcRainfall(VIS_list);
	#circos.genomicDensity(VIS_list[[i]], col = c("#600000"), track.height = 0.1,window.size=win_size,ylim.force = FALSE);
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

rep.col<-function(x,n){
   matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
}

reorder_cormat <- function(cormat){
	# Use correlation between variables as distance
	dd <- as.dist((1-cormat)/2)
	hc <- hclust(dd)
	cormat <-cormat[hc$order, hc$order]
	return(cormat);
}


simplify_homer_annotation <- function(u_VIS_merge_homer_annot,M){
	print("load hg19 annotation");
	if (missing(M)){
		load("./hg19.basic.annotation.update.Rdata");
	}
	
	x<-u_VIS_merge_homer_annot$Annotation;
	x<-gsub("3' UTR","3_UTR",x);
	x<-gsub("5' UTR","5_UTR",x);

	X<-M[x,c(8,9)];
	u_VIS_merge_homer_annot_modify<-cbind(u_VIS_merge_homer_annot,X);
	return(u_VIS_merge_homer_annot_modify);
}

get_geneFreq <- function(u_VIS_merge_homer_annot_modify){
	print('we drop sites located at the intergenic region');
	gl<-u_VIS_merge_homer_annot_modify[which(u_VIS_merge_homer_annot_modify$region_anot!='Intergenic'),]$Gene.Name;
	geneFreq<-table(gl)[order(table(gl),decreasing = TRUE)];
	return(geneFreq);
}

generate_GeneCloud <- function(geneFreq,out_file){
	#library(tagcloud);
	#let's just take the top genes cover 25% lof all VIS..
	library(wordcloud)
	set.seed(1234) # for reproducibility 
	df <- data.frame(gene = names(geneFreq),freq=geneFreq);
	df<-df[,c(1,3)];
	colnames(df)<-c("gene","freq");
	df$gene<-as.character(df$gene);
	cf<-cumsum(df$freq)/sum(df$freq);
	df2<-df[1:which(cf>.25)[1],];
	top_n_gene<-dim(df2)[1];
	minf<-df[top_n_gene,2]-1;
	max_font<-df2[1,2]/sum(df$freq);
	max_font<-round(max_font*500,1);
	min_font<-df2[top_n_gene,2]/sum(df$freq);
	min_font<-round(min_font*500,1);
	dev.new(width=5, height=5, unit="in")
	wordcloud(words = df2$gene, freq = df2$freq,min.freq=minf,max.words=top_n_gene, random.order=FALSE,rot.per=0,fixed.asp=TRUE,colors=brewer.pal(8, "Dark2"),scale=c(max_font,min_font))
	#tagcloud(names(geneFreq),geneFreq,col="red",sel=1:top_n_gene,algorithm= "oval");
	dev.copy2pdf(file=out_file, out.type= "cairo" );
}

run_fgsea_from_geneFreq<-function(geneFreq,m_list){

	#m_list can be generated by msigdb
	#H	hallmark gene sets  are coherently expressed signatures derived by aggregating many MSigDB gene sets to represent well-defined biological states or processes.
	#C1	positional gene sets  for each human chromosome and cytogenetic band.
	#C2	curated gene sets  from online pathway databases, publications in PubMed, and knowledge of domain experts.
	#C3	motif gene sets  based on conserved cis-regulatory motifs from a comparative analysis of the human, mouse, rat, and dog genomes.
	#C4	computational gene sets  defined by mining large collections of cancer-oriented microarray data.
	#C5	GO gene sets  consist of genes annotated by the same GO terms.
	#C6	oncogenic gene sets  defined directly from microarray gene expression data from cancer gene perturbations.
	#C7	immunologic gene sets  defined directly from microarray gene expression data from immunologic studies.
	#library(msigdbr);
	#m_df = msigdbr(species = "Homo sapiens", category = "C5");
	#m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

	library(fgsea);
  	#library(annotables);
  	fgseaRes <- fgsea(m_list,geneFreq,minSize=15, maxSize=1000, nperm=10000);
  	fgseaRes<-data.frame(fgseaRes);

  	#fgseaRes<-fgseaRes[order(fgseaRes$ES,decreasing=TRUE),]
  	rownames(fgseaRes)<-fgseaRes$pathway;
  	return(fgseaRes);  
}


get_hotspots_thres<-function(X,nVIS){
#X is the output from get_genomicDensity
	window=X$end[1]-X$start[1]+1;
	genome_size=3.235e9;
	n<-genome_size/window;
	lambda<-nVIS/n;
	count=X$pct*window;
	all_Pvalue<-1-ppois(count,lambda);
	X$P<-all_Pvalue;
	X$P_Bonferroni<-X$P*dim(X)[1];
	X$P_Bonferroni[X$P_Bonferroni>1]<-1;
	X<-X[order(X$P),];
	library(sgof);
	res<-BH(X$P);
	X$P.adjust<-res$Adjusted.pvalues;
	rownames(X)<-paste0(X[,1],":",X[,2],"-",X[,3]);
	library(bedr);
	sort.regions <- bedr.sort.region(rownames(X));
	Xout<-X[sort.regions,];
	return(Xout);
}

convert_bed_format<-function(P1_hotspots.sort){
	library(dplyr)
	library(tidyr)
	df<-data.frame(P1_hotspots.sort);
	colnames(df)<-'region';
	df <- df %>% separate(region, c("chr", "tmp"),":")
	df <- df %>% separate(tmp, c("start", "end"),"-")
	return(df);
}

generate_cormat_heatmap<-function(cormat,orient){

	library(reshape2)
	cormat <- reorder_cormat(cormat);
	upper_tri <- get_upper_tri(cormat)
	# Melt the correlation matrix
	melted_cormat <- melt(upper_tri, na.rm = TRUE)
	melted_cormat$value2<-round(melted_cormat$value*100)/100;

	if (orient=='upper'){
		# Create a ggheatmap
		ggheatmap <- ggplot(melted_cormat, aes(Var1, Var2, fill = value))+
 		geom_tile(color = "white")+
 		scale_fill_gradient2(low = "white", high = "red", mid = "orange", 
   		midpoint = .7, limit = c(0.15,1), space = "Lab", 
    		name="Pearson\nCorrelation") +
  		theme_minimal()+ # minimal theme
 		theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    		size = 14, hjust = 1))+
 		coord_fixed()
		# Print the heatmap
		# print(ggheatmap)
		ggheatmap<-ggheatmap + 
		geom_text(aes(Var1, Var2, label = value2), color = "black", size = 5) +
		theme(
	  	axis.title.x = element_blank(),
	  	axis.title.y = element_blank(),
	  	panel.grid.major = element_blank(),
	  	panel.border = element_blank(),
	  	panel.background = element_blank(),
	  	axis.ticks = element_blank(),
	  	legend.justification = c(1, 0),
	  	legend.position = c(0.8, 0.2),
	  	legend.direction = "horizontal")+
	  	guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
        	        title.position = "top", title.hjust = 0.5));
	} else if (orient=='lower'){
		ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
 		geom_tile(color = "white")+
 		scale_fill_gradient2(low = "white", high = "red", mid = "orange", 
   		midpoint = .7, limit = c(0.15,1), space = "Lab", 
    		name="Pearson\nCorrelation") +
  		theme_minimal()+ # minimal theme
 		theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    		size = 14, hjust = 1))+
 		coord_fixed()
		# Print the heatmap
		# print(ggheatmap)
		ggheatmap<-ggheatmap + 
		geom_text(aes(Var2, Var1, label = value2), color = "black", size = 5) +
		theme(
	  	axis.title.x = element_blank(),
	  	axis.title.y = element_blank(),
	  	panel.grid.major = element_blank(),
	  	panel.border = element_blank(),
	  	panel.background = element_blank(),
	  	axis.ticks = element_blank(),
	  	legend.justification = c(1, 0),
	  	legend.position = c(0.6, 0.7),
	  	legend.direction = "horizontal")+
	  	guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
        	        title.position = "top", title.hjust = 0.5));
	}

	return(ggheatmap);
}

run_bedtools_closest<-function(overlap_hotspots,all_SE_bed){
	aux1<-convert_bed_format(overlap_hotspots);
	aux2<-convert_bed_format(all_SE_bed);
	write.table(aux1,file='hotspots.bed',sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE);
	write.table(aux2,file="SE.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE);
	cmd<-paste0("bedtools closest -a hotspots.bed -b SE.bed > res_out");
	system(cmd);
	#bedtools closest -a vis.bed -b SE.bed > res_out
	res<-read.table('res_out');
	tmp<-paste0(res[,1],":",res[,2],"-",res[,3]);
	d1<-abs(res[,2]-res[,5]);
	d2<-abs(res[,2]-res[,6]);
	d3<-abs(res[,3]-res[,5]);
	d4<-abs(res[,3]-res[,6]);
	out<-cbind(res,minD=pmin(d1,d2,d3,d4));
	is.overlap <- in.region(tmp,all_SE_bed);
	out[which(is.overlap),]$minD<-0;
	#res<-cbind(res,overlap=is.overlap);
	iz<-which(!duplicated(tmp));
	out<-out[iz,];
	rownames(out)<-paste0(out[,1],":",out[,2],"-",out[,3]);
	return(out);
}




