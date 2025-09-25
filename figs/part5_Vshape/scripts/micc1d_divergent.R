library(ggplot2)
library(dplyr)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)
library("dendextend") 
library(data.table)
library(IRanges)
library(eulerr)
library(ggsci)
library(UpSetR)
library(ggsignif)
library(patchwork)
library(GenomicRanges)
library(BSgenome) 
library(Biostrings)
library(httpgd)
npg_colors <- pal_npg()(10)
hgd()
setwd('/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/figs/part5_Vshape/scripts')
options(bedtools.path = "/research/xieyeming1/software/Miniconda/envs/fyt_py311/bin/")

MICC_1d<-fread('../files/hekMiccEcoG_3lanes.narrowpeak',sep='\t',header = F)
divergent<-fread('../files/divergent_peaks_w_header.bed',sep='\t',header = T)
<<<<<<< HEAD
# colnames(divergent)[9]<-'p_value_vs_Local'
head(divergent)
# divergent$p_value_vs_Local<- -log10(divergent$p_value_vs_Local)
# summary(divergent$p_value_vs_Local)
# quantile_95 <-quantile(divergent$p_value_vs_Local, probs = 0.95)

# adjust all val > quantile_95 to quantile_95
# divergent$p_value_vs_Local[divergent$p_value_vs_Local > quantile_95] <- quantile_95
=======
colnames(divergent)[9]<-'p_value_vs_Local'

divergent$p_value_vs_Local<- -log10(divergent$p_value_vs_Local)
summary(divergent$p_value_vs_Local)
quantile_95 <-quantile(divergent$p_value_vs_Local, probs = 0.95)

# adjust all val > quantile_95 to quantile_95
divergent$p_value_vs_Local[divergent$p_value_vs_Local > quantile_95] <- quantile_95
>>>>>>> ccbc3adffb4561a004fd608049033dd22442bd07
divergent$focus_ratio<-divergent$focus_ratio+1
# plot density of divergent$p_value_vs_Local
ggplot(divergent, aes(x = focus_ratio)) +
  geom_density(fill = "indianred", alpha = 0.5) +
  theme_bw() +
  labs(title = "Density Plot of p_value_vs_Local", x = "p_value_vs_Local")

head(divergent)
head(MICC_1d)
divergent_MiccOverlap<-bedtoolsr::bt.intersect(a=divergent,b=MICC_1d,wa=TRUE,u=TRUE,sorted=TRUE)
divergent_NoMiccOverlap<-bedtoolsr::bt.intersect(a=divergent,b=MICC_1d,wa=TRUE,sorted=TRUE,v=TRUE)
summary(divergent$focus_ratio)
head(divergent_MiccOverlap)
dim(divergent_MiccOverlap)
dim(divergent_NoMiccOverlap)

source("../../../scripts/R_func/violin_compare.R")
feature_cols<-colnames(divergent)
<<<<<<< HEAD
cols_vec<-c(5,7,8)
=======
cols_vec<-c(5,7,8,9)
>>>>>>> ccbc3adffb4561a004fd608049033dd22442bd07
feature_col<-feature_cols[cols_vec]
plot_list <- list()

for (i in seq_along(cols_vec)) {
    col_idx <- cols_vec[i]
    col_name <- feature_cols[col_idx]
    
    p <- plot_comparison(
        data1_lab = 'divergent_MiccOverlap',
        data2_lab = 'divergent_NoMiccOverlap',
        overlap_data1 = divergent_MiccOverlap,
        overlap_data2 = divergent_NoMiccOverlap,
        plot_title = paste0("divergent_vs_Micc1D_", col_name),
        y_axis = paste0("log2_", col_name),
        feature_col = col_idx,
        y_limits_fold = c(0.5, 1.1),
        transform_type='log2',
        args_bw = 0.1,
        fill_colors = c("indianred", "lightblue")
    )
    plot_list[[i]] <- p
}

panel_plot <- wrap_plots(plot_list, ncol = 2, guides = "collect") & 
    theme(legend.position = "bottom")
print(panel_plot)




hek_atac_overlap<-bedtoolsr::bt.intersect(a=hek_atac_dt,b=hek_micc_dt,sorted=TRUE)
head(hek_atac_overlap)
hek_atac_non_overlap<-bedtoolsr::bt.intersect(a=hek_atac_dt,b=hek_micc_dt,sorted=TRUE,v=TRUE)
p_violin_hek_micc_vs_atac<-plot_comparison(hek_atac_overlap, hek_atac_non_overlap, 
                          plot_title = "hek_micc_vs_atac", data1_lab='hek_atac_overlap',data2_lab='hek_atac_non_overlap',
                          y_axis = "log2_value",transform_type='log2',
                          feature_col = 10, args_bw = 0.3, fill_colors = c("indianred", "lightblue"),
                          y_limits_fold = c(0.8, 1.3)) 

################### 250918 end ###################

################### prepare feature table ###################
# rna_exp<-fread('/research/xieyeming1/db/attribute_table/hg19_gene_exp/hg19_rna_seq.bed',sep='\t',header = F)
# head(rna_exp)
# dim(rna_exp)
# rna_exp$transcript<-gsub('\\.\\d+$','',rna_exp$V4)
# rna_exp$tx_length<-rna_exp$V3-rna_exp$V2

# rna_exp_clean<-rna_exp[,c(9,10,6,7)]
# colnames(rna_exp_clean)<-c('transcript','tx_length','rna_exp','gene_id')

MICC_1d<-fread('/research/xieyeming1/proj_2025/MICC_paper/peak_hg19/peak_overlap_1d/hekMiccEcoG_3lanes/hekMiccEcoG_3lanes.narrowpeak',sep='\t',header = F)
divergent_<-fread('divergent_peaks_w_header.bed',sep='\t',header = T)
head(divergent_)
# colnames(pausing_)<-c('transcript','chr','start','end','strand','Length','Copies','annotation','pausing_ratio','Promoter_reads','GeneBody_reads')
# keep the substring after last | in annotation
# pausing_$annotation<-gsub('.*\\|','',pausing_$annotation)
# table(pausing_$annotation)
dim(divergent_)
promoter_w<-200
# pausing<-merge(pausing_,rna_exp_clean,by='transcript')
# head(pausing[,c('chr','start','end','transcript')])
# dim(pausing)
# add genebody cpg density
# genebody_cpg<-calculate_motif_density(pausing[,c('chr','start','end','transcript')],motif='CG')
# colnames(genebody_cpg)<-c('chr','start','end','transcript','genebody_cpg_count','genebody_CpG')

################### paused/noPaused tss pausing index ###################
# percentile_threshold <- quantile(pausing$GeneBody_reads, 0.25)
# percentile_threshold
# pause_thres <- quantile(pausing$pausing_ratio, 0.25)
# pause_thres
# p0<-ggplot(pausing, aes(x=log10(GeneBody_reads))) + geom_density() + theme_bw() +
#  geom_vline(xintercept = log10(percentile_threshold), color = "red")
# p1<-ggplot(pausing, aes(x=log10(pausing_ratio))) + geom_density() + theme_bw() + 
# geom_vline(xintercept = log10(pause_thres), color = "red")+xlim(-1,3)

# library(gridExtra)
# grid.arrange(p0, p1, ncol = 2, widths = c(1, 1))
# # Print or save the combined plot
# combined_plot <- p1 + p0 + 
#     plot_layout(ncol = 2, widths = c(1, 1))  # Adjust widths as needed
# # Print or save the combined plot
# print(combined_plot)
head(divergent_bed)
pausing_bed<-pausing[,c('chr','start','end','transcript','pausing_ratio','strand','GeneBody_reads','Promoter_reads','tx_length','annotation','rna_exp')]
divergent_bed<-divergent_
divergent_tss_bed<-divergent_bed
pausing_tss_bed$end<-ifelse(pausing_tss_bed$strand=='+',pausing_tss_bed$start+promoter_w,pausing_tss_bed$end+promoter_w)
pausing_tss_bed$start<-ifelse(pausing_tss_bed$strand=='+',pausing_tss_bed$start-promoter_w,pausing_tss_bed$end-promoter_w)
pausing_tss_bed<-pausing_tss_bed[order(pausing_tss_bed$chr,as.numeric(pausing_tss_bed$start)),]
dim(pausing_tss_bed)
head(pausing_tss_bed)
table(pausing_tss_bed$annotation)
promoter_cpg<-calculate_motif_density(pausing_tss_bed[,c('chr','start','end','transcript')],motif='CG')
colnames(promoter_cpg)<-c('chr','start','end','transcript','promoter_cpg_count','promoter_CpG')
promoter_TATA<-calculate_motif_density(pausing_tss_bed[,c('chr','start','end','transcript')],motif='TATA')
colnames(promoter_TATA)<-c('chr','start','end','transcript','promoter_TATA_count','promoter_TATA')

pausing_tss_bed<-merge(pausing_tss_bed,promoter_cpg[,c(4,5,6)],by='transcript')
pausing_tss_bed<-merge(pausing_tss_bed,promoter_TATA[,c(4,5,6)],by='transcript')
pausing_tss_bed<-merge(pausing_tss_bed,genebody_cpg[,c(4,5,6)],by='transcript')
length(colnames(pausing_tss_bed))
pausing_tss_bed<-pausing_tss_bed[,c(2,3,4,1,5:17)]
head(pausing_tss_bed)
dim(pausing_tss_bed)

paused_tss_bed<-pausing_tss_bed[pausing_tss_bed$pausing_ratio>pause_thres&pausing_tss_bed$GeneBody_reads>percentile_threshold,]
paused_tss_bed<-paused_tss_bed[order(paused_tss_bed$chr,as.numeric(paused_tss_bed$start)),]
dim(paused_tss_bed)
noPaused_tss_bed<-pausing_tss_bed[pausing_tss_bed$pausing_ratio<pause_thres&pausing_tss_bed$GeneBody_reads>percentile_threshold,]
noPaused_tss_bed<-noPaused_tss_bed[order(noPaused_tss_bed$chr,as.numeric(noPaused_tss_bed$start)),]
dim(noPaused_tss_bed)

################### paused/noPaused tss overlapping with Micc1D ###################
multi_overlap<-bedtoolsr::bt.multiinter(list(paused_tss_bed,noPaused_tss_bed,MICC_1d),cluster=TRUE)
head(multi_overlap)
multi_overlap_mat<-multi_overlap[multi_overlap$V5!='1,2,3',6:length(colnames(multi_overlap))]
cols<-c('paused_tss_bed','noPaused_tss_bed','MICC_1d')
colnames(multi_overlap_mat)<-cols
head(multi_overlap_mat)
combinations <- apply(multi_overlap_mat, 1, function(row) {
    paste(colnames(multi_overlap_mat)[row == 1], collapse = "&")})
counts <- table(combinations)

set_counts <- as.numeric(counts)
names(set_counts) <- names(counts)
fit2 <- euler(set_counts)

p3<-plot(fit2, quantities = TRUE,     fills = npg_colors[c(1,9,2)])
p3

# pdf('euler_hekMiccEcoG_3lanes.pdf')
# print(p3)
# dev.off()

upset_data <- multi_overlap_mat
p2<-upset( data = as.data.frame(upset_data),
  sets = cols,
  order.by = "freq",
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Set Size")
p2

# pdf('upset_tss_pausing.pdf',width=5,height=5, onefile=F)
# print(p2)
# dev.off()

################### paused vs Micc1D numeric feature ###################
paused_tss_MiccOverlap<-bedtoolsr::bt.intersect(a=paused_tss_bed,b=MICC_1d,wa=TRUE,u=TRUE,sorted=TRUE)
paused_tss_NoMiccOverlap<-bedtoolsr::bt.intersect(a=paused_tss_bed,b=MICC_1d,wa=TRUE,sorted=TRUE,v=TRUE)
dim(paused_tss_MiccOverlap)
dim(paused_tss_NoMiccOverlap)
head(paused_tss_MiccOverlap)
head(paused_tss_NoMiccOverlap)

feature_cols<-colnames(paused_tss_bed)
cols_vec<-c(5,7,8,9,11,13,15,17)

feature_col<-feature_cols[cols_vec]
plot_list <- list()

for (i in seq_along(cols_vec)) {
    col_idx <- cols_vec[i]
    col_name <- feature_cols[col_idx]
    
    p <- plot_comparison(
        data1_lab = 'paused_tss_MiccOverlap',
        data2_lab = 'paused_tss_NoMiccOverlap',
        overlap_data1 = paused_tss_MiccOverlap,
        overlap_data2 = paused_tss_NoMiccOverlap,
        plot_title = paste0("paused_Micc1D_overlap_", col_name),
        y_axis = paste0("log2_", col_name),
        feature_col = col_idx,
        y_limits = c(-5, 10),
        transform_type='log2',
        args_bw = 0.3,
        fill_colors = c("indianred", "lightblue")
    )
    plot_list[[i]] <- p
}

panel_plot <- wrap_plots(plot_list, ncol = 4, guides = "collect") & 
    theme(legend.position = "bottom")
print(panel_plot)


################### noPaused vs Micc1D numeric feature ###################
noPaused_tss_MiccOverlap<-bedtoolsr::bt.intersect(a=noPaused_tss_bed,b=MICC_1d,wa=TRUE,u=TRUE,sorted=TRUE)
noPaused_tss_NoMiccOverlap<-bedtoolsr::bt.intersect(a=noPaused_tss_bed,b=MICC_1d,wa=TRUE,sorted=TRUE,v=TRUE)
feature_cols<-colnames(noPaused_tss_bed)
cols_vec<-c(5,7,8,9,11,13,15,17)
feature_col<-feature_cols[cols_vec]
plot_list <- list()

for (i in seq_along(cols_vec)) {
    col_idx <- cols_vec[i]
    col_name <- feature_cols[col_idx]
    
    p <- plot_comparison(
        data1_lab = 'noPaused_tss_MiccOverlap',
        data2_lab = 'noPaused_tss_NoMiccOverlap',
        overlap_data1 = noPaused_tss_MiccOverlap,
        overlap_data2 = noPaused_tss_NoMiccOverlap,
        plot_title = paste0("noPaused_Micc1D_overlap_", col_name),
        y_axis = paste0("log2_", col_name),
        feature_col = col_idx,
        y_limits = c(-5, 10),
        args_bw = 0.5,
        fill_colors = c("indianred", "lightblue")
    )
    plot_list[[i]] <- p
}

panel_plot <- wrap_plots(plot_list, ncol = 4, guides = "collect") & 
    theme(legend.position = "bottom")
print(panel_plot)

################### tss vs Micc1D numeric feature ###################
pausing_tss_bed<-pausing_tss_bed[order(pausing_tss_bed$chr,as.numeric(pausing_tss_bed$start)),]

tss_MiccOverlap<-bedtoolsr::bt.intersect(a=pausing_tss_bed,b=MICC_1d,wa=TRUE,u=TRUE,sorted=TRUE)
tss_NoMiccOverlap<-bedtoolsr::bt.intersect(a=pausing_tss_bed,b=MICC_1d,wa=TRUE,sorted=TRUE,v=TRUE)

dim(tss_MiccOverlap)
dim(tss_NoMiccOverlap)
colnames(pausing_tss_bed)
head(pausing_tss_bed)
feature_cols<-colnames(pausing_tss_bed)
cols_vec<-c(5,7,8,9,11,13,15,17)
feature_col<-feature_cols[cols_vec]
plot_list <- list()

for (i in seq_along(cols_vec)) {
    col_idx <- cols_vec[i]
    col_name <- feature_cols[col_idx]
    
    p <- plot_comparison(
        data1_lab = 'allTss_MiccOverlap',
        data2_lab = 'allTss_NoMiccOverlap',
        overlap_data1 = tss_MiccOverlap,
        overlap_data2 = tss_NoMiccOverlap,
        plot_title = paste0("allTss_Micc1D_overlap_", col_name),
        y_axis = paste0("log2_", col_name),
        feature_col = col_idx,
        y_limits_fold = c(0.5, 1.3),
        args_bw = 0.5,
        fill_colors = c("indianred", "lightblue")
    )
    plot_list[[i]] <- p
}

panel_plot <- wrap_plots(plot_list, ncol = 4, guides = "collect") & 
    theme(legend.position = "bottom")
print(panel_plot)

