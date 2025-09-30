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
library(bedtoolsr)
library(cowplot)
library(gridGraphics) 

npg_colors <- pal_npg()(10)
hgd()
setwd('/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/figs/part1_Micc1D/scripts')
options(bedtools.path = "/research/xieyeming1/software/Miniconda/envs/fyt_py311/bin/")

hg19_10kBin<-'../../../db/peaks/hg19_windows_10000.bed'
hek_micc<-'../../../db/peaks/hekMiccHia5_3lanes.narrowpeak'
hek_dnase<-'../../../db/peaks/hek_dnase.narrowpeak.conservative'
hek_nipbl<-'../../../db/peaks/hek_nipbl.narrowpeak.conservative'

hek_micc_10kBin<-bedtoolsr::bt.intersect(a=hg19_10kBin,b=hek_micc,wa=TRUE,u=TRUE,sorted=TRUE)
hek_nipbl_10kBin<-bedtoolsr::bt.intersect(a=hg19_10kBin,b=hek_nipbl,wa=TRUE,u=TRUE,sorted=TRUE)
hek_dnase_10kBin<-bedtoolsr::bt.intersect(a=hg19_10kBin,b=hek_dnase,wa=TRUE,u=TRUE,sorted=TRUE)

source("../../../scripts/R_func/eulerr.R")
multi_overlap<-bedtoolsr::bt.multiinter(list(hek_micc_10kBin,hek_nipbl_10kBin,hek_dnase_10kBin),cluster=TRUE)
head(multi_overlap)
cols<-c('hek_micc_10kBin','hek_nipbl_10kBin','hek_dnase_10kBin')
p_euler_upset_hek <- plot_overlaps(multi_overlap,set_names=cols)
p_euler_upset_hek$euler

mm10_10kBin<-'../../../db/peaks/mm10_windows_10000.bed'
NIH3T3_micc<-'../../../db/peaks/3t3MiccHia5_3lanes.narrowpeak'
NIH3T3_dnase<-'../../../db/peaks/NIH3T3_dnase.narrowpeak'
NIH3T3_nipbl<-'../../../db/peaks/NIH3T3_nipbl.narrowpeak'

NIH3T3_micc_10kBin<-bedtoolsr::bt.intersect(a=mm10_10kBin,b=NIH3T3_micc,wa=TRUE,u=TRUE,sorted=TRUE)
NIH3T3_nipbl_10kBin<-bedtoolsr::bt.intersect(a=mm10_10kBin,b=NIH3T3_nipbl,wa=TRUE,u=TRUE,sorted=TRUE)
NIH3T3_dnase_10kBin<-bedtoolsr::bt.intersect(a=mm10_10kBin,b=NIH3T3_dnase,wa=TRUE,u=TRUE,sorted=TRUE)

source("../../../scripts/R_func/eulerr.R")
multi_overlap<-bedtoolsr::bt.multiinter(list(NIH3T3_micc_10kBin,NIH3T3_nipbl_10kBin,NIH3T3_dnase_10kBin),cluster=TRUE)
head(multi_overlap)
cols<-c('NIH3T3_micc_10kBin','NIH3T3_nipbl_10kBin','NIH3T3_dnase_10kBin')
p_euler_upset_NIH3T3 <- plot_overlaps(multi_overlap,set_names=cols)

##############################
hek_micc_dt<-fread('../../../db/peaks/hekMiccHia5_3lanes.narrowpeak')
# hek_nipbl_dt<-fread('/research/xieyeming1/db/hek293/tf_chip/nibpl/align/nibpl/nibpl_peaks.narrowPeak')
hek_nipbl_dt<-fread('/research/xieyeming1/proj_2025/MICC_paper/peak_hg19/deeptools/hekMiccEcoG_3lanes/NIPBL/NIPBL_hek.bed')

hg19_10kBin<-'../../../db/peaks/hg19_windows_10000.bed.gz'

hek_nipbl_dt_10kBin<-bedtoolsr::bt.intersect(a=hg19_10kBin,b=hek_nipbl_dt,wa=TRUE,u=TRUE,sorted=TRUE)
hek_micc_dt_10kBin<-bedtoolsr::bt.intersect(a=hg19_10kBin,b=hek_micc_dt,wa=TRUE,u=TRUE,sorted=TRUE)
dim(hek_nipbl_dt_10kBin)
merge(hek_nipbl_dt_10kBin,hek_micc_dt_10kBin,by=c('V1','V2','V3'))

source("../../../scripts/R_func/eulerr.R")
multi_overlap<-bedtoolsr::bt.multiinter(list(hek_micc_dt_10kBin,hek_nipbl_dt_10kBin),cluster=FALSE)
head(multi_overlap)
cols<-c('hek_micc_1d','hek_nipbl_peak')
p_euler_upset_hek <- plot_overlaps(multi_overlap,set_names=cols)
p_euler_upset_hek$euler
p_euler_upset_hek$upset

upset_hek_grob <- grid.grabExpr(draw(p_euler_upset_hek$upset))
upset_hek_grob





hek_dnase_overlap<-bedtoolsr::bt.intersect(a=hek_dnase_dt,b=hek_micc_dt,sorted=TRUE)
head(hek_dnase_overlap)
hek_dnase_non_overlap<-bedtoolsr::bt.intersect(a=hek_dnase_dt,b=hek_micc_dt,sorted=TRUE,v=TRUE)
p_violin_hek_micc_vs_dnase<-plot_comparison(hek_dnase_overlap, hek_dnase_non_overlap, 
                          plot_title = "hek_micc_vs_dnase", data1_lab='hek_dnase_overlap',data2_lab='hek_dnase_non_overlap',
                          y_axis = "log2_value",transform_type='log2',
                          feature_col = 9, args_bw = 0.3, fill_colors = c("indianred", "lightblue"),
                          y_limits_fold = c(0.8, 1.3)) 

upset_hek_grob <- grid.grabExpr(draw(p_euler_upset_hek$upset))
upset_NIH3T3_grob <- grid.grabExpr(draw(p_euler_upset_NIH3T3$upset))

png('../files/hek_3t3_overlap_gs.png',width = 18,height = 7,units = 'in',res = 300)
plot_grid(plot_grid(p_euler_upset_hek$euler, p_euler_upset_NIH3T3$euler, upset_hek_grob, upset_NIH3T3_grob, ncol = 2, labels = c("A", "B")),
  plot_grid(p_violin_hek_micc_vs_nipbl, p_violin_hek_micc_vs_dnase, ncol = 2, labels = c("C", "D")),nrow = 1)
dev.off()

pdf('../files/hek_3t3_overlap_gs.pdf',width = 18,height = 7)
plot_grid(plot_grid(p_euler_upset_hek$euler, p_euler_upset_NIH3T3$euler, upset_hek_grob, upset_NIH3T3_grob, ncol = 2, labels = c("A", "B")),
  plot_grid(p_violin_hek_micc_vs_nipbl, p_violin_hek_micc_vs_dnase, ncol = 2, labels = c("C", "D")),nrow = 1)
dev.off()



