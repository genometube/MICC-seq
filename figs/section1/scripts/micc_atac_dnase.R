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

npg_colors <- pal_npg()(10)
# hgd()
setwd('/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/figs/section1/scripts')
options(bedtools.path = "/research/xieyeming1/software/Miniconda/envs/fyt_py311/bin/")

hg19_w10k<-'../../../db/peaks/hg19_windows_10000.bed'
hek_micc<-'../../../db/peaks/HekMiccEcoG2h.narrowpeak'
hek_dnase<-'../../../db/peaks/hek_dnase.narrowpeak.conservative'
hek_atac<-'../../../db/peaks/hek_atac.narrowpeak.conservative'

hek_micc_w10k<-bedtoolsr::bt.intersect(a=hg19_w10k,b=hek_micc,wa=TRUE,u=TRUE,sorted=TRUE)
hek_atac_w10k<-bedtoolsr::bt.intersect(a=hg19_w10k,b=hek_atac,wa=TRUE,u=TRUE,sorted=TRUE)
hek_dnase_w10k<-bedtoolsr::bt.intersect(a=hg19_w10k,b=hek_dnase,wa=TRUE,u=TRUE,sorted=TRUE)

source("../../../scripts/R_func/eulerr.R")
multi_overlap<-bedtoolsr::bt.multiinter(list(hek_micc_w10k,hek_atac_w10k,hek_dnase_w10k),cluster=TRUE)
head(multi_overlap)
cols<-c('hek_micc_w10k','hek_atac_w10k','hek_dnase_w10k')
p_euler_upset_hek <- plot_overlaps(multi_overlap,set_names=cols)
p_euler_upset_hek$euler

mm10_w10k<-'../../../db/peaks/mm10_windows_10000.bed'
NIH3T3_micc<-'../../../db/peaks/3t3MiccHia5_3lanes.narrowpeak'
NIH3T3_dnase<-'../../../db/peaks/NIH3T3_dnase.narrowpeak'
NIH3T3_atac<-'../../../db/peaks/NIH3T3_atac.narrowpeak'

NIH3T3_micc_w10k<-bedtoolsr::bt.intersect(a=mm10_w10k,b=NIH3T3_micc,wa=TRUE,u=TRUE,sorted=TRUE)
NIH3T3_atac_w10k<-bedtoolsr::bt.intersect(a=mm10_w10k,b=NIH3T3_atac,wa=TRUE,u=TRUE,sorted=TRUE)
NIH3T3_dnase_w10k<-bedtoolsr::bt.intersect(a=mm10_w10k,b=NIH3T3_dnase,wa=TRUE,u=TRUE,sorted=TRUE)

source("../../../scripts/R_func/eulerr.R")
multi_overlap<-bedtoolsr::bt.multiinter(list(NIH3T3_micc_w10k,NIH3T3_atac_w10k,NIH3T3_dnase_w10k),cluster=TRUE)
head(multi_overlap)
cols<-c('NIH3T3_micc_w10k','NIH3T3_atac_w10k','NIH3T3_dnase_w10k')
p_euler_upset_NIH3T3 <- plot_overlaps(multi_overlap,set_names=cols)
p_euler_upset_NIH3T3$euler
##############################
hek_micc_dt<-fread('../../../db/peaks/HekMiccEcoG2h.narrowpeak')
hek_atac_dt<-fread('../../../db/peaks/hek293_atac_medium_depth_peaks.narrowPeak')
hek_dnase_dt<-fread('../../../db/peaks/hek_dnase.narrowpeak')

head(hek_micc_dt)
source("../../../scripts/R_func/violin_compare.R")
hek_atac_overlap<-bedtoolsr::bt.intersect(a=hek_atac_dt,b=hek_micc_dt,sorted=TRUE)
head(hek_atac_overlap)
hek_atac_non_overlap<-bedtoolsr::bt.intersect(a=hek_atac_dt,b=hek_micc_dt,sorted=TRUE,v=TRUE)
p_violin_hek_micc_vs_atac<-plot_comparison(hek_atac_overlap, hek_atac_non_overlap, 
                          plot_title = "hek_micc_vs_atac", data1_lab='hek_atac_overlap',data2_lab='hek_atac_non_overlap',
                          y_axis = "log2_value",transform_type='log2',
                          feature_col = 10, args_bw = 0.3, fill_colors = c("indianred", "lightblue"),
                          y_limits_fold = c(0.8, 1.3)) 
p_violin_hek_micc_vs_atac
hek_dnase_overlap<-bedtoolsr::bt.intersect(a=hek_dnase_dt,b=hek_micc_dt,sorted=TRUE)
head(hek_dnase_overlap)
hek_dnase_non_overlap<-bedtoolsr::bt.intersect(a=hek_dnase_dt,b=hek_micc_dt,sorted=TRUE,v=TRUE)
p_violin_hek_micc_vs_dnase<-plot_comparison(hek_dnase_overlap, hek_dnase_non_overlap, 
                          plot_title = "hek_micc_vs_dnase", data1_lab='hek_dnase_overlap',data2_lab='hek_dnase_non_overlap',
                          y_axis = "log2_value",transform_type='log2',
                          feature_col = 9, args_bw = 0.3, fill_colors = c("indianred", "lightblue"),
                          y_limits_fold = c(0.8, 1.3)) 
p_violin_hek_micc_vs_dnase

p_euler_upset_hek$upset
p_euler_upset_NIH3T3$upset
plot_grid(
  plot_grid(p_euler_upset_hek$euler, p_euler_upset_NIH3T3$euler, test1,test2, ncol = 2, labels = c("A", "B")),
  plot_grid(p_violin_hek_micc_vs_atac, p_violin_hek_micc_vs_dnase, ncol = 2, labels = c("C", "D")),
  nrow = 2
)
test2

png('../files/hek_3t3_overlap_gs.png',width = 8,height = 11,units = 'in',res = 300)
grid.newpage()

pushViewport(viewport(layout = grid.layout(3, 2, heights = unit(c(3, 3, 5), 'in'))))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_euler_upset_hek$euler, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p_euler_upset_NIH3T3$euler, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(p_euler_upset_hek$upset, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
draw(p_euler_upset_NIH3T3$upset, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
print(p_violin_hek_micc_vs_atac, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 2))
print(p_violin_hek_micc_vs_dnase, newpage = FALSE)
popViewport()

dev.off()



pdf('../files/hek_3t3_overlap_gs.pdf',width = 10,height = 15)
grid.newpage()

pushViewport(viewport(layout = grid.layout(3, 2, heights = unit(c(4, 3, 5), 'in'))))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_euler_upset_hek$euler, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p_euler_upset_NIH3T3$euler, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(p_euler_upset_hek$upset, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
draw(p_euler_upset_NIH3T3$upset, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
print(p_violin_hek_micc_vs_atac, newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 2))
print(p_violin_hek_micc_vs_dnase, newpage = FALSE)
popViewport()

dev.off()
