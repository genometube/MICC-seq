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
# hgd()
setwd('/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/figs/section1/scripts')
options(bedtools.path = "/research/xieyeming1/software/Miniconda/envs/fyt_py311/bin/")

hg19_10kBin<-'../../../db/peaks/hg19_windows_10000.bed'
hekMicc1h<-'../../../db/peaks/HekMiccEcoG1h.narrowpeak'
hekMicc2h<-'../../../db/peaks/HekMiccEcoG2h.narrowpeak'
hekMicc3h<-'../../../db/peaks/HekMiccEcoG3h.narrowpeak'

hek_dnase<-'../../../db/peaks/hek_dnase.narrowpeak.conservative'
hek_atac<-'../../../db/peaks/hek_atac.narrowpeak.conservative'

hekMicc1h_10kBin<-bedtoolsr::bt.intersect(a=hg19_10kBin,b=hekMicc1h,wa=TRUE,u=TRUE,sorted=TRUE)
hekMicc2h_10kBin<-bedtoolsr::bt.intersect(a=hg19_10kBin,b=hekMicc2h,wa=TRUE,u=TRUE,sorted=TRUE)
hekMicc3h_10kBin<-bedtoolsr::bt.intersect(a=hg19_10kBin,b=hekMicc3h,wa=TRUE,u=TRUE,sorted=TRUE)

hek_atac_10kBin<-bedtoolsr::bt.intersect(a=hg19_10kBin,b=hek_atac,wa=TRUE,u=TRUE,sorted=TRUE)
hek_dnase_10kBin<-bedtoolsr::bt.intersect(a=hg19_10kBin,b=hek_dnase,wa=TRUE,u=TRUE,sorted=TRUE)

source("../../../scripts/R_func/eulerr.R")
multi_overlap<-bedtoolsr::bt.multiinter(list(hekMicc1h_10kBin,hekMicc2h_10kBin,hekMicc3h_10kBin,hek_atac_10kBin,hek_dnase_10kBin),cluster=FALSE)
head(multi_overlap)
cols<-c('hekMicc1h_10kBin','hekMicc2h_10kBin','hekMicc3h_10kBin','hek_atac_10kBin','hek_dnase_10kBin')
p_euler_upset_hek <- plot_overlaps(multi_overlap,set_names=cols)
p_euler_upset_hek$euler
p_euler_upset_hek$upset

upset_hek_grob <- grid.grabExpr(draw(p_euler_upset_hek$upset))

png('../files/hekMicc1h2h3h_overlap.png',width = 15,height = 8,units = 'in',res = 300)
plot_grid(p_euler_upset_hek$euler, upset_hek_grob, ncol = 2, labels = c("A", "B"))
dev.off()

pdf('../files/hekMicc1h2h3h_overlap.pdf',width = 15,height = 8)
plot_grid(p_euler_upset_hek$euler, upset_hek_grob, ncol = 2, labels = c("A", "B"))
dev.off()


