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
npg_colors <- pal_npg()(10)
hgd()
setwd('/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/figs/section1/scripts')
options(bedtools.path = "/research/xieyeming1/software/Miniconda/envs/fyt_py311/bin/")

hg19_w10k<-'/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/db/peaks/hg19_windows_10000.bed'
hek_micc<-'/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/db/peaks/HekMiccEcoG2h.narrowpeak'
hek_dnase<-'/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/db/peaks/hek_dnase.narrowpeak'
hek_atac<-'/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/db/peaks/hek_atac.narrowpeak'

hek_micc_w10k<-bedtoolsr::bt.intersect(a=hg19_w10k,b=hek_micc,wa=TRUE,u=TRUE,sorted=TRUE)
hek_atac_w10k<-bedtoolsr::bt.intersect(a=hg19_w10k,b=hek_atac,wa=TRUE,u=TRUE,sorted=TRUE)
hek_dnase_w10k<-bedtoolsr::bt.intersect(a=hg19_w10k,b=hek_dnase,wa=TRUE,u=TRUE,sorted=TRUE)

source("/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/scripts/R_func/eulerr.R")
multi_overlap<-bedtoolsr::bt.multiinter(list(hek_micc_w10k,hek_atac_w10k,hek_dnase_w10k),cluster=TRUE)
head(multi_overlap)
# multi_overlap_mat<-multi_overlap[multi_overlap$V5!='1,2,3',6:length(colnames(multi_overlap))]
cols<-c('hek_micc_w10k','hek_atac_w10k','hek_dnase_w10k')
plots <- plot_overlaps(multi_overlap,set_names=cols)
plots$euler
plots$upset

