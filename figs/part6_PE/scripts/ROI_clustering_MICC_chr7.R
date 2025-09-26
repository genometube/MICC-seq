library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(zoo)
library(pheatmap)
library(Matrix)
library(data.table)
# library(strawr)
# library(LiebermanAidenHiC2009)
library(pheatmap)
library("dendextend")
library(ComplexHeatmap)
library(IRanges)
library(httpgd)
# library(umap)
library(circlize)  
library(cowplot)
library(gridGraphics) 
library(bedtoolsr)
options(bedtools.path = "/research/xieyeming1/software/Miniconda/envs/fyt_py311/bin/")

# hgd()
setwd('/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/figs/part6_PE/files/micc')

####################### 2D mtx #######################
micc_pos_count<-fread('MiccHEK_chr7.tsv.gz')

head(micc_pos_count)
mtx_max<-max(micc_pos_count$pos1,micc_pos_count$pos2)
mtx_min<-min(micc_pos_count$pos1,micc_pos_count$pos2)

positions <- c(mtx_min:mtx_max)
cartesian_join <- CJ(pos1 = positions, pos2 = positions)
mtx_dt <- merge(cartesian_join, micc_pos_count, by.x = c("pos1", "pos2"), 
                    by.y = c("pos1", "pos2"), all.x = TRUE)
mtx_dt[is.na(mtx_dt$count), "count"] <- 0
tail(mtx_dt)
ROI_L<-14800
ROI_R<-15000
micc_ROI<-mtx_dt[mtx_dt$pos1>=ROI_L & mtx_dt$pos1<=ROI_R & mtx_dt$pos2>=ROI_L & mtx_dt$pos2<=ROI_R,]
micc_ROI_df<-dcast(micc_ROI, pos1~pos2, value.var='count')
head(micc_ROI_df)

micc_ROI_mtx<-as.matrix(micc_ROI_df[,-1])
row.names(micc_ROI_mtx)<-micc_ROI_df$pos1
head(micc_ROI_mtx)
dim(micc_ROI_mtx)
# use complexheatmap to plot log10(micc_ROI_mtx+1)
log_mtx <- log10(micc_ROI_mtx + 1)

# Plot with ComplexHeatmap
p1<-Heatmap(log_mtx,
        show_heatmap_legend = FALSE,
        col = colorRamp2(c(0, 2), c("white", "red")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        heatmap_legend_param = list(title =paste0(ROI_L,'_',ROI_R,"\nlog10(count+1)")))
# show_heatmap_legend = FALSE,

# subsample 100000 rows from mtx_dt
dt_subsample<-mtx_dt[sample(1:nrow(mtx_dt),1000000),]

# ggplot density distribution of mtx_dt$count
p2<-ggplot(dt_subsample, aes(x = log10(count))) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Density Plot of count", y = "Density") +
  theme_minimal()

####################### 1D peak #######################
hek_miccPeak<-fread('../../../../db/peaks/hekMiccHia5_3lanes.narrowpeak')
colnames(hek_miccPeak)[1]<-c('V1')
hek_miccPeak_chr7<-hek_miccPeak[hek_miccPeak$V1=='chr7',]

### map promoter enhancer
promoter_293t_chr7<-fread('../../../../db/peaks/promoter_293t_chr7.bed')
enhancer_293t_chr7<-fread('../../../../db/peaks/enhancer_293t_chr7.bed')
head(hek_miccPeak)

hek_miccPeak_p<-bedtoolsr::bt.intersect(a=hek_miccPeak_chr7,b=promoter_293t_chr7,wao=TRUE,sorted=TRUE)
hek_miccPeak_pe<-bedtoolsr::bt.intersect(a=hek_miccPeak_p,b=enhancer_293t_chr7,wao=TRUE,sorted=TRUE)
head(hek_miccPeak_pe,100)
# remove duplicated row is V1 V2 V3 is duplicated
dim(hek_miccPeak_pe)
hek_miccPeak_pe<-hek_miccPeak_pe[!duplicated(hek_miccPeak_pe[,c('V1','V2','V3')]),]
# remove overlap both p and e
hek_miccPeak_pe<-hek_miccPeak_pe[!(hek_miccPeak_pe$V11=="chr7" & hek_miccPeak_pe$V15=="chr7"),]
# remove overlap neither p nor e
hek_miccPeak_pe<-hek_miccPeak_pe[!(hek_miccPeak_pe$V11=="." & hek_miccPeak_pe$V15=="."),]
head(hek_miccPeak_pe)

### adjacent loop pair, mapped with count
pos_pe<-cbind.data.frame(pos=(hek_miccPeak_pe$V2+hek_miccPeak_pe$V3)/2,PE_type=ifelse(hek_miccPeak_pe$V11=="chr7",'P','E'))
pos_pe$pos<-floor(pos_pe$pos/10000)
head(pos_pe)
pos_pe_adj<-cbind.data.frame(pos1=pos_pe$pos[-length(pos_pe$pos)],pos2=pos_pe$pos[-1],
                 PE_type=paste0(pos_pe$PE_type[-length(pos_pe$PE_type)],pos_pe$PE_type[-1]))
head(pos_pe_adj)
pos_pe_adj[pos_pe_adj$PE_type=='EP',3]<-'PE'
# head(pos_pe_adj_dt)
pos_pe_adj_dt<-as.data.table(pos_pe_adj)
pos_pe_adj_count<-mtx_dt[pos_pe_adj_dt, on = c('pos1', 'pos2'), nomatch = NULL]

p3<-ggplot(pos_pe_adj_count, aes(x = log10(count))) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "chr7 micc loop signal", y = "Density") +
  theme_classic()
pos_pe_adj_count_sig<-pos_pe_adj_count[pos_pe_adj_count$count>10,]

p4<-ggplot(pos_pe_adj_count_sig, aes(x = log10(count))) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "chr7 micc loop signal\n(log10(count+1) > 1)", y = "Density") +
  theme_classic()

#summarize features, loop span, count, h3k4me1, h3k27ac, pausing, nibpl, h3k36me33em63k3h ,lpbin ,gnisua cg ,ca72k3h ,1em4k3h ,tnuoc ,naps pool ,serutaef ezirammus#
head(pos_pe_adj_count_sig)
pos_pe_adj_count_sig$loop_span<-pos_pe_adj_count_sig$pos2-pos_pe_adj_count_sig$pos1+1

PE_type_freq<-as.data.frame(table(pos_pe_adj_count_sig$PE_type))

### plot features
# plot PE_type_freq
p5<-ggplot(PE_type_freq, aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  labs(title = "chr7 micc loop signal\n(log10(count+1) > 1)", y = "Frequency") +
  theme_classic()

# plot loop_span
p6<-ggplot(pos_pe_adj_count_sig, aes(x = log2(loop_span),group=PE_type,fill=PE_type)) +
  geom_density( alpha = 0.5) +
  labs(title = "chr7 micc loop span\n(log10(count+1) > 1)", y = "Density") +
  theme_classic()

# plot count
p7<-ggplot(pos_pe_adj_count_sig, aes(x = log2(count),group=PE_type,fill=PE_type)) +
  geom_density( alpha = 0.5) +
  labs(title = "chr7 micc count signal\n(log10(count+1) > 1)", y = "Density") +
  theme_classic()

p1_grob <- grid.grabExpr(draw(p1))  # 捕获绘图对象

dev.off()
png('micc_pe_signal.png',width=700,height=1000)
plot_grid(p1_grob, p2, p3, p4, p5, p6, p7, ncol = 2, 
          labels = c("A", "B", "C", "D", "E","F",'G'))
dev.off()
