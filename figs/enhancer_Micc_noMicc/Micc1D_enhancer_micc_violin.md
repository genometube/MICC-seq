Micc1D_enhancer_micc_violin
================
<yemingxie@gmail.com>
Thu Aug 21 11:05:04 2025

``` r
knitr::opts_chunk$set(echo = TRUE)

library(ggplot2)
library(ggsci)
library(httpgd)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyr)
library(ggsignif)
npg_colors <- pal_npg()(10)

args = commandArgs(trailingOnly=TRUE)
file_path=getwd()

setwd('/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/figs/enhancer_Micc_noMicc')
plot_title='Micc1D_enhancer_micc_violin'
y_axis='log2_bedgraph_val'
in_file1='Micc1D_h3k4me1_Micc.bdg'
in_file2='Micc1D_h3k4me1_noMicc.bdg'
in_file_1 <- read.table(paste0(file_path,'/',in_file1), sep="\t", header=F)
in_file_1$subgroup<-in_file1
in_file_2 <- read.table(paste0(file_path,'/',in_file2), sep="\t", header=F)
in_file_2$subgroup<-in_file2

# combine all data
in_file_1_<-in_file_1[sample(1:nrow(in_file_1),6000),]
# random sample rows from Micc1D_enhancer_micc
in_file_2_<-in_file_2[sample(1:nrow(in_file_2),6000),]
# combine all data
all_data <- rbind(in_file_1_, in_file_2_)
# all_data$feature_val <- ifelse(all_data$V4 > 0, log2(all_data$V4), NA)
all_data$feature_val <- log2(all_data$V4+1)
wilcox_result <- wilcox.test(in_file_2$V4, in_file_1$V4, 
                           alternative = "two.sided")
head(all_data)
```

    ##         V1        V2        V3       V4                subgroup feature_val
    ## 7012  chr6 150311118 150311691 172.3333 Micc1D_h3k4me1_Micc.bdg    7.437405
    ## 822   chr1 200285984 200286297  96.0000 Micc1D_h3k4me1_Micc.bdg    6.599913
    ## 2896 chr16  22824939  22825255 135.5000 Micc1D_h3k4me1_Micc.bdg    7.092757
    ## 3488 chr17  42787860  42788095  50.0000 Micc1D_h3k4me1_Micc.bdg    5.672425
    ## 3846 chr18  46479462  46479694 159.0000 Micc1D_h3k4me1_Micc.bdg    7.321928
    ## 719   chr1 157108902 157109151  77.0000 Micc1D_h3k4me1_Micc.bdg    6.285402

``` r
summ <- all_data %>%
  group_by(subgroup) %>%
  dplyr::summarize(n = n(), mean = round(mean(feature_val, na.rm = TRUE),2),
    max_val = round(max(feature_val, na.rm = TRUE),2),
    sd = sd(feature_val, na.rm = TRUE))
summ
```

    ## # A tibble: 2 × 5
    ##   subgroup                      n  mean max_val    sd
    ##   <chr>                     <int> <dbl>   <dbl> <dbl>
    ## 1 Micc1D_h3k4me1_Micc.bdg    6000  6.82    15.6 0.930
    ## 2 Micc1D_h3k4me1_noMicc.bdg  6000  4.78    11.7 0.753

``` r
levels(factor(all_data$subgroup))
```

    ## [1] "Micc1D_h3k4me1_Micc.bdg"   "Micc1D_h3k4me1_noMicc.bdg"

``` r
p1<-ggplot(all_data, aes(x=subgroup, y=feature_val, fill=subgroup)) +
  geom_violin(trim=FALSE, bw=0.3, na.rm = TRUE) +
  geom_boxplot(width=0.1, outlier.shape = NA, na.rm = TRUE) +
  geom_text(aes(label = paste0('N=',n), y = max(mean, na.rm = TRUE)), 
            data = summ, size=4, vjust = -2, hjust = 2) +
  geom_text(aes(label = paste0('mean=',mean), y = max(mean, na.rm = TRUE)), 
            data = summ, size=4, vjust = 4, hjust = 2) +
  scale_fill_manual(values=c("lightblue","indianred")) +
  theme_classic() +  geom_signif(comparisons = list(levels(factor(all_data$subgroup))), 
              test = "wilcox.test", map_signif_level = TRUE,
              y_position = max(all_data$feature_val, na.rm = TRUE) * 1.2) +
  labs(title = plot_title,
    subtitle = paste0("Wilcoxon p-value = ", format.pval(wilcox_result$p.value)),
    y = y_axis) +
  theme(axis.text = element_text(face='bold'),
    axis.title = element_text(face="bold"),plot.title = element_text(face="bold"),
    legend.position="top") + ylim(0,10)
print(p1)
```

    ## Warning: Removed 30 rows containing non-finite outside the scale range
    ## (`stat_signif()`).

    ## Warning: Removed 3 rows containing missing values or values outside the scale
    ## range (`geom_signif()`).

![](Micc1D_enhancer_micc_violin_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
pdf(paste0(plot_title,'.pdf'))
print(p1)
```

    ## Warning: Removed 30 rows containing non-finite outside the scale range (`stat_signif()`).
    ## Removed 3 rows containing missing values or values outside the scale range (`geom_signif()`).

``` r
dev.off()
```

    ## png 
    ##   2
