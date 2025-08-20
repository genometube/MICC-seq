RPII_enhancer_micc_violin
================
<yemingxie@gmail.com>
Wed Aug 20 18:15:46 2025

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
plot_title='RPII_enhancer_micc_violin'
y_axis='log2_MACS_q_val'
in_file1='RPII_h3k4me1_Micc.bdg'
in_file2='RPII_h3k4me1_noMicc.bdg'
in_file_1 <- read.table(paste0(file_path,'/',in_file1), sep="\t", header=F)
in_file_1$subgroup<-in_file1
in_file_2 <- read.table(paste0(file_path,'/',in_file2), sep="\t", header=F)
in_file_2$subgroup<-in_file2

# combine all data
in_file_1_<-in_file_1[sample(1:nrow(in_file_1),6000),]
# random sample rows from RPII_enhancer_micc
in_file_2_<-in_file_2[sample(1:nrow(in_file_2),6000),]
# combine all data
all_data <- rbind(in_file_1_, in_file_2_)
# all_data$feature_val <- ifelse(all_data$V4 > 0, log2(all_data$V4), NA)
all_data$feature_val <- all_data$V4
wilcox_result <- wilcox.test(in_file_2$V4, in_file_1$V4, 
                           alternative = "two.sided")
head(all_data)
```

    ##        V1        V2        V3        V4              subgroup feature_val
    ## 5001 chr2 183579990 183580361  1.266667 RPII_h3k4me1_Micc.bdg    1.266667
    ## 4880 chr2  99772144  99772587  4.092308 RPII_h3k4me1_Micc.bdg    4.092308
    ## 650  chr1 153650526 153651016  0.000000 RPII_h3k4me1_Micc.bdg    0.000000
    ## 635  chr1 151300493 151300884 21.555000 RPII_h3k4me1_Micc.bdg   21.555000
    ## 6099 chr5  14663961  14664277  8.690000 RPII_h3k4me1_Micc.bdg    8.690000
    ## 7064 chr6 158982185 158982621  2.850000 RPII_h3k4me1_Micc.bdg    2.850000

``` r
summ <- all_data %>%
  group_by(subgroup) %>%
  dplyr::summarize(n = n(), mean = round(mean(feature_val, na.rm = TRUE),2),
    max_val = round(max(feature_val, na.rm = TRUE),2),
    sd = sd(feature_val, na.rm = TRUE))
summ
```

    ## # A tibble: 2 × 5
    ##   subgroup                    n  mean max_val    sd
    ##   <chr>                   <int> <dbl>   <dbl> <dbl>
    ## 1 RPII_h3k4me1_Micc.bdg    6000  8.97   1386. 26.8 
    ## 2 RPII_h3k4me1_noMicc.bdg  6000  2.43    139.  3.71

``` r
levels(factor(all_data$subgroup))
```

    ## [1] "RPII_h3k4me1_Micc.bdg"   "RPII_h3k4me1_noMicc.bdg"

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
    legend.position="top") + ylim(0,20)
p1
```

    ## Warning: Removed 423 rows containing non-finite outside the scale range
    ## (`stat_signif()`).

    ## Warning: Removed 3 rows containing missing values or values outside the scale
    ## range (`geom_signif()`).

![](/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/figs/enhancer_Micc_noMicc/RPII_enhancer_micc_violin_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
pdf(paste0(plot_title,'.pdf'))
print(p1)
```

    ## Warning: Removed 423 rows containing non-finite outside the scale range (`stat_signif()`).
    ## Removed 3 rows containing missing values or values outside the scale range (`geom_signif()`).

``` r
dev.off()
```

    ## png 
    ##   2
