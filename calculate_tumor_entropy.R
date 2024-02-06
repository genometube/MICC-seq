library(dplyr)

### 1. read hic files
ncm_df <- strawr::straw("KR","/research/fangyitong/Project2022/wj_dpnIdpnII/data_share/0116_NCM_DLP_LOVO/NCM_200U_IP_30.hic","1","1","BP",5000) %>% na.omit() 
dld_df <- strawr::straw("KR","/research/fangyitong/Project2022/wj_dpnIdpnII/data_share/0116_NCM_DLP_LOVO/DLP_200U_30.hic","1","1","BP",5000) %>% na.omit()
lovo_df <- strawr::straw("KR","LOVO.hic","1","1","BP",5000) %>% na.omit()

### 2. calculate contact frequence of each pair of interaction Pi
ncm_df <- ncm_df %>% group_by(x) %>% summarise(count_1D=sum(counts)) %>% ungroup() %>% mutate(Pi=count_1D/sum(count_1D))
dld_df <- dld_df %>% group_by(x) %>% summarise(count_1D=sum(counts)) %>% ungroup() %>% mutate(Pi=count_1D/sum(count_1D))
lovo_df <-lovo_df %>% group_by(x) %>% summarise(count_1D=sum(counts)) %>% ungroup() %>% mutate(Pi=count_1D/sum(count_1D))
#ncm_df <- ncm_df %>% mutate(Pi=counts/sum(counts))
#dld_df <- dld_df %>% mutate(Pi=counts/sum(counts))
#lovo_df <-lovo_df %>% mutate(Pi=counts/sum(counts))

### 3. calculate entropy
ncm_entropy <- -sum(ncm_df$Pi * log2(ncm_df$Pi))
dld_entropy <- -sum(dld_df$Pi * log2(dld_df$Pi))
lovo_entropy <- -sum(lovo_df$Pi * log2(lovo_df$Pi))

### 4. plot barplot
dd <- data.frame(tech=c("NCM460","DLD1","LOVO"),chr1_entropy=c(ncm_entropy, dld_entropy, lovo_entropy))
library(ggplot2)
pdf("0_tumor_entropy_chr1.pdf")
ggplot(data=dd, aes(x=tech, y=chr1_entropy, color=tech, fill=tech)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=format(round(chr1_entropy, 2), nsmall = 2)), vjust=1.6, color="black", size=3.5)+
  theme_classic()
dev.off()