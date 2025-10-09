library(ggplot2)
library(gridExtra)
library(ggsci)
library(reshape2)
mypal <- pal_npg("nrc", alpha = 0.7)(9)
color_palette<-c(mypal[1],mypal[2],mypal[3],mypal[4],mypal[5],mypal[7],mypal[8],mypal[9],'darkred','darkgreen','darkblue')

raw_df<-read.table('juicer_metrics.xls',header=TRUE)
#raw_df<-raw_df[c(1,2,4,5,7,8),1:15]

# total reads QC
#stack_df<-rbind.data.frame(raw_df,raw_df_pcr)
stack_df<-raw_df
stack_df_total_reads<-raw_df[,c(1:5)]
p1<-ggplot(stack_df_total_reads,aes(x=sample,y=log10(Total_pairs_processed)))+
  geom_bar(stat = 'identity',alpha=0.8)+
  theme_classic()+scale_fill_manual(values='grey')+coord_flip()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#p1

stack_df_total_reads[,c(3:5)]<-stack_df_total_reads[,c(3:5)]/stack_df_total_reads$Total_pairs_processed
stack_df_total_reads<-stack_df_total_reads[,-2]
total_reads_melted <- melt(stack_df_total_reads, id.vars = 'sample')

p2<-ggplot(total_reads_melted, aes(x = sample, y = value, fill = variable)) +
  geom_bar(stat = 'identity',alpha=0.8) +
  labs(x = 'Samples', y = 'Percentage', title = 'Total reads QC') +
  theme_classic()+scale_fill_manual(values=color_palette)+
  coord_flip()+theme(legend.position = "top")+guides(fill = guide_legend(nrow = 2))

# p2

# reported pair QC
stack_df_validpair<-stack_df[,c(1,5:9,12)]
p3<-ggplot(stack_df_validpair,aes(x=sample,y=log10(Reported_pairs)))+
  geom_bar(stat = 'identity',alpha=0.8)+
  theme_classic()+scale_fill_manual(values='grey')+coord_flip()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# p3
stack_df_validpair[,c(3:7)]<-stack_df_validpair[,c(3:7)]/stack_df_validpair$Reported_pairs
stack_df_validpair<-stack_df_validpair[,-2]
validpair_melted <- melt(stack_df_validpair, id.vars = 'sample')

p4<-ggplot(validpair_melted, aes(x = sample, y = value, fill = variable)) +
  geom_bar(stat = 'identity',alpha=0.8) +
  labs(x = 'Samples', y = 'Percentage', title = 'Reported pair QC') +
  theme_classic()+scale_fill_manual(values=color_palette)+
  coord_flip()+theme(legend.position = "top")+guides(fill = guide_legend(nrow = 3))
# p4

# cis trans
stack_df_cistrans<-stack_df[,c(1,12:15)]
p5<-ggplot(stack_df_cistrans,aes(x=sample,y=log10(valid_interaction_rmdup)))+
  geom_bar(stat = 'identity',alpha=0.8)+
  theme_classic()+scale_fill_manual(values='grey')+coord_flip()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# p5
stack_df_cistrans[,c(3:5)]<-stack_df_cistrans[,c(3:5)]/stack_df_cistrans$valid_interaction_rmdup
stack_df_cistrans<-stack_df_cistrans[,-2]
cistrans_melted <- melt(stack_df_cistrans, id.vars = 'sample')

p6<-ggplot(cistrans_melted, aes(x = sample, y = value, fill = variable)) +
  geom_bar(stat = 'identity',alpha=0.8) +
  labs(x = 'Samples', y = 'Percentage', title = 'cis trans QC') +
  theme_classic()+scale_fill_manual(values=color_palette)+
  coord_flip()+theme(legend.position = "top")+guides(fill = guide_legend(nrow = 3))
# p6

plotlist<-list(p2,p1,p4,p3,p6,p5)
library(gridExtra)
pdf('bulk_QC_metrics.pdf',width=10,height=10)
grid.arrange(grobs=plotlist,nrow=3,ncol=2,widths=c(2,1))
dev.off()

png('bulk_QC_metrics.png',width=800,height=800)
grid.arrange(grobs=plotlist,nrow=3,ncol=2,widths=c(2,1))
dev.off()
