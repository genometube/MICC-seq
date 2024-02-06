#!/usr/bin/env Rscript

#USAGE: Rscript find_extrusion.R -i merged_hg38_chr7.hic -f hic -c 7 -r 5000
#USAGE: Rscript find_extrusion.R -i hg38chr7sim.matrix -f matrix -c 7 -5 5000

#ACCEPT INPUT FORMAT: .hic or contact matrix (chrom1  start1  end1    chrom2  start2  end2    count),which is generated from cooler dump
#args = commandArgs(trailingOnly=TRUE)
library("optparse")
 
option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
              help=".hic file name or contact matrix file name (assumed to be the format same as cooler dump output)", metavar="character"),
    make_option(c("-f", "--format"), type="character", default="hic", 
              help="input file format, default is hic", metavar="character"),
    make_option(c("-c", "--chromosome"), type="character", default="1", 
              help="read a chromosome from hic file, attention: some hic files stores chromosome without the \"chr\" character, default is 1",metavar="character"),
    make_option(c("-r", "--resolution"), type="integer", default=5000,
	     help="resolution to read .hic file", metavar="number")
); 
 
opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser);
my_chr = args$chromosome
my_res = args$resolution
##### default hg19 chromosome sizes #####
chr_size <- list("1"="249250621","2"="243199373","3"="198022430","4"="191154276","5"="180915260","6"="171115067","7"="159138663","8"="146364022",
                 "9"="141213431","10"="135534747","11"="135006516","12"="133851895","13"="115169878","14"="107349540","15"="102531392",
                 "16"="90354753","17"="81195210","18"="78077248","19"="59128983","20"="63025520","21"="48129895","22"="51304566",
                 "M"="16571","X"="155270560","Y"="59373566")
##### Load libraries #####
library(strawr)
library(Rcpp)
library(dplyr)
library(ggplot2)
library(data.table)

##### Read hic file #####
if (args$format=="hic") {
	hic_df <- strawr::straw("KR", args$input, my_chr, my_chr, "BP", my_res)
	} else {
	hic_df <- read.csv(args$input,sep="\t",header=T) %>% select(start1,start2,count) %>% rename(x=start1,y=start2,counts=count)} 

############### filter contact distances ###############
df <- hic_df %>% filter(x != y & y-x<2000000)

##### identify statistically significant extrusions #####
threshold_q <- 0.95 
x_peak <- df %>% group_by(x) %>% summarise(n = n()) %>% filter(n>qnorm(threshold_q,mean(n),sd(n)))
y_peak <- df %>% group_by(y) %>% summarise(n = n()) %>% filter(n>qnorm(threshold_q,mean(n),sd(n)))

##### generate x extrusions #####
num_bin <- 20
bin_size <- my_res*num_bin
x_extr_contacts <- hic_df %>% inner_join(x_peak,by="x") %>% na.omit() %>% #overlap contact matrix and x peaks 
                   mutate(bin=cut(y,breaks=seq(0,as.numeric(chr_size[my_chr]),by=bin_size))) %>% #eg. bin the chr7(length 159,138,663bp) for every 100,000bp
                   group_by(x,bin) %>% summarise(count_y_in_bin=n()) %>% ungroup() %>% #count the number of contacts for each x peak
                   mutate(bin_end=as.numeric(bin)*bin_size) #replace bin range with bin end
# print(chr_size$my_chr)
##### calculate lengths and densities of extrusions #####
x_extr_len <- data.frame(chr1=character(),start1=character(),end1=character(),
			 chr2=character(),start2=character(),end2=character(),
                         extr_len=character(),extr_dens=character()) #create an empty dataframe

for (each_x in unique(x_extr_contacts$x)){                  #for each x extrusion
    num <- subset(x_extr_contacts,x==each_x)$count_y_in_bin #subset x_peak_contacts
    dens <- c()
    
    for (j in seq(1,length(num))){                        #loop through each x_peak_contact
        dens <- c(dens,num[j]/num_bin)
        #if there are less than 50% contacts in two consecutive bins, we think the extrusion ends
        #if j reaches the last element, we think the extrusion ends
        if ((j==length(num)) | j!=1 & (num[j]<(num_bin*0.5) & num[j+1]<(num_bin*0.5))){    
            start_tmp <- each_x                                 #extrusion start 
            end_tmp <- subset(x_extr_contacts,x==each_x)$bin_end[j] #extrusion end 
            len_tmp <- end_tmp - start_tmp                      #extrusion length 
            dens_tmp <- mean(dens)
            x_extr_len <- x_extr_len %>% add_row(chr1=my_chr,
                                                 start1=format(start_tmp,scientific=F,trim=T),
                                                 end1=format(start_tmp+my_res,scientific=F,trim=T),
						 chr2=my_chr,
                                                 start2=format(end_tmp,scientific=F,trim=T), 
                                                 end2=format(end_tmp+my_res,scientific=F,trim=T), 
                                                 extr_len=format(len_tmp,scientific=F,trim=T),
                                                 extr_dens=format(dens_tmp,scientific=F,trim=T))
            break
        }
    }
}

x_extr_len <- x_extr_len %>% mutate_at(c('start1','end1','start2','end2','extr_len'), as.integer) %>%
              mutate_at(c('extr_dens'), as.numeric) %>% na.omit() %>% filter(extr_len>500000)# | extr_dens>=0.5)
write.table(x_extr_len, file=paste0("extrusion_x_",my_chr,".tsv"), sep="\t", quote=F, row.names=F, col.names=F)

##### generate y extrusions #####
num_bin <- 20
bin_size <- my_res*num_bin
y_extr_contacts <- hic_df %>% inner_join(y_peak,by="y") %>% na.omit() %>% #overlap contact matrix and y peaks 
                   mutate(bin=cut(x,breaks=seq(0,as.numeric(chr_size[my_chr]),by=bin_size))) %>% #bin the chr7(length 159,138,663bp) for every 100,000bp
                   group_by(y,bin) %>% summarise(count_x_in_bin=n()) %>% ungroup() %>% #count the number of contacts for each y peak
                   mutate(bin_end=as.numeric(bin)*bin_size) #replace bin range with bin end

##### calculate lengths and densities of extrusions #####
y_extr_len <- data.frame(chr1=character(),start1=character(),end1=character(),
                         chr2=character(),start2=character(),end2=character(),
                         extr_len=character(),extr_dens=character()) #create an empty dataframe

for (each_y in unique(y_extr_contacts$y)){                  #for each y extrusion
    num <- subset(y_extr_contacts,y==each_y)$count_x_in_bin #subset y_peak_contacts
    dens <- c()
    
    for (j in seq(length(num),1)){                        #loop through each y_peak_contact
        dens <- c(dens,num[j]/num_bin)

        #if j reaches the last element, we think the extrusion ends
        if (j==1) {
            start_tmp <- each_y                                 #extrusion start 
            end_tmp <- subset(y_extr_contacts,y==each_y)$bin_end[j] #extrusion end 
            len_tmp <- start_tmp - end_tmp                      #extrusion length 
            dens_tmp <- mean(dens)
            y_extr_len <- y_extr_len %>% add_row(chr1=my_chr,
                                                 start1=format(start_tmp,scientific=F,trim=T),
                                                 end1=format(start_tmp+my_res,scientific=F,trim=T),
						 chr2=my_chr,
                                                 start2=format(end_tmp,scientific=F,trim=T), 
                                                 end2=format(end_tmp+my_res,scientific=F,trim=T), 
                                                 extr_len=format(len_tmp,scientific=F,trim=T),
                                                 extr_dens=format(dens_tmp,scientific=F,trim=T))
            break
          #if there are less than 50% contacts in two consecutive bins, we think the extrusion ends
        } else if (j!=length(num) & num[j]<(num_bin*0.5) & num[j-1]<(num_bin*0.5)){    
            start_tmp <- each_y                                 #extrusion start 
            end_tmp <- subset(y_extr_contacts,y==each_y)$bin_end[j] #extrusion end 
            len_tmp <- start_tmp - end_tmp                      #extrusion length 
            dens_tmp <- mean(dens)
            y_extr_len <- y_extr_len %>% add_row(chr1=my_chr,
                                                 start1=format(start_tmp,scientific=F,trim=T),
                                                 end1=format(start_tmp+my_res,scientific=F,trim=T),
						 chr2=my_chr,
                                                 start2=format(end_tmp,scientific=F,trim=T), 
                                                 end2=format(end_tmp+my_res,scientific=F,trim=T), 
                                                 extr_len=format(len_tmp,scientific=F,trim=T),
                                                 extr_dens=format(dens_tmp,scientific=F,trim=T))
            break
            }
    }
}
y_extr_len <- y_extr_len %>% mutate_at(c('start1','end1','start2','end2','extr_len'), as.integer) %>%
	      mutate_at(c('extr_dens'), as.numeric) %>% na.omit() %>% filter(extr_len>500000)# | extr_dens>=0.5)  
write.table(y_extr_len, file=paste0("extrusion_y_",my_chr,".tsv"), sep="\t", quote=F, row.names=F, col.names=F)

#xy_extr <- bind_rows(x_extr_len, y_extr_len) # %>% distinct(chr,start,end)
xy_extr <- bind_rows(x_extr_len, y_extr_len)  %>% arrange(chr1,start1,end1,chr2,start2,end2)
write.table(xy_extr, file=paste0("extrusion_xy_",my_chr,".tsv"), sep="\t", quote=F, row.names=F, col.names=F)
