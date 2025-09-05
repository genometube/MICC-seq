#!/bin/bash
echo ==========start at : `date` ==========
wd="/ifswh7/BC_PBG_R/fangyitong/wj_EcoGdpn_hic/0108bioC_linker/"
sample="HEK293T2h"
samtools="/ifswh4/BC_PUB_T1/Pipeline/RNA/RNA_RNAref/RNA_RNAref_2018a/software_new/samtools"
bamToBed="/ifswh7/BC_PBG_R/fangyitong/Software/anaconda3/envs/binf/bin/bamToBed"
ln -s /ifswh7/BC_PBG_R/fangyitong/archive/wj_EcoGdpn_reference/ reference &&\
ln -s /ifswh7/BC_PBG_R/fangyitong/archive/wj_dpnII_restriction/ restriction_sites &&\

###1. run juicer hic pipeline
bash /ifswh7/BC_PBG_R/fangyitong/Software/juicer/CPU/juicer.sh -g hg19 -d "$wd/$sample/" -s GATC -p "$wd/$sample/restriction_sites/hg19.chrom.sizes" -y "$wd/$sample/restriction_sites/hg19_DpnII.txt" -z "$wd/$sample/reference/hg19.fa" -D /ifswh7/BC_PBG_R/fangyitong/Software/juicer -b ligation -t 12 &&\
cd aligned &&\
rm abnormal.sam collisions_nodups.txt collisions.txt dups.txt merged_sort.txt opt_dups.txt unmapped.sam &&\

cd $wd/$sample/splits &&\
rm $sample.fastq.gz.sort.txt $sample.fastq.gz_unmapped.sam &&\
$samtools view -bS -h $sample.fastq.gz.sam > $sample.bam &&\
$samtools sort -o $sample.sort.bam $sample.bam &&\
$samtools index $sample.sort.bam &&\
rm $sample.fastq.gz.sam $sample.fastq.gz_abnorm.sam $sample.bam &&\
/ifswh7/BC_PBG_R/fangyitong/Software/anaconda3/envs/deeptools/bin/bamCoverage -b $sample.sort.bam -o $sample.bdg100 --outFileFormat bedgraph -bs 100 &&\
sort -k1,1 -k2,2n $sample.bdg100 > tmp &&\
mv tmp $sample.bdg100 &&\
/ifswh7/BC_PBG_R/fangyitong/Software/anaconda3/envs/deeptools/bin/bedGraphToBigWig $sample.bdg100 /ifswh7/BC_PBG_R/fangyitong/archive/hg19_index/hg19.fa.fai $sample.bdg100.bw &&\

echo ==========end at : `date` ========== && \
echo Still_waters_run_deep 1>&2 && \
echo Still_waters_run_deep > hic.sh.sign
