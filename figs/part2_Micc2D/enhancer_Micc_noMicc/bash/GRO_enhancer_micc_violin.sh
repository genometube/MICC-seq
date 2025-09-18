bedtools=/research/xieyeming1/software/Miniconda/envs/xieyeming1/bin/bedtools

enhancer_MICC_stripe=/research/xieyeming1/proj_2025/MICC_paper/peak_hg19/peak_overlap_1d/hekMiccEcoG_3lanes/hekMiccEcoG_3lanes.narrowpeak
h3k4me1=/research/xieyeming1/proj_2025/MICC_paper/peak_hg19/misc/enhancer_gro_stripe_h3k4me1/h3k4me1.bed_

${bedtools} intersect -a $h3k4me1 -b $enhancer_MICC_stripe -v -wa -sorted > h3k4me1_noMicc.bed
${bedtools} intersect -a $h3k4me1 -b $enhancer_MICC_stripe -wa -sorted > h3k4me1_Micc.bed

cat /research/fangyitong/Project2022/wj_dpnIdpnII/GEOdata_newRNA/GROseq_merge.bdg | awk '{if(index($1,"_")==0) print $0}'|sort -k1,1 -k2,2n > GROseq_merge.bdg

${bedtools} map -a h3k4me1_Micc.bed -b GROseq_merge.bdg -c 4 -o max| awk '{if ($4 == ".") $4 = "0"; print}'|sed 's/ /\t/g' > GRO_h3k4me1_Micc.bdg
${bedtools} map -a h3k4me1_noMicc.bed -b GROseq_merge.bdg -c 4 -o max| awk '{if ($4 == ".") $4 = "0"; print}'|sed 's/ /\t/g' > GRO_h3k4me1_noMicc.bdg

source /mnt/software/anaconda3/bin/activate R4_4
plot_title=GRO_enhancer_micc_violin
outdir=/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/figs/enhancer_Micc_noMicc/
Rscript -e "rmarkdown::render('${plot_title}.Rmd', output_dir = '${outdir}')"

