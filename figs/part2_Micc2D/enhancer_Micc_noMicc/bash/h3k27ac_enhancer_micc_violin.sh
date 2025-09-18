bedtools=/research/xieyeming1/software/Miniconda/envs/xieyeming1/bin/bedtools
enhancer_MICC_stripe=/research/xieyeming1/proj_2025/MICC_paper/peak_hg19/peak_overlap_1d/hekMiccEcoG_3lanes/hekMiccEcoG_3lanes.narrowpeak
h3k4me1=/research/xieyeming1/proj_2025/MICC_paper/peak_hg19/misc/enhancer_gro_stripe_h3k4me1/h3k4me1.bed_

${bedtools} intersect -a $h3k4me1 -b $enhancer_MICC_stripe -v -wa -sorted > h3k4me1_noMicc.bed
${bedtools} intersect -a $h3k4me1 -b $enhancer_MICC_stripe -wa -sorted > h3k4me1_Micc.bed

h3k27ac_bw=/research/xieyeming1/db/hek293/histone_chip/h3k27ac/hek_h3k27ac.bw
/research/zhangchen/software/UCSC_application/bigWigToBedGraph $h3k27ac_bw h3k27ac.bedgraph
cat h3k27ac.bedgraph|awk '{print $1"\t"$2"\t"$3"\t"$4/($3-$2)}'|awk '$3-$2<1000' > h3k27ac.bedgraph_

cat hek_h3k27ac.200.bedgraph|sort -k1,1 -k2,2n > hek_h3k27ac.200.bedgraph_sorted
${bedtools} map -a h3k4me1_Micc.bed -b hek_h3k27ac.100.bedgraph_sorted -c 4 -o sum > h3k27ac_h3k4me1_Micc.bdg
${bedtools} map -a h3k4me1_noMicc.bed -b hek_h3k27ac.100.bedgraph_sorted -c 4 -o sum > h3k27ac_h3k4me1_noMicc.bdg

source /mnt/software/anaconda3/bin/activate R4_4
plot_title=h3k27ac_enhancer_micc_violin
outdir=/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/figs/enhancer_Micc_noMicc/
Rscript -e "rmarkdown::render('${plot_title}.Rmd', output_dir = '${outdir}')"

