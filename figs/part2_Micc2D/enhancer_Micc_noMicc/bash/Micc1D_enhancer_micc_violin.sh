bedtools=/research/xieyeming1/software/Miniconda/envs/xieyeming1/bin/bedtools
enhancer_MICC_stripe=/research/xieyeming1/proj_2025/MICC_paper/peak_hg19/peak_overlap_1d/hekMiccEcoG_3lanes/hekMiccEcoG_3lanes.narrowpeak
h3k4me1=/research/xieyeming1/proj_2025/MICC_paper/peak_hg19/misc/enhancer_gro_stripe_h3k4me1/h3k4me1.bed_

# use h3k4me1 peak overlaped enhancer
${bedtools} intersect -a $h3k4me1 -b $enhancer_MICC_stripe -v -wa -sorted > h3k4me1_noMicc.bed
${bedtools} intersect -a $h3k4me1 -b $enhancer_MICC_stripe -wa -sorted > h3k4me1_Micc.bed

bdg=/research/xieyeming1/proj_2025/MICC_paper/peak_hg19/peak_overlap_1d/hekMiccEcoG_3lanes/hekMiccEcoG_3lanes.bedgraph
${bedtools} map -a h3k4me1_Micc.bed -b ${bdg} -c 4 -o mean > Micc1D_h3k4me1_Micc.bdg
${bedtools} map -a h3k4me1_noMicc.bed -b ${bdg} -c 4 -o mean > Micc1D_h3k4me1_noMicc.bdg

source /mnt/software/anaconda3/bin/activate R4_4
plot_title=Micc1D_enhancer_micc_violin
outdir=/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/figs/enhancer_Micc_noMicc/
Rscript -e "rmarkdown::render('${plot_title}.Rmd', output_dir = '${outdir}')"

