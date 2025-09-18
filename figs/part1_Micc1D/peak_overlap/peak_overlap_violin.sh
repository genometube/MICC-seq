# narrowpeaks
zcat /research/xieyeming1/db/hek293/dnase/encode/ENCFF013WVF.bed.gz |awk '$9>15' > hek_dnase.narrowpeak
ln -s /research/xieyeming1/proj_2025/MICC_paper/peak_hg19/peak_overlap_1d/HekMiccEcoG2h/HekMiccEcoG2h.narrowpeak .

bedtools="/research/xieyeming1/software/Miniconda/envs/xieyeming1/bin/bedtools"
${bedtools} intersect -a hek_dnase.narrowpeak -b HekMiccEcoG2h.narrowpeak -v |cut -f1-3,9 > hek_dnase_non_overlap_pval.bdg
${bedtools} intersect -a hek_dnase.narrowpeak -b HekMiccEcoG2h.narrowpeak |cut -f1-3,9 > hek_dnase_overlap_pval.bdg

cat hek293_atac_medium_depth_peaks.narrowPeak > hek_atac.narrowpeak_
${bedtools} intersect -a hek_atac.narrowpeak_ -b HekMiccEcoG2h.narrowpeak -v|cut -f1-3,9 > hek_atac_non_overlap_pval.bdg
${bedtools} intersect -a hek_atac.narrowpeak_ -b HekMiccEcoG2h.narrowpeak |cut -f1-3,9 > hek_atac_overlap_pval.bdg

source /mnt/software/anaconda3/bin/activate R4_4
Rscript -e "rmarkdown::render('hek_atac_overlap_pval.Rmd', 
output_dir = '/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/figs/section1/peak_overlap')"

Rscript -e "rmarkdown::render('hek_dnase_overlap_pval.Rmd', 
output_dir = '/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/figs/section1/peak_overlap')"
