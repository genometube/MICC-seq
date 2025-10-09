date
cat out_xls/*xls |awk '$2>0' > juicer_metrics.xls.tmp
echo "sample Total_pairs_processed Unmapped_pairs Chimeric Reported_pairs PCR_dup Intra_fragment_pairs low_MAPQ Optical Chimeric_paired valid_pairs valid_interaction_rmdup trans_interaction cis_shortRange cis_longRange"| tr [:blank:] \\t > juicer_metrics.header
cat juicer_metrics.header juicer_metrics.xls.tmp > juicer_metrics.xls

Rscript='/mnt/software/anaconda3/envs/Bioinfo/bin/Rscript'
${Rscript} juicer.R
date
