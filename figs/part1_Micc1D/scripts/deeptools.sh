computeMatrix=/mnt/software/anaconda3/envs/fyt/bin/computeMatrix
plotHeatmap=/mnt/software/anaconda3/envs/fyt/bin/plotHeatmap
# sample="$1"

sample=HekMiccEcoG1h2h3h
feature=atac
cat /research/xieyeming1/proj_2025/MICC_paper/peak_hg19/peak_overlap_1d/HekMiccEcoG1h/hek_atac.narrowpeak|cut -f1-3 > hek_atac.bed
peak=hek_atac.bed

ln -s /research/xieyeming1/proj_2025/MICC_paper/peak_hg19/juicer_track_vs_gs/HekMiccEcoG1h/HekMiccEcoG1h.bdg.bw HekMiccEcoG1h.bw
ln -s /research/xieyeming1/proj_2025/MICC_paper/peak_hg19/juicer_track_vs_gs/HekMiccEcoG2h/HekMiccEcoG2h.bdg.bw HekMiccEcoG2h.bw
ln -s /research/xieyeming1/proj_2025/MICC_paper/peak_hg19/juicer_track_vs_gs/HekMiccEcoG3h/HekMiccEcoG3h.bdg.bw HekMiccEcoG3h.bw
arr=($(ls HekMiccEcoG*.bw))
$computeMatrix reference-point \
-S ${arr[@]} \
-R ${peak} \
--referencePoint center \
-a 5000 -b 5000 --missingDataAsZero --skipZeros \
-out ${sample}_${feature}.tab.gz

$plotHeatmap \
 -m ${sample}_${feature}.tab.gz\
 -out ${sample}_${feature}.pdf \
 --heatmapHeight 15  \
 --refPointLabel hek.atac.center \
 --regionsLabel hek.atac \
 --plotTitle ${sample}_${feature}--colorMap Reds 

#  --zMax 400 --yMax 400 

#  /research/xieyeming1/db/hek293/atac/analysis/hek293_atac_medium_depth_peaks.narrowPeak.chr7.bed
# feature=dnase
# cat /research/xieyeming1/proj_2025/MICC_paper/peak_hg19/peak_overlap_1d/HekMiccEcoG1h/hek_dnase.narrowpeak|cut -f1-3 > hek_dnase.bed
# peak=hek_dnase.bed

# $computeMatrix reference-point \
# -S ../../juicer_track_vs_gs/${sample}/${sample}.bdg.bw \
# -R ${peak} \
# --referencePoint center \
# -a 5000 -b 5000 --missingDataAsZero \
# -out ${sample}_${feature}.tab.gz

# $plotHeatmap \
#  -m ${sample}_${feature}.tab.gz\
#  -out ${sample}_${feature}.pdf \
#  --heatmapHeight 15  \
#  --refPointLabel hek.atac.center \
#  --regionsLabel hek.atac \
#  --plotTitle ${feature} --colorMap Reds 

