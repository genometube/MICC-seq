# batchfile=($(ls -d /research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/bulk_biotin_hic/HekMiccEcoG*/|rev|cut -c2-|rev))
# sample_id=(${batchfile[@]##*/})
#sample_id=(${sample_id_tmp[@]%_1.fq.gz})
sample_id=($(cat id_mapping.txt|awk '{print $1}'))
echo ${sample_id[@]}

echo ${sample_id[@]}|sed 's/\ /\n/g'|parallel "sh juicer-metrics.sh {}"

