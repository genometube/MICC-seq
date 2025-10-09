# proj_id=wj_WtFlaTripto_20241225
# proj_dir="/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/bulk_biotin_hic/${proj_id}"

echo ==========START at `date`==========
#  bedtools="/research/xieyeming1/software/Miniconda/envs/xieyeming1/bin/bedtools"
#  gatc_bed="/research/xieyeming1/db/hic_pro_ref/hg19/hg19_gatc.bed"
#  frag_bed="/research/xieyeming1/db/hic_pro_ref/hg19/DpnII_resfrag_hg19.bed"

mkdir -p out_bed out_xls
sample="$1"
#juicer_results="/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/bulk_biotin_hic/${sample}/juicer/${sample}/aligned/inter_30.txt"
juicer_results="aligned/${sample}/inter_30.txt"
date

# mapping
Total_pairs_processed=($(cat ${juicer_results}|grep Sequenced|awk '{print $4}'|sed 's/\,//g'))
Unmapped_pairs=($(cat ${juicer_results}|grep Unmapped|awk '{print $2}'|sed 's/\,//g'))
Chimeric=($(cat ${juicer_results}|grep "Chimeric Ambiguous:"|awk '{print $3}'|sed 's/\,//g'))
Reported_pairs=($(cat ${juicer_results}|grep "Alignable"|awk '{print $4}'|sed 's/\,//g'))

# pair filter
PCR_dup=($(cat ${juicer_results}|grep "PCR"|awk '{print $3}'|sed 's/\,//g'))
Optical=($(cat ${juicer_results}|grep "Optical"|awk '{print $3}'|sed 's/\,//g'))
Intra_fragment_pairs=($(cat ${juicer_results}|grep "Intra-fragment"|awk '{print $3}'|sed 's/\,//g'))
low_MAPQ=($(cat ${juicer_results}|grep MAPQ|awk '{print $4}'|sed 's/\,//g'))
valid_interaction_rmdup=($(cat ${juicer_results}|grep "Contacts"|awk '{print $3}'|sed 's/\,//g'))
Chimeric_paired=($(cat ${juicer_results}|grep "Chimeric Paired:"|awk '{print $3}'|sed 's/\,//g'))

# dup, range
valid_pairs=($(awk "BEGIN {printf \"%.0f\n\", $valid_interaction_rmdup + $PCR_dup }"))

trans_interaction=($(cat ${juicer_results}|grep "Inter"|awk '{print $2}'|sed 's/\,//g'))
cis_shortRange=($(cat ${juicer_results}|grep "Short"|awk '{print $4}'|sed 's/\,//g'))
cis_longRange=($(cat ${juicer_results}|grep ">20Kb"|awk '{print $4}'|sed 's/\,//g'))

mapped_id=($(cat id_mapping.txt|awk -v a="${sample}" '$1==a'|awk '{print $2}'))
echo "${mapped_id} ${Total_pairs_processed} ${Unmapped_pairs} ${Chimeric} ${Reported_pairs} ${PCR_dup} ${Intra_fragment_pairs} ${low_MAPQ} ${Optical} ${Chimeric_paired} ${valid_pairs} ${valid_interaction_rmdup} ${trans_interaction} ${cis_shortRange} ${cis_longRange}"| tr [:blank:] \\t > out_xls/${mapped_id}_metrics.xls

echo ==========END at `date`==========
