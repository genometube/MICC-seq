## Predict MICC from ATAC

### Dataset access
#### HEK293/NIH3T3 ATAC-seq
[HEK293_atac.narrowPeak](hek293_atac_medium_depth_peaks.narrowPeak) <br>
[HEK293_atac.bw](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2902628)<br>
[NIH3T3_atac.narrowPeak](../../db/peaks/NIH3T3_atac.narrowpeak) <br>
[NIH3T3_atac.bw](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2796909)

#### HEK293/NIH3T3 MICC-seq
[HEK293_micc.parquet](files/MICC_HEKwt_hg19_10000.parquet) <br>
[NIH3T3_micc.parquet](files/MICC_NIH3T3_mm10_10000.parquet)

---
### Dataset description
The .parquet files contain chromatin interaction data converted from Juicer's merged_nodups.txt format. 
<br>Each record represents a non-duplicate valid pair, with pos1 and pos2 indicating genomic loci (binned at 10kb resolution).

| chr |  pos1 |  pos2 | count |
|-----|-------|-------|-------|
|   1 |  2053 | 11780 |     1 |
|   1 |  8025 |  8365 |     1 |
|   1 | 22452 | 24178 |     1 |

<br>
Quick check of .parquet file:<br>
To generate .mtx file from .parquet file, run: python scripts/micc_mtx.py <br>
To generate ROI heatmap from .mtx file, run: python scripts/mtx_heatmap.py <br>

<img src="files/MICC_HEKwt_chr7_10000_14800_15000.png" alt="ROI heatmap" width="300">

---
### Reference
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012136
这篇用扩散模型，用HFFc6细胞的Hi-C (x) 和 Micro-C (y) 数据做训练，然后在6个不同细胞中验证（输入hic，输出microC）

https://www.science.org/doi/10.1126/sciadv.adr8265
这篇用扩散模型，用12878的dnase seq和单细胞hic数据做训练，然后输入IMR90的dnase seq，可以预测IMR90的hic

---
### Rationale
MICC 1D signal shows high correlation with ATAC-seq, which suggests that ATAC-seq can be used to predict MICC. 
<img src="files/HekMiccEcoG_fullBlood_snHiChew_HEK293Twt_10000_148M_150M.png" alt="MICC vs Hi-C Comparison" width="500">
<img src="files/hek_3t3_overlap_gs.png" alt="MICC vs ATAC/DNase" width="500">