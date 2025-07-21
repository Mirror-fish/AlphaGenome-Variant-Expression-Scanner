# AlphaGenome Variant Expression Scanner
This is an established tool using Alphagenome to screen the possible effects of variants on expression patterns.

## 摘要  Summary 

**中文:** 输入一个含变异列表的文件 (TSV/CSV/VCF)。对指定 UBERON 组织调用 AlphaGenome DNA 模型的 RNA 预测，在变异中心 ±`scan_span` bp 区域内，用 `window_size` 滑窗计算 ALT/REF−1；根据 `threshold`、`min_length`、`merge_distance` 判定显著表达改变区段；输出结果汇总表 + (可选) 图像。

**English:** Input a file containing a list of variants (TSV/CSV/VCF). Perform RNA prediction using the AlphaGenome DNA model for the specified UBERON organ, within the ±`scan_span` bp region around the variant center. Use a sliding window of size `window_size` to calculate ALT/REF−1. Significant expression change segments are determined based on `threshold`, `min_length`, and `merge_distance`. The results are output in a summary table + (optional) images.

## 使用示例 Case of using

```bash
python alphagenome_sv_expression_scan.py \
  --variants SVs_hg38.txt \
  --organs UBERON:0000992 UBERON:0000955 \
  --threshold 0.2 --min-length 2000 --merge-distance 300 \
  --window-size 100 --scan-span 100000 \
  --output-table results.csv --output-dir plots \
  --api-key your_API_key \
  --plot-non-sig --scan-all-tracks
```

## 开始之前的调试 Before you start

**Install Package:** alphagenome, pandas, numpy, matplotlib
```
pip install alphagenome
```

**API key:** Please see [Alphagenome's Github](https://github.com/google-deepmind/alphagenome?tab=readme-ov-file) before you start, you can get the API key there.

**UBERON organs:** We have set UBERON:0000992 (ovary),**UBERON:0002371 (bone marrow, NOTICE: no tracks avaliable)**, UBERON:0000948 (heart), UBERON:0000955 (brain), UBERON:0001264 (pancreas), UBERON:0001134 (skeletal muscle tissue); Find more on [Ontology Search](https://www.ebi.ac.uk/ols4/)

---

## 1. Parameters（CLI）

| Parameter            | Type/Default                             | Meaning                                              | Notes               |
| -------------------  | ---------------------------------------- | ---------------------------------------------------- | ------------------- |
| `--variants`         | path                                     | Variant table TSV/CSV/VCF. Must contain CHROM, POS (1-based), REF, ALT | Required            |
| `--organs`           | list                                     | UBERON codes; default uses the built-in example list | Multiple selection  |
| `--threshold`        | float=0.5                                | Determines significance: ALT/REF−1 > threshold. Lower values are more sensitive but increase false positives | Lower values are more sensitive but increase false positives |
| `--min-length`       | int=1000                                  | Segment length threshold (bp, calculated after window) | Denoising           |
| `--merge-distance`   | int=300                                   | Maximum merge distance for adjacent candidate segments (bp) | Fills small gaps    |
| `--window-size`      | int=100                                   | Sliding window size (bp)                             | Smoothness vs Resolution |
| `--scan-span`        | int=50000                                 | Number of bp to scan on either side of the variant center | Scan range          |
| `--plot-non-sig`     | flag                                     | Plot even when no significant results are found       | Default does not plot |
| `--scan-all-tracks`  | flag                                     | Scan all tracks (default summarizes only, scans all when this flag is present; keeps backward compatibility) | Recommended to enable |
| `--epsilon`          | float=1e-8                                | Adds a small value to REF to avoid division by zero   | Numerical stability |
| `--output-table`     | str=alphagenome_scan_results.csv         | Summary table path (extension determines format)      | Supports csv/tsv/xlsx |
| `--output-dir`       | str=alphagenome_scan_plots               | Output directory for plots                           | Automatically created |
| `--api-key`          | str                                      | AlphaGenome API key                                  | see [Alphagenome API](https://deepmind.google.com/science/alphagenome/) |
| `--gtf`              | str=GENCODE v46 feather                  | Annotation file                                       | Customizable        |
| `--chrom-col` etc.   | str                                      | Input column name mapping                             | Compatible with multiple source files |

## 2. Input

Input file format can be TSV/CSV/VCF. You can select which column to use in `--chrom-col`, `--pos-col`, `--ref-col`,`--alt-col`.

We recommend using TSV:
```
CHROM	POS	REF	ALT
chr1	908963	A	ATCG
chr1	49084104	GAGTC A
...
```

## 3. Output

The outcome is a **detailed table** (`--output-table`) and a **variant×organ summary table**（`_variant_organ_summary.*`）。

**detailed table columns：**

- `chrom`, `pos`, `ref`, `alt`
- `ontology` (UBERON)
- `track_name`
- `is_significant` (bool)
- `n_regions` (significant regions)
- `region_index` (0..n-1; no significant region: -1)
- `rel_start_bp`, `rel_end_bp` (relative position, bp)
- `abs_start_bp`, `abs_end_bp` (hg38 coodinate; pos + rel)
- `mean_score`, `max_score`, `min_score`
- `direction` (up/down/mixed)
- `plot_file` (if plot generated)

**variant×organ summary table columns：** every row showed variant×organ results; If any track is significant then `is_significant_any=True`，you can check `tracks_significant` col to see the significant tracks。

---
## 工具原理 How it works
This tool is built on top of AlphaGenome, with additional scanning methods organized to process the prediction results. 

We mainly use a sliding window (`--window-size`, default 100bp) to scan the two predicted tracks for change ratios. 

<img width="3832" height="1876" alt="image" src="https://github.com/user-attachments/assets/081b671a-3683-459b-a585-eeb7ace3a005" />

After filtering out noise from the significant change signals (increase/decrease > `--threshold`), we evaluate whether the variant significantly alters expression in a specific UBERON organ by checking if the significant change signal region exceeds `--min-length`.

<img width="3817" height="1257" alt="image" src="https://github.com/user-attachments/assets/5d4522d7-3110-4248-aa50-f7ecc98c86ba" />


---
## 未来可能的改进方面 Future directions:
The parameters are **not validated** at this time, so we strongly suggest you to test a few rounds before using it on large dataset. We are now working on the literature review for the best parameters.

Right now this tool is only applied on RNA-SEQ data prediction, we'll add more datatypes once we tested them.

We corrected the coordinated shift caused by large variant alternation (see [Alphagenome forum post](https://www.alphagenomecommunity.com/t/coordiate-change-when-predicting-long-variants/445)) with deletion on the reference track, so far this method will not affect prediction. We'll follow AlphaGenome's fix in future.

If there's any idea or suggestion, please contact Jingyu (zengjingyu@genomics.cn)

