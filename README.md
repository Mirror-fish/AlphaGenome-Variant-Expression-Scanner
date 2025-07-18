# Alphagenome----Variant-Effect-Screening-Method
This is an established tool using Alphagenome to screen the possible effects of variants on expression patterns.

## 摘要  Summary 

**中文:** 输入一个含变异列表的文件 (TSV/CSV/VCF)。对指定 UBERON 器官调用 AlphaGenome DNA 模型的 RNA 预测，在变异中心 ±`scan_span` bp 区域内，用 `window_size` 滑窗计算 ALT/REF−1；根据 `threshold`、`min_length`、`merge_distance` 判定显著表达改变区段；输出结果汇总表 + (可选) 图像。

**English:** Input a file containing a list of variants (TSV/CSV/VCF). Perform RNA prediction using the AlphaGenome DNA model for the specified UBERON organ, within the ±`scan_span` bp region around the variant center. Use a sliding window of size `window_size` to calculate ALT/REF−1. Significant expression change segments are determined based on `threshold`, `min_length`, and `merge_distance`. The results are output in a summary table + (optional) images.

## 使用示例 Case of using

```bash
python alpha_variant_scan.py \
  --variants SVs_hg38.txt \
  --organs UBERON:0000992 UBERON:0000955 \
  --threshold 0.5 --min-length 1000 --merge-distance 300 \
  --window-size 100 --scan-span 50000 \
  --output-table results.csv --output-dir plots \
  --plot-non-sig --scan-all-tracks
```
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
| `--api-key`          | str                                      | AlphaGenome API key                                  | Optional, auto retrieved from env/colab |
| `--gtf`              | str=GENCODE v46 feather                  | Annotation file                                       | Customizable        |
| `--chrom-col` etc.   | str                                      | Input column name mapping                             | Compatible with multiple source files |

## 2. Output

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


