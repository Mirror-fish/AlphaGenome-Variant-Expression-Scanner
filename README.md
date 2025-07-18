# Alphagenome----Variant-Effect-Screening-Method
This is a established tool using Alphagenome to screen the possible effects of variants to expression patterns.
## 0. 快速摘要 / Quick Summary / Resumen breve / 概要

**中文:** 输入一个含变异列表的文件 (TSV/CSV/VCF)。对指定 UBERON 器官调用 AlphaGenome DNA 模型的 RNA 预测，在变异中心 ±`scan_span` bp 区域内，用 `window_size` 滑窗计算 ALT/REF−1；根据 `threshold`、`min_length`、`merge_distance` 判定显著表达改变区段；输出结果汇总表 + (可选) 图像。

---

## 1. 可配置参数（CLI 选项）

| 参数                  | 类型/默认                              | 含义                                               | 备注               |      |            |
| ------------------- | ---------------------------------- | ------------------------------------------------ | ---------------- | ---- | ---------- |
| `--variants`        | path                               | 变异表 TSV/CSV/VCF。需含 CHROM, POS(1‑based), REF, ALT | 必填               |      |            |
| `--organs`          | list                               | UBERON 代码；缺省使用内置示例列表                             | 可多选              |      |            |
| `--threshold`       | float=0.5                          | 判定显著：                                            | ALT/REF−1        | > 阈值 | 越低越敏感但假阳性↑ |
| `--min-length`      | int=1000                           | 区段长度阈值 (bp, 以窗口计后换算)                             | 去噪               |      |            |
| `--merge-distance`  | int=300                            | 相邻候选区段最大合并间隔 (bp)                                | 填补小空洞            |      |            |
| `--window-size`     | int=100                            | 滑窗长度 (bp)                                        | 平滑度 vs 分辨率       |      |            |
| `--scan-span`       | int=50000                          | 扫描变异中心两侧 bp 数                                    | 扫描范围             |      |            |
| `--plot-non-sig`    | flag                               | 无显著结果时仍绘图                                        | 默认不绘             |      |            |
| `--scan-all-tracks` | flag                               | 扫描所有轨道 (默认仅汇总，无此 flag 时仍扫描全部；保留向后兼容)             | 建议开              |      |            |
| `--epsilon`         | float=1e-8                         | REF 加小值避免除零                                      | 数值稳定性            |      |            |
| `--output-table`    | str=alphagenome\_scan\_results.csv | 汇总表路径 (扩展名决定格式)                                  | 支持 csv/tsv/xlsx  |      |            |
| `--output-dir`      | str=alphagenome\_scan\_plots       | 绘图输出目录                                           | 自动创建             |      |            |
| `--api-key`         | str                                | AlphaGenome API key                              | 可留空自动取 env/colab |      |            |
| `--gtf`             | str=GENCODE v46 feather            | 注释文件                                             | 可自定义             |      |            |
| `--chrom-col` 等     | str                                | 输入列名映射                                           | 兼容多源文件           |      |            |


## 2. 输出结果表字段说明

生成 **track 级细粒度表** (`--output-table`)；并自动生成一个 **variant×organ 汇总表**（同目录，文件名追加 `_variant_organ_summary.*`）。

**细粒度表列：**

- `chrom`, `pos`, `ref`, `alt`
- `ontology` (UBERON)
- `track_name`
- `is_significant` (bool)
- `n_regions` (本轨道显著区段数)
- `region_index` (0..n-1；若无显著则-1)
- `rel_start_bp`, `rel_end_bp` (相对变异位置，bp；负=上游)
- `abs_start_bp`, `abs_end_bp` (hg38 绝对坐标；近似 = pos + rel)
- `mean_score`, `max_score`, `min_score`
- `direction` (up/down/mixed)
- `plot_file` (若生成器官级覆盖图)

**汇总表列：** 每个 variant×organ 一行；若任一轨道显著则 `is_significant_any=True`，并附 `tracks_significant` 列 (逗号分隔轨道名)。

---

## 3. 使用示例

```bash
python alpha_variant_scan.py \
  --variants SVs_hg38.txt \
  --organs UBERON:0000992 UBERON:0000955 \
  --threshold 0.5 --min-length 1000 --merge-distance 300 \
  --window-size 100 --scan-span 50000 \
  --output-table results.csv --output-dir plots \
  --plot-non-sig --scan-all-tracks
```
