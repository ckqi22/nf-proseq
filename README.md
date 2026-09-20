# nf-proseq — PRO-seq 分析流程（Nextflow DSL2）

PRO-seq 全流程自动化：从 fastq 出发，完成质控、比对、spike-in 归一化、基因体定量、
差异表达、Pol II 活性位点（单碱基 5' 端 / full-read 覆盖度）分布、TSS metagene、
暂停指数（Pausing Index）与报告打包。

---

## 1. 环境要求

| 依赖 | 说明 |
|---|---|
| Nextflow | DSL2，需 v25.x（`into/tap/multiMap/-log` 已移除的算子不再使用） |
| 软件 | bowtie2、samtools、bedtools、deeptools、featureCounts、Rscript、python3、fastp、fastq_screen |
| 运行方式 | 宿主裸跑（软件用绝对路径，见 `nextflow.config` 的 `params`）；可选 `singularity` profile |

软件路径、数据库路径统一在 [nextflow.config](nextflow.config) 的 `params` 里配置（`feature_counts`、`r`、`python`、`fastp` 等），
deeptools 等其余工具路径在各模块内写死绝对路径。

---

## 2. 快速开始

```bash
# 1) 填写配置文件
cp /workplace/chenkai/nf-proseq/params.yml /path/to/project_dir/

# 2) 激活nextflow环境
mamba activate nextflow

# 3) 运行
nextflow run /workplace/chenkai/nf-proseq -params-file params.yml

# 4) 断点续跑（nextflow.config 已默认 resume = true）
nextflow run /workplace/chenkai/nf-proseq -params-file params.yml -resume
```

---

## 3. 输入文件说明

### 3.1 samplesheet.csv

CSV 表头 `id,sample,group,r1,r2`。`r2` 为空即单端；`id`只是标识，无实际作用；`sample`为展示在图/表/报告的样本名；`group` 缺失默认 `unknown`；`id` 列可省略。

```csv
id,sample,group,r1,r2
SRR28785827,Ints11_R2,Ints11,/path/to/SRR28785827.fastq.gz
SRR28785828,Ints11_R1,Ints11,/path/to/SRR28785828.fastq.gz
SRR28785829,WT_R2,WT,/path/to/SRR28785829.fastq.gz
SRR28785830,WT_R1,WT,/path/to/SRR28785830.fastq.gz
```

- 必填列：`sample`、`r1`（缺则报错）。
- `r2`：双端填 read2 路径；单端留空。
- `group`：分组列，驱动 metagene 组图、差异分析、暂停指数分组。

### 3.2 参考基因组与注释（数据库）

流程通过 [bin/parse_config.py](bin/parse_config.py) 解析两个 INI 数据库文件（路径在 `nextflow.config`）：

- `species_config`：`/workplace/database/database_species_taxid_config.txt` —— species → 构建名映射。
- `information_config`：`/workplace/database/database_information_config.txt` —— 每个构建名 section 下给出
  `gtf`、`bowtie2_index`、`genome_fasta`、`gene_annotation`、`rRNA_index`。

`params.yml` 里 `species` + `build` 指定用哪套参考；`gtf`/`genome_fasta` 可直接给路径覆盖数据库。

### 3.3 spike-in（可选）

| 参数 | 作用 |
|---|---|
| `spike_genome` | `information_config` 里的 section 名（如 `dm6`），自动查 spike fasta |
| `spike_fasta` | 直接给 spike fasta 路径，覆盖 `spike_genome` 查找 |
| `spike_index` | 预构建的「主+spike 合并」bowtie2 索引前缀，给则跳过 concat+build |
| `spike_chroms` | spike 染色体名单文件（配合 `spike_index` 必填） |
| `spike_min_fraction` | spike 占比下限（`spike_count / total_mapped`），`<=` 此值报错退出；`0.0` = 仅拒绝 0 spike reads |

留空 = 完全跳过 spike-in（库大小 CPM 归一化）。

> `spike_count` 只数 read1（与下游信号口径一致，PE 不再把 read2 算进去）。开 spike 时任一占比不达标的样本会直接 fail 退出，而不是被下游 join 静默丢弃。

---

## 4. 关键参数（params.yml）

| 参数 | 默认 | 说明 |
|---|---|---|
| `species` / `build` | mouse / null | 参考基因组；`build` 覆盖 species 映射 |
| `sample_sheet` | samplesheet.csv | 样本表（相对启动目录） |
| `adapter` | `I` | 接头类型：`I`/`UMI`/`HT`/`SP` |
| `signal_mode` | single | 分析信号：`single`=单碱基 5' 端 / `full`=full-read 覆盖度 / `both`=都算 |
| `normalize_methods` | cpm,fpkm | profile 表归一化（cpm/fpkm/rpkm，逗号分隔） |
| `report` | true | 是否打包报告 |
| `compared_groups` | — | 差异比较组列表：`A, B, FC, P, paired/unpaired` |
| `spike_min_fraction` | 0.0 | spike-in reads 占比下限（`spike_count/total_mapped`），低于即报错 |

QC 阈值、TSS/暂停指数、metagene 窗口、差异阈值、信号表阈值等分别在 `threshold:`、`tss:`、`metagene:`、`diff:`、`signal_table:` 段配置。

`signal_table:` 段（逐碱基信号表阈值）：`active_frac`（活跃基因 genebody 密度下限）、`peak_frac`（位点阈值）、`noise_quantile`（噪声分位，0=关）、`min_reps`（组内重复数下限）；`min_genebody_length` 复用 `tss.min_genebody_length`。

> **信号与归一化命名约定**（bigWig/bedGraph）：`{sample}_[single|full]_[plus|minus]_[cpm|spike].{bedgraph|bigWig}`
> - `plus`/`minus` = 基因链方向信号（PRO-seq reverse 建库约定，已处理链向）。
> - `single` = 单碱基覆盖度；`full` = full-read 覆盖度。
> - `cpm` = 1e6/total_mapped（total_mapped = read1 mapped 数）；`spike` = 1e6/spike_count（开 spike 时并行多出一轨，不覆盖 cpm）。raw 无后缀（整数，供 PI 计数）。

---

## 5. 工作流步骤

0. **配置解析 / 基因组准备**：`parse_config` → `prepare_genome`（代表转录本 GTF + promoter/genebody/gene BED + genebody-union SAF + bowtie2 索引）。
1. **质控**：Fastp 去接头 + FastQC（raw/trimmed）。
2. **比对**：Bowtie2（read1 单端化供覆盖度分析）。
3. **spike-in**（条件）：统计每样本 spike reads，生成归一化因子。
4. **定量**：featureCounts 基因体计数（genebody union）。
5. **差异表达**（条件）：DESeq2（依赖 `compared_groups`）。
6. **富集**（条件）：GO/KEGG/GSEA。
7. **Pol II 覆盖度**：bedtools genomecov → bedGraph/bigWig（按 `signal_mode` 门控：single=单碱基 5' 端 / full=full-read / both=两者都算；各出 raw + cpm，开 spike 时再并行出 spike 轨）。
8. **TSS metagene**：deeptools computeMatrix + plotProfile/plotHeatmap（每样本 + 组平均；开 spike 用 spike 轨、否则 CPM，命名随信号）。
9. **暂停指数**：promoter / genebody 单碱基计数比值。
10. **信号表**（`signal_mode != full`）：逐碱基 Pol II 活性位点（活跃基因 pause 窗口，`Count`=组内 raw 末端数求和、`Signal`=组内归一化均值，归一化轨随 spike/CPM 切换）→ `pol2_signal_table.tsv`。
11. **报告打包**（条件）：Report.R 汇总（信号表落入报告 `14.Pol_II_active_site/14.1.PROSeq_profiling/`）。

---

## 6. 输出目录

| 目录 | 内容 |
|---|---|
| `01.Info/` | 参考基因组信息、代表转录本 GTF、promoter/genebody/gene/tss BED、genebody-union SAF、chrom.sizes、spike 染色体名单 |
| `02.Fastqc/` | FastQC raw/trimmed（zip + html） |
| `03.Data_QC/` | Fastp 报告、trimmed reads、碱基质量图、统计 |
| `04.Alignment/` | BAM/BAI、比对率 |
| `05.Quantification/` | 基因体计数矩阵、profile 表（`.Fpkm/.Cpm/.Rpkm/.Spike`）、spike 因子表 |
| `06.Differential_Expression/` | DESeq2 结果、CheckDE、火山图、PCA |
| `07.enrich/` | GO/KEGG/GSEA 富集结果 |
| `08.Pol2_coverage/` | 覆盖度 bedGraph/bigWig（raw/cpm/spike × single/full）、组平均 bigWig、`pol2_signal_table.tsv`（逐碱基信号表） |
| `09.TSS_Metagene/` | metagene profile/heatmap 图与矩阵 |
| `10.Pausing_Index/` | 暂停指数表、boxplot |
| `11.Report/` | 打包报告 |

---

## 7. 常见问题
