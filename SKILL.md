---
name: vcf-target-variant-mapper
description: 从全基因组重测序 VCF 中检索目标基因（如 DELLA/GRAS 转录因子家族）的变异并完成"变异→蛋白→结构域→三维结构"全链路注释与可视化；并支持多性状汇总统计的共有基因发现（geneGPS 概念+MAGMA 基因级联合评分、colocboost/HyPrColoc 贝叶斯共定位、功能优先级排序，含截断数据模式切换）。This skill should be used when the user wants to find a named gene's variants in a resequencing VCF, annotate them with a custom SnpEff database, filter HIGH/MODERATE effects, map missense changes onto the reference protein and domain architecture (e.g. DELLA motif, TVHYNP, GRAS LHRI/VHIID/LHRII/PFY/SAW), compare mutant vs reference protein sequences, compute conservation (Shannon entropy) and solvent accessibility, inject scores into AlphaFold structures as B-factors, and render publication-ready PyMOL figures; or discover pleiotropic genes shared across multiple traits from GWAS summary statistics (MAGMA gene-level multi-trait scoring, Bayesian colocalization with colocboost/HyPrColoc/POEMColoc, functional prioritization with truncation-aware mode switching). Triggers: "VCF 找基因"、"靶基因变异"、"DELLA 变异"、"SnpEff 注释"、"变异映射结构域"、"突变体结构可视化"、"variant mapping"、"missense on protein structure"、"多性状共有基因"、"共定位"、"colocalization"、"MAGMA 基因分析"、"截断数据"、"pleiotropic genes"。
agent_created: true
---

# VCF Target Variant Mapper — VCF 靶基因变异注释与结构可视化

## Overview

把重测序 VCF 中**指定目标基因**（示例：紫花苜蓿 DELLA 家族 Msa085132）的变异，沿着
「变异检出 → 效应注释 → 蛋白映射 → 结构域定位 → 三维结构 → 保守性/溶剂可及性 → 出版级渲染」
九个阶段推进到最终成图。核心产出：目标基因变异 CSV、效应分布图、蛋白差异报告、结构域标注图、
PyMOL 保守性着色渲染图 + 图注。

适用输入：已完成 SnpEff 注释的 `.ann.vcf`（或可自行建库注释的 raw VCF + GFF3 + 基因组 FASTA）。

## Workflow Decision Tree

```
有多性状汇总统计 / 截断信号表（供共有基因发现）？
 ├─ 有 → 【可选分支】Stage M1–M3 多性状共有基因发现（截断数据自动模式切换）
 └─ 无 → 跳过 Stage M
有 SnpEff 注释的 VCF？
 ├─ 有 → Stage 2 直接提取目标基因变异
 └─ 无 → Stage 1 先建自定义 SnpEff 数据库并注释
之后顺序执行 Stage 3–9；仅需表格结果可止于 Stage 4；仅需结构域图可止于 Stage 6。
```

**Stage M 触发与完整性判定（可配置）**：提供多性状 summary stats（或 p<1e-5 截断
信号表）即触发；`example_config.yaml` 的 `data_mode: auto` 自动检测完整性——
任一性状 max(p) ≤ p_trunc×1.5 判为截断并切换降级路径（详见 Stage M 与 README
"截断数据使用说明"）。未提供多性状数据时 v1 映射链行为完全不变。

## Stage 1 — 构建 SnpEff 自定义数据库（无注释 VCF 时）

1. 建目录 `data/<genome>/`，放入 `sequences.fa`（基因组）与 `genes.gff`（GFF3 注释）。
2. 在 `snpEff.config` 追加一行：`<genome>.genome : <genome>`。
3. 构建：`snpEff build -gff3 -v <genome>`
4. 注释：`snpEff ann <genome> input.vcf > annotated.vcf`
5. 快速验证：`bcftools view -i 'INFO/ANN ~ "Msa085132"' annotated.vcf | head`

细节与常见问题见 `references/snpeff_custom_database.md`。

## Stage 2 — 提取目标基因变异

```bash
python scripts/extract_target_gene_variants.py \
  --vcf annotated.vcf \
  --gene Msa085132 \
  --impact high-moderate \
  --out-csv output/Msa085132_variants.csv
```

- 用 cyvcf2 逐条解析 `INFO/ANN`；ANN 为逗号分隔多条记录、竖线分隔字段。
- **ANN 字段索引（0-based）**：1=EFFECT，2=IMPACT，3=GENE，10=HGVS.p（见 Key Pitfalls）。
- `--impact high-moderate` 只保留 HIGH/MODERATE；需要同义变异等全部记录用 `--impact all`。
- 同时输出各效应类型计数汇总，供 Stage 3 画图。

## Stage 3 — 效应分布图（ggplot2）

```bash
Rscript scripts/plot_effect_distribution.R output/Msa085132_variants.csv output/effect_distribution.pdf
```

按 EFFECT 类型（missense_variant、synonymous_variant、intron_variant…）计数作条形图，
Nature 风格（despine、无网格、Arial、600 dpi PDF）。

## Stage 4 — 错义变异映射到参考蛋白

```bash
python scripts/map_missense_to_protein.py \
  --csv output/Msa085132_variants.csv \
  --ref-fasta data/protein/Msa085132_protein.faa \
  --out output/Msa085132_protein_map.tsv
```

- HGVS.p 三字母格式正则：`p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})`
- 内置 `AA_3_TO_1` 三字母→单字母映射；**必须**与参考蛋白 FASTA 校验（位点氨基酸一致才接受）。
- 输出：位点编号、参考/突变氨基酸、是否落在已知基序内（先验知识可选填）。

## Stage 5 — 突变体 vs 参考蛋白序列差异报告

```bash
python scripts/align_and_diff.py \
  --ref data/protein/Msa085132_protein.faa \
  --mut data/protein/Msa085132_mutant.faa \
  --out output/align_diff.txt
```

Biopython `PairwiseAligner` 全局比对，输出逐位差异（含 InDel），并检查是否非 3 倍数移码。

## Stage 6 — 变异映射到结构域

```bash
python scripts/map_variants_to_domains.py \
  --map output/Msa085132_protein_map.tsv \
  --out output/domain_map.pdf
```

- DELLA 家族先验结构域（紫花苜蓿参考编号，可按实际蛋白调整）：
  N 端可变区、DELLA 基序 56–70、TVHYNP 基序 81–95；
  GRAS 区：LHRI、VHIID、LHRII、PFY、SAW（详见 `references/della_domain_architecture.md`）。
- 输出线性结构域图 + 变异位点标记。

## Stage 7 — 三维结构保守性与溶剂可及性

```bash
python scripts/structure_conservation.py \
  --pdb data/structure/Msa085132.pdb \
  --msa data/msa/DELLA_family.aa.msa \
  --out output/entropy_asa.tsv
```

- ASA：Shrake-Rupley 算法，**< 20 Å² 判为埋藏位点**。
- Shannon 熵：由多序列比对逐位计算，**< 0.2 严格保守，< 1.0 中度保守**。
- PDB 残基编号与蛋白序列编号必须先做偏移校验（AlphaFold PDB 常从 1 开始但需实测）。

## Stage 8 — 熵值注入 B-factor + PyMOL 渲染

```bash
python scripts/entropy_bfactor.py --pdb data/structure/Msa085132.pdb \
  --entropy output/entropy_asa.tsv --out output/Msa085132_entropy.pdb
python scripts/generate_pymol_scripts.py --pdb output/Msa085132_entropy.pdb \
  --out output/render_conservation.pml
pymol -cq output/render_conservation.pml
```

- **熵值写入 B-factor 会覆盖 AlphaFold 的 pLDDT**——仅对副本操作，原始 PDB 永不改动。
- 渲染参数：`spectrum b, blue_white_red`，`ray 3000,2250`，白色背景，输出 PNG。

## Stage 9 — 图例与拼版、图注

```bash
python scripts/legend_and_panel.py --png output/render.png \
  --out output/figure1_panel.png --caption output/figure1_caption.txt
```

matplotlib `bwr_r` 色条图例 + PIL 拼版；图注模板含结构域标注、保守性阈值与渲染参数说明。

## Stage M —（可选分支）多性状共有基因发现

当用户提供多性状 GWAS 汇总统计（完整或截断）时启用；三个脚本位于
`modules/multitrait_gene_discovery/`，配置统一走 `example_config.yaml`。

### M1 — 方法 A：geneGPS 概念 + MAGMA 基因级多性状联合

```bash
python modules/multitrait_gene_discovery/run_magma_multitrait.py \
  --sumstats data/ph.tsv data/mf.tsv data/sh.tsv \
  --trait-names PH MF SH \
  --genes data/genes.tsv \
  --out-prefix output/multitrait/gene_scores \
  --config example_config.yaml
```

- SNP 信号按基因窗口（**上游 2 kb / 下游 1 kb**，链向感知）聚合为 GPS-N（显著 SNP
  计数）/ GPS-M（信号幅度）+ MAGMA SNP-wise mean 基因级 Z，跨性状 Stouffer 联合。
  ⚠️ geneGPS 公开检索无对应仓库——按概念实现，不绑定第三方包；严格基因级检验由
  MAGMA 承担（`magma.binary` 配置路径；缺失时内置 emulated 引擎兜底并显式标注）。
- **截断数据（如 p<1e-5）自动切换"显著信号基因级富集模式"**：仅对携带信号的基因
  排序，结果为优先级而非基因级 p，每行强制标注"截断数据模式，基因级显著性不可直接解释"。

### M2 — 方法 B：贝叶斯共定位确认

```bash
Rscript modules/multitrait_gene_discovery/run_coloc.R \
  --panels-dir output/multitrait/panels --ld-dir output/multitrait/ld \
  --traits PH,MF,SH --out-prefix output/multitrait/coloc \
  --config example_config.yaml
```

- 主引擎 **colocboost**（StatFunGen/colocboost），备选 **HyPrColoc**
  （jrs95/hyprcoloc，RcppEigen 锁定 0.3.3.9.3）；coloc SuSiE 仅显式开启且 LD 充分时补跑。
- **截断数据不得静默运行**：完整性检查不过 → 降级路径 (a) 警告+建议完整区域统计 /
  (b) LD 充分改用 POEMColoc / (c) 跳过 coloc 用方法 A 证据；降级结果逐行标注。

### M3 — 功能优先级排序

```bash
python modules/multitrait_gene_discovery/functional_prioritization.py \
  --genes-a output/multitrait/gene_scores.ranked.tsv \
  --genes-b output/multitrait/coloc.coloc_results.tsv \
  --homologs data/mt_at_homologs.tsv \
  --gene-desc data/gene_desc.tsv \
  --out output/multitrait/prioritized_candidates.tsv
```

- 同源映射 *M. truncatula*/拟南芥；**优先保留** GA/IAA 激素、叶片发育、细胞分裂、
  转录调控（GRAS/TCP/WOX）；**直接剔除**转座子与未知功能（同源有明确描述除外）。
- 输出强制包含：共定位证据（或联合评分）、同源 ID、功能描述、过滤理由。

## Key Pitfalls

1. **ANN 字段布局**：逗号分隔多条注释、竖线分隔字段；**0-based 索引 1=EFFECT、2=IMPACT、3=GENE、10=HGVS.p**。解析前先打印一条真实 ANN 验证布局。
2. **HGVS.p 三字母格式**：必须用 `p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})` 校验，并与参考蛋白 FASTA 交叉验证位点氨基酸；不一致即弃用。
3. **PDB 编号偏移**：结构残基号 ≠ 蛋白序列号时先求偏移量再映射，否则结构域全错位。
4. **B-factor 覆盖 pLDDT**：熵值注入只操作副本文件。
5. **VCF ≠ GFF**：VCF 是变异坐标表，不能当注释用，也不能直接 BLAST 跨物种比对（跨物种用 MUMmer/LASTZ/Minimap2）。
6. **坐标大小写**：GWAS/VCF 常为 `chr4_1`，GFF3 为 `Chr4_1`，匹配前统一 `.lower()`。
7. **重测序不交付自定义 GFF**：注释以参考基因组 GFF3 为准；自定义 SnpEff 库见 `references/snpeff_custom_database.md`。
8. **截断汇总统计（p<1e-5）≠ 完整数据**：截断模式产出一律为"基于显著信号的优先级筛选"，禁止表述为基因级全基因组显著、禁止与完整数据 PP.H4 等同解释；置信度分级 完整 MAGMA/colocboost > 截断 POEMColoc 降级 > 截断 GPS 排序。
9. **MAGMA 双引擎**：`magma.engine=auto` 时二进制缺失回退 emulated 近似（LD 未校正），`run_report.tsv` 有记录——引用 Z 值前必查引擎标注；正式分析请提供真实 MAGMA + `--bfile` 参考基因型。

## Resources

### scripts/
九个可执行脚本：`extract_target_gene_variants.py`、`plot_effect_distribution.R`、
`map_missense_to_protein.py`、`align_and_diff.py`、`map_variants_to_domains.py`、
`structure_conservation.py`、`entropy_bfactor.py`、`generate_pymol_scripts.py`、`legend_and_panel.py`。

### modules/multitrait_gene_discovery/
`run_magma_multitrait.py`（方法 A：GPS-N/GPS-M + MAGMA + 截断切换）、
`run_coloc.R`（方法 B：colocboost/HyPrColoc/POEMColoc + 截断降级 a/b/c）、
`functional_prioritization.py`（同源映射 + 功能过滤 + 证据整合）。

### 配置与依赖
`example_config.yaml`（全部参数）、`requirements.txt`（Python）、`DESCRIPTION`
（R 依赖清单，RcppEigen 锁定 0.3.3.9.3）、`tests/`（回归 + 双模式冒烟）。

### references/
- `snpeff_custom_database.md` — SnpEff 自定义库构建详解
- `della_domain_architecture.md` — DELLA/GRAS 结构域坐标与基序
- `faq_vcf_concepts.md` — VCF 概念 FAQ（VCF 与 GFF、BLAST、找基因等）
