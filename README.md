# vcf-target-variant-mapper

从全基因组重测序 VCF 中定位**目标基因（如 DELLA/GRAS 家族）的功能变异**，并一路产出蛋白水平与 3D 结构水平的论文级可视化图件；**v2.0 起**支持多性状汇总统计的**共有基因发现**（geneGPS 概念 + MAGMA 基因级联合、colocboost/HyPrColoc 贝叶斯共定位、功能优先级排序，含截断数据模式切换）。

**English**: Map functional variants of a named target gene from resequencing VCF to publication-ready protein and 3D-structure conservation figures; discover pleiotropic genes shared across multiple traits from GWAS summary statistics with truncation-aware mode switching.

## 功能总览

```
VCF ──► SnpEff 自定义库注释 ──► cyvcf2 靶基因提取（HIGH/MODERATE 过滤）
      ──► 效应分布图（ggplot2）──► 错义变异映射到参考蛋白（HGVS.p 校验）
      ──► 双序列比对差异报告 ──► 结构域/基序定位（DELLA·TVHYNP·GRAS 五亚域）
      ──► ASA 埋藏/暴露 + Shannon 熵保守性 ──► 熵写入 PDB B-factor
      ──► PyMOL 谱色渲染（300 DPI）──► 图例 + 图版拼版 + 图注

多性状（可选分支 Stage M）──► 截断检测/模式切换
      ──► 方法 A：GPS-N/GPS-M + MAGMA 基因级多性状联合（窗口 2kb/1kb）
      ──► 方法 B：colocboost 主 / HyPrColoc 备 / POEMColoc 降级（截断铁律）
      ──► 同源映射（M. truncatula/拟南芥）──► 功能优先级排序（GA/IAA·叶片·细胞分裂·GRAS/TCP/WOX）
      ──► 高置信度多效性候选基因清单（证据 + 同源 + 描述 + 过滤理由）
```

## 目录结构

```
vcf-target-variant-mapper/
├── SKILL.md                 # 工作流入口（AI Agent 技能定义，含 Stage M 可选分支）
├── README.md
├── CHANGELOG_v2.0.md        # v2.0 变更记录
├── example_config.yaml      # 多性状模块配置（模式/截断/窗口/MAGMA/降级路径）
├── requirements.txt         # Python 依赖
├── DESCRIPTION              # R 依赖清单（RcppEigen 锁定 0.3.3.9.3）
├── scripts/                 # v1：VCF→蛋白→结构域→3D 九阶段（9 个脚本，v1.0 不变）
├── modules/multitrait_gene_discovery/
│   ├── run_magma_multitrait.py      # 方法 A：GPS-N/GPS-M + MAGMA + 截断模式切换
│   ├── run_coloc.R                  # 方法 B：colocboost/HyPrColoc/POEMColoc + 截断降级
│   └── functional_prioritization.py # 同源映射 + 功能过滤 + 证据整合
├── tests/                   # 回归测试 + 完整/截断双模式冒烟测试
├── tools/                   # 外部二进制放置处（MAGMA 自行下载后填入配置）
└── references/
    ├── snpeff_custom_database.md    # SnpEff 自定义库构建五步
    ├── della_domain_architecture.md # DELLA/GRAS 结构域知识 + 图注模板
    └── faq_vcf_concepts.md          # VCF 概念 FAQ（GFF 转换/跨物种 BLAST 等）
```

## 快速开始

```bash
# 1. 注释（已有注释 VCF 可跳过）
java -jar snpEff.jar ann -v msativa_v2 input.vcf > annotated.vcf

# 2. 提取靶基因高/中影响变异
python scripts/extract_target_gene_variants.py \
    --vcf annotated.vcf --gene DELLA --impact high-moderate \
    --out-csv della_high_moderate.csv

# 3. 映射到参考蛋白并生成突变序列
python scripts/map_missense_to_protein.py \
    --csv della_high_moderate.csv --ref-fasta della_reference.faa \
    --out-fasta della_mutated.faa

# 4. 差异位点报告
python scripts/align_and_diff.py \
    --ref-fasta della_reference.faa --mut-fasta della_mutated.faa \
    --out-csv protein_differences_report.csv

# 5. 结构域定位 + 结构/保守性分析
python scripts/map_variants_to_domains.py --diff-csv protein_differences_report.csv
python scripts/structure_conservation.py --diff-csv protein_differences_report.csv \
    --pdb della_model.pdb --msa della_homologs_aln.fasta

# 6. 保守性着色渲染 + 图版
python scripts/entropy_bfactor.py --pdb della_model.pdb --msa della_homologs_aln.fasta
python scripts/generate_pymol_scripts.py --diff-csv protein_differences_report.csv
pymol render_conservation.pml
python scripts/legend_and_panel.py --structure-img della_conservation_highres.png --vmin 0 --vmax 1.5
```

## 多性状共有基因发现（Stage M）快速开始

```bash
# M1 方法 A：GPS-N/GPS-M + MAGMA 基因级多性状联合（自动截断检测）
python modules/multitrait_gene_discovery/run_magma_multitrait.py \
    --sumstats ph.tsv mf.tsv sh.tsv --trait-names PH MF SH \
    --genes genes.tsv --out-prefix out/gene_scores --config example_config.yaml

# M2 方法 B：贝叶斯共定位（colocboost 主 / HyPrColoc 备 / 截断降级 a/b/c）
Rscript modules/multitrait_gene_discovery/run_coloc.R \
    --panels-dir out/panels --ld-dir out/ld --traits PH,MF,SH \
    --out-prefix out/coloc --config example_config.yaml

# M3 功能优先级排序（同源映射 + 功能过滤 + 证据整合）
python modules/multitrait_gene_discovery/functional_prioritization.py \
    --genes-a out/gene_scores.ranked.tsv --genes-b out/coloc.coloc_results.tsv \
    --homologs mt_at_homologs.tsv --gene-desc gene_desc.tsv \
    --out out/prioritized_candidates.tsv
```

## 截断数据使用说明

数据方（如华智生物）常只交付 **p < 1e-5 的显著信号位点**（截断汇总统计）。此类数据
丢失全部非显著位点信息，**统计性质与完整汇总统计根本不同**，本包处理规则如下：

### 1. 自动检测与模式切换

- `data_mode: auto` 时逐性状检测：**max(p) ≤ p_trunc × 1.5（默认阈值 1e-5）→ 判为截断**；
  任一性状截断即整体按截断处理。检测证据写入 `<out-prefix>.truncation_check.tsv`。
- 也可在配置/命令行强制指定 `data_mode: truncated | complete`。

### 2. 方法 A（geneGPS + MAGMA）截断行为

- 自动切换**"显著信号基因级富集模式"**：仅对携带截断信号的基因做 GPS-N/GPS-M
  评分排序，**不运行**全基因组基因级检验（无非显著位点时检验不可解释）。
- 结果定位为**优先级排序**，每行强制携带标注：
  **"截断数据模式，基因级显著性不可直接解释"**。

### 3. 方法 B（colocboost/HyPrColoc）截断铁律 —— 不得静默运行

截断数据不满足贝叶斯共定位的统计前提（需要完整 p 值分布做先验校准），模块拒绝
静默输出 PP 值，按 `method_b.truncation.action` 选择降级路径：

| 路径              | 触发条件                                | 行为                                                 |
| :---------------- | :-------------------------------------- | :--------------------------------------------------- |
| **(a) warn**      | 默认（或 LD 不充分）                    | 明确警告 + 建议获取完整区域统计数据，分析中止        |
| **(b) poemcoloc** | LD 信息充分（区间 LD 矩阵覆盖全部 SNP） | 改用 **POEMColoc**（面向 top-SNP/不完整数据，s=0.5） |
| **(c) skip**      | 配置 `action: skip`                     | 跳过共定位，方法 A 评分作为主要证据                  |
| auto              | —                                       | 有充分 LD 走 (b)，否则 (a)                           |

降级结果在 `degraded_path` 与 `confidence_label` 列逐行显式标注。

### 4. 置信度分级（引用时必须遵守）

```
完整数据 MAGMA / colocboost 证据  >  截断模式 POEMColoc 降级证据  >  截断模式 GPS 评分排序
```

⚠️ **任何截断模式结果不得表述为"基因级全基因组显著"，不得与完整数据 PP.H4 等同
解释**；论文/报告中应表述为"基于显著信号的优先级筛选"。如需严格统计结论，请向数据
方索取**完整（未截断）区域级汇总统计**（含非显著位点）后重跑标准流程。

## 依赖

- Java 11+ · SnpEff · bcftools
- Python 3.9+：`pip install -r requirements.txt`（cyvcf2 biopython numpy scipy pandas pyyaml matplotlib pillow）
- R 4.2+：`yaml`、`colocboost`（StatFunGen/colocboost）、`hyprcoloc`（**jrs95/hyprcoloc**，
  **RcppEigen 锁定 0.3.3.9.3**）、`POEMColoc`、`coloc`/`susieR`（可选）——安装顺序见 `DESCRIPTION`
- MAGMA v1.10（可选但推荐）：从 <https://ctg.cncr.nl/software/magma/> 下载，
  在 `example_config.yaml` 的 `magma.binary` 或环境变量 `MAGMA_BIN` 指定路径；
  缺失时内置 emulated 引擎兜底（LD 未校正近似，仅用于流程验证/排序参考）
- PyMOL 2.x（渲染）

## 使用须知

- 所有脚本使用相对路径/命令行参数，不硬编码绝对路径。
- HGVS.p 为三字母氨基酸码，脚本会对照参考蛋白 FASTA 校验坐标（防止转录本/蛋白异构体错配）。
- 熵写入 B-factor 会覆盖 AlphaFold pLDDT，务必保留原始 PDB 副本。
- VCF 不能直接转 GFF、不能用于跨物种 BLAST——详见 `references/faq_vcf_concepts.md`。
- 截断汇总统计的置信度标注规则见上文"截断数据使用说明"，引用结果前必读。

## License

MIT
