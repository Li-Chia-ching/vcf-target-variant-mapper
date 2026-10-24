# CHANGELOG — vcf-target-variant-mapper v2.0

发布日期：2026-10-24
定位升级：**VCF 靶基因变异注释工具包 → 面向多性状共有基因发现的全流程分析包**
（原 VCF 映射链完整保留，新增可选多性状分支，v1 全部功能回归通过）

---

## ✨ 新增：`modules/multitrait_gene_discovery/` 多性状共有基因发现模块

### 方法 A（第一代）— `run_magma_multitrait.py`
- **geneGPS 概念基因级多性状联合评分**：GPS-N（基因窗口内显著 SNP 计数）+ GPS-M
  （窗口内信号幅度 Σ−log10 p），跨性状 Stouffer 联合 Z + GPS 复合优先级分。
  ⚠️ geneGPS 公开检索无对应仓库，本模块按概念实现、不绑定第三方包（遵用户嘱托）。
- **MAGMA 基因级检验**：SNP-wise mean 模型（默认），基因窗口 **上游 2 kb / 下游 1 kb**
  （按链向展开）；MAGMA 可执行文件路径可在配置 `magma.binary` / 环境变量
  `MAGMA_BIN` 中指定；二进制缺失时自动回退内置 emulated 引擎并在
  `run_report.tsv` 中显式标注（Z 为 LD 未校正近似，不与真 MAGMA 混淆）。
- **截断数据模式切换（核心）**：自动检测（任一性状 max(p) ≤ 阈值×1.5 即判截断；
  华智 p<1e-5 数据命中）→ 自动切换**"显著信号基因级富集模式"**：
  - 只对携带截断信号的基因打分排序，不做全基因组基因级检验；
  - 结果定位为**优先级排序**，非严格基因级 p；
  - 每行强制标注 `截断数据模式，基因级显著性不可直接解释`。
- 输出：`.gene_scores.tsv` / `.ranked.tsv` / `.run_report.tsv` / `.truncation_check.tsv`。

### 方法 B（第二代）— `run_coloc.R`
- **主引擎 colocboost**（StatFunGen/colocboost 最新版）：sumstat 列表
  （variant/z 或 beta+sebeta/n）+ LD 矩阵，提取 `cos_details$cos$cos_index`
  共定位信号。
- **备选 HyPrColoc**（**jrs95/hyprcoloc**，按用户指定仓库；**RcppEigen 锁定 0.3.3.9.3**）。
- coloc v5/6 SuSiE 仅在用户显式要求（`method_b.susie: true`）且区间 LD 充分时补跑。
- **截断数据铁律 —— 不得静默运行**：先做完整性检查，不满足 colocboost/HyPrColoc
  统计前提时按配置选择降级路径：
  - **(a) warn**：醒目警告 + 建议获取完整区域统计数据，分析中止；
  - **(b) poemcoloc**：LD 信息充分时改用 POEMColoc（面向 top-SNP/不完整数据，s=0.5）；
  - **(c) skip**：跳过共定位，方法 A 评分作为主要证据。
  - `auto`：有充分 LD 走 (b)，否则 (a)。降级结果在 `degraded_path` 与
    `confidence_label` 列显式标注，绝不与完整数据 PP 值等同解释。
- 输出：`.coloc_results.tsv` / `.run_report.tsv`。

### 功能优先级排序 — `functional_prioritization.py`
- 同源映射（*M. truncatula* / 拟南芥，tblastn/BLAST 表格直接输入）。
- **优先保留**：GA/IAA 激素信号、叶片发育、细胞分裂/扩展、转录调控（GRAS/TCP/WOX/
  bZIP/MYB/WRKY/NAC/bHLH/homeobox/MADS）。
- **直接剔除**：转座子/逆转座子元件、未知功能蛋白——例外：同源基因有明确描述时
  可转 neutral 或按同源描述重分类（规则写死在脚本，剔除理由逐行输出）。
- 输出表**强制包含**：共定位证据（或联合评分）、同源基因 ID、功能描述、过滤理由、
  总置信度标注。

---

## 📦 结构更新
- `modules/multitrait_gene_discovery/`：`run_magma_multitrait.py`、`run_coloc.R`、
  `functional_prioritization.py`。
- `example_config.yaml`：全部可调参数（数据模式、截断检测、基因窗口、MAGMA 路径、
  降级路径、PP 阈值、功能注释输入）。
- `requirements.txt`（Python）+ `DESCRIPTION`（R 依赖清单，含 **RcppEigen 0.3.3.9.3
  锁定**与安装顺序说明）。
- `tests/`：回归测试 + 双模式冒烟测试（完整/截断），详见 README"截断数据使用说明"。

## 🔌 主流程可选分支
SKILL.md 工作流新增 **Stage M（可选）**：当用户提供多性状汇总统计（或截断信号表）时
触发；触发条件与完整性判定逻辑均可配置（`data_mode: auto`）。未提供多性状数据时，
v1 的 Stage 1–9 映射链行为完全不变。

## ⚠️ 截断数据置信度声明（合规要求）
截断模式全部输出均定位为**"基于显著信号的优先级筛选"**，置信度分级：
完整数据 MAGMA/colocboost 证据 > 截断模式 POEMColoc 降级证据 > 截断模式 GPS 评分
排序。**任何截断模式结果不得表述为"基因级全基因组显著"或与完整数据 PP 值等同**。

## 🔧 已知限制
- MAGMA 官方站点（ctg.cncr.nl）本环境下反复 TLS 断连（curl exit 56），二进制未随包
  分发；请从 <https://ctg.cncr.nl/software/magma/> 自行下载并在 `magma.binary`
  中指定路径。未指定时 emulated 引擎仅供流程验证/排序参考。
- hyprcoloc 编译要求 RcppEigen **0.3.3.9.3**（已锁定；与 coloc/susieR/colocboost/
  POEMColoc 无冲突，经依赖图核查）。
- 方法 B 格式 2（逐性状散表）需先聚类为区间面板（参照
  `02_analysis/08_GWAS/coloc_stage1/build_coloc_inputs.py` 思路），模块直接消费
  面板目录，避免重复实现聚类。

## ✅ 回归与冒烟
- v1 全部 9 个脚本回归通过（参数解析 + 最小输入冒烟）。
- 模拟双性状 summary stats 冒烟：
  - 完整模式：MAGMA（emulated 引擎）+ colocboost 正常执行；
  - 截断模式：正确触发降级路径（a/b/c）并输出清晰警告，结果带截断标注。
