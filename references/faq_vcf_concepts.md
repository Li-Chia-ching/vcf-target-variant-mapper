# VCF 核心概念 FAQ

## Q1. 能用 VCF 替代 GFF 做注释吗？

**不能。** VCF 是"变异坐标表"（哪条染色体、哪个位置、什么碱基变成了什么），GFF 是"基因结构注释表"（基因/转录本/外显子的坐标与层级）。两者信息维度不同。变异要获得基因级效应注释，必须走 SnpEff/VEP 等注释工具 + GFF3。

## Q2. 能直接拿 VCF 去 BLAST 跨物种比对吗？

**不能。** VCF 只含差异位点，不含完整序列，BLAST 无从下手。跨物种共线性/序列比对请用 MUMmer（nucmer）、LASTZ 或 Minimap2，输入为基因组 FASTA。

## Q3. 怎么在 VCF 里找一个已知基因的变异？

三种方法（按优先级）：
1. **坐标区间筛**（最快）：用 GFF3 查到基因起止坐标后
   ```bash
   bcftools view -r chr4_1:77259347-77261462 input.vcf
   ```
2. **注释字段筛**（需已注释 VCF）：
   ```bash
   bcftools view -i 'INFO/ANN ~ "Msa085132"' annotated.vcf
   ```
3. **脚本解析**：本 Skill 的 `scripts/extract_target_gene_variants.py`（支持 IMPACT 过滤与 HGVS 提取）。

## Q4. 重测序公司会给自定义 GFF 吗？

一般**不给**。交付物是 VCF、群体统计表、GWAS 信号表和以参考基因组为准的注释（SnpEff/ANNOVAR）。如果参考基因组（如 ZM4_V2）的 GFF3 由你自己持有，注释一致性以你手里的 GFF3 为准；发现注释基因 ID 对不上时，先核对参考版本号。

## Q5. VCF 里的 ANN 字段怎么读？

`INFO/ANN` 由 SnpEff 生成：**逗号**分隔多条注释（一个变异可影响多个基因/转录本），每条内部用**竖线**分隔字段。常用字段（0-based 索引）：1=EFFECT、2=IMPACT、3=GENE、9=HGVS.c、10=HGVS.p。完整表见 `snpeff_custom_database.md`。

## Q6. HIGH/MODERATE/LOW/MODIFIER 怎么取舍？

- **HIGH**：终止获得/丢失、移码、剪接位点破坏 —— 必看。
- **MODERATE**：错义变异、部分 in-frame InDel —— 功能研究主力层。
- **LOW**：同义变异 —— 仅在做密码子使用/剪接预测时需要。
- **MODIFIER**：内含子、UTR、基因间区 —— 群体遗传分析（如 GWAS 信号定位）时需要，基因功能分析常过滤掉。

## Q7. MAF 过滤会丢掉什么？

以 212 样本四倍体为例：848 份等位拷贝，`--maf 0.05` 要求 MAC≥43。**私有/稀有变异（单样本特异，MAF≤0.47%）100% 被剔除**；加上 `--max-missing 0.9`、双等位过滤后，保留率可能仅 ~7%（48.4M → 3.6M）。这是行业常规，但意味着**克隆验证的稀有功能等位变异可能根本不在交付 VCF 里**——此时必须向公司索取未过滤 VCF/gVCF 才能定论。
