# SnpEff 自定义数据库构建详解

适用：参考基因组不在 SnpEff 官方库里（如自建的紫花苜蓿 ZM4_V2 四倍体基因组），需要对重测序 VCF 做基因级效应注释。

## 目录结构

```
snpEff/
├── snpEff.config
└── data/
    └── <genome>/
        ├── sequences.fa   # 基因组（染色体名必须与 VCF 的 CHROM 完全一致）
        └── genes.gff      # GFF3 注释
```

## 五步流程

1. **准备数据**：把基因组 FASTA 与 GFF3 放入 `data/<genome>/`。文件名必须是 `sequences.fa` 与 `genes.gff`。
2. **注册基因组**：在 `snpEff.config` 中追加一行：
   ```
   <genome>.genome : <genome>
   ```
3. **构建数据库**：
   ```bash
   snpEff build -gff3 -v <genome>
   ```
4. **注释 VCF**：
   ```bash
   snpEff ann <genome> input.vcf > annotated.vcf
   ```
5. **验证**：
   ```bash
   bcftools view -i 'INFO/ANN ~ "Msa085132"' annotated.vcf | head
   ```

## 常见问题

- **注释为空**：99% 是染色体命名不一致。VCF 用 `chr4_1`、GFF3 用 `Chr4_1` 时 SnpEff 匹配不到。构建前统一大小写（推荐都转小写）。
- **GFF3 须含 `ID`/`Parent` 层级**：gene → mRNA → CDS/exon 链完整，否则 SnpEff 只能注释到基因级。
- **四倍体/多拷贝基因组**：每个亚基因组染色体（如 chr4_1..chr4_4）独立命名即可，SnpEff 不感知倍性，坐标正确即可。
- **重测序公司不交付自定义 GFF**：以参考基因组 GFF3 为准；VCF 是变异坐标表，不能当注释使用。
- **只想快速找基因而不建库**：若 VCF 已带 `INFO/ANN`（如华智生物交付的 `.ann.vcf`），可跳过 Stage 1 直接 Stage 2；或用 bcftools 按坐标区间筛：
  ```bash
  bcftools view -r chr4_1:77259347-77261462 input.vcf
  ```

## ANN 字段速查（0-based）

| 索引 | 字段 | 示例 |
|---|---|---|
| 0 | Annotation (EFFECT) | missense_variant |
| 1 | EFFECT | missense_variant |
| 2 | IMPACT | MODERATE |
| 3 | GENE | Msa085132 |
| 4 | GeneID | Msa085132 |
| 5 | FeatureType | transcript |
| 6 | FeatureID | Msa085132.1 |
| 7 | Biotype | protein_coding |
| 8 | Rank | 1/2 |
| 9 | HGVS.c | c.383G>A |
| 10 | HGVS.p | p.Gly128Glu |

> 多条注释以逗号分隔（一个变异影响多个基因/转录本）；字段内以竖线 `|` 分隔。解析前先打印一条真实 ANN 验证布局。
