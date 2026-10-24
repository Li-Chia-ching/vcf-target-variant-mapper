#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stage 2 — 从 SnpEff 注释 VCF 中提取目标基因变异并导出 CSV。

用法:
    python extract_target_gene_variants.py --vcf annotated.vcf \
        --gene Msa085132 --impact high-moderate --out-csv out.csv

ANN 字段布局（0-based，逗号分隔多条、竖线分隔字段）:
    0=Annotation  1=EFFECT  2=IMPACT  3=GENE  4=GeneID  5=FeatureType
    6=FeatureID  7=Biotype  8=Rank  9=HGVS.c  10=HGVS.p
"""
import argparse
import csv
import re
import sys
from collections import Counter

try:
    from cyvcf2 import VCF
except ImportError:
    sys.exit("[ERROR] 缺少 cyvcf2，请先: pip install cyvcf2")

IMPACT_ORDER = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}


def parse_ann(info_ann):
    """拆解 INFO/ANN 字段 -> list[dict]，非法条目跳过。"""
    records = []
    for entry in (info_ann or "").split(","):
        fields = entry.split("|")
        if len(fields) < 11:
            continue
        records.append({
            "effect": fields[1],
            "impact": fields[2],
            "gene": fields[3],
            "gene_id": fields[4],
            "feature": fields[6],
            "biotype": fields[7],
            "hgvs_c": fields[9],
            "hgvs_p": fields[10],
        })
    return records


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", required=True, help="SnpEff 注释后的 VCF")
    ap.add_argument("--gene", required=True, help="目标基因名/ID，如 Msa085132")
    ap.add_argument("--impact", default="high-moderate",
                    choices=["high-moderate", "all"],
                    help="high-moderate=仅 HIGH/MODERATE；all=全部")
    ap.add_argument("--out-csv", required=True, help="输出 CSV 路径")
    ap.add_argument("--allow-partial", action="store_true",
                    help="基因名允许部分匹配（默认精确匹配）")
    args = ap.parse_args()

    vcf = VCF(args.vcf)
    rows, eff_counter = [], Counter()
    gene_re = re.compile(re.escape(args.gene)) if args.allow_partial else None

    for var in vcf:
        ann_raw = var.INFO.get("ANN")
        if not ann_raw:
            continue
        anns = parse_ann(str(ann_raw))
        if args.allow_partial:
            hits = [a for a in anns if a["gene"] and gene_re.search(a["gene"])]
        else:
            hits = [a for a in anns if a["gene"] == args.gene]
        if not hits:
            continue
        for a in hits:
            if args.impact == "high-moderate" and a["impact"] not in ("HIGH", "MODERATE"):
                continue
            eff_counter[a["effect"]] += 1
            rows.append({
                "CHROM": var.CHROM,
                "POS": var.POS,
                "REF": var.REF,
                "ALT": ",".join(var.ALT),
                "QUAL": var.QUAL,
                "IMPACT": a["impact"],
                "EFFECT": a["effect"],
                "GENE": a["gene"],
                "FEATURE": a["feature"],
                "BIOTYPE": a["biotype"],
                "HGVS_c": a["hgvs_c"],
                "HGVS_p": a["hgvs_p"],
            })

    if not rows:
        print(f"[WARN] 未找到基因 {args.gene} 的变异（--impact={args.impact}）。"
              f"可尝试 --impact all 或 --allow-partial。")

    fieldnames = ["CHROM", "POS", "REF", "ALT", "QUAL", "IMPACT", "EFFECT",
                  "GENE", "FEATURE", "BIOTYPE", "HGVS_c", "HGVS_p"]
    with open(args.out_csv, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print(f"[OK] 共 {len(rows)} 条记录 -> {args.out_csv}")
    for eff, n in eff_counter.most_common():
        print(f"    {eff}: {n}")


if __name__ == "__main__":
    main()
