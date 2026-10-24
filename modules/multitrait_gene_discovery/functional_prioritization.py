#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
functional_prioritization.py — 功能注释与优先级排序（方法 B 第二步 / 收敛步）
====================================================================

【职责】
  对方法 A（geneGPS+MAGMA 联合评分）与方法 B（colocboost/HyPrColoc/POEMColoc）
  产出的候选基因做同源映射 + 功能分类 + 优先级裁定，输出高置信度多效性候选清单。

【分类规则（与用户规则一致）】
  优先保留（priority）：植物激素信号转导（GA/IAA）、叶片发育、细胞分裂、
    转录调控（GRAS/TCP/WOX 等）；
  直接剔除（force-drop）：转座子/逆转座子元件、未知功能蛋白
    ——例外：若同源基因有明确功能描述（hit_desc 信息量足够），未知注释可转为
    neutral/priority 按同源描述再分类；
  其余：neutral（保留但不优先）。

【输出（每行强制包含）】
  共定位证据（PP.H4/cos 或联合评分）、同源基因 ID、功能描述、过滤理由。

【用法】
  python functional_prioritization.py \
      --genes-a output/multitrait/gene_scores.ranked.tsv \
      [--genes-b output/multitrait/coloc.coloc_results.tsv] \
      --homologs data/mt_at_homologs.tsv \
      [--gene-desc data/gene_desc.tsv] \
      --out output/multitrait/prioritized_candidates.tsv \
      [--config example_config.yaml]

  --homologs 表列：gene_id, hit_id, hit_desc, pident, evalue, [source]
  --gene-desc 表列：gene_id, description, [pfams]（GFF3/InterPro 注释，可选）

依赖: pandas, numpy
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from typing import Dict, List, Optional, Sequence

import pandas as pd

# ----------------------------------------------------------------------------
# 功能分类正则（先验知识；可被配置覆盖）
# ----------------------------------------------------------------------------
TE_PAT = re.compile(
    r"transpos|retrotranspos|retroelement|\bRVT_1\b|\bRVT_2\b|\brve\b|\bRVP_2\b"
    r"|Retrotrans_gag|gag_pre-integrs|zf-CCHC|DUF4283|integrase|\bTy1\b|\bcopia\b"
    r"|gypsy|\bHelitron\b|\bMutator\b|\bHarbinger\b|\bPIF\b.*transpos|逆转座|转座",
    re.I)

UNKNOWN_PAT = re.compile(
    r"^(no annotation|unknown( protein)?|uncharacterized( protein)?"
    r"|hypothetical( protein)?|expressed protein|unnamed protein|\s*)$", re.I)

PRIORITY_PATTERNS = {
    "激素信号(GA/IAA)": re.compile(
        r"DELLA|GID1|gibberell|auxin|\bIAA\b|\bSAUR\b|\bARF\b|Aux/IAA|indole-3-acetic"
        r"|赤霉素|生长素", re.I),
    "叶片发育": re.compile(
        r"\bleaf\b|leafy|\bKNOX\b|knotted|\bTCP\b|\bWOX\b|\bNAM\b|\bCUC\b|\bYABBY\b"
        r"|blade|leaflet|叶片|叶发育", re.I),
    "细胞分裂/扩展": re.compile(
        r"cyclin|\bCDK\b|kinesin|expansin|\bXET\b|xyloglucan|cell division|cell cycle"
        r"|cytokinesis|细胞分裂|细胞扩展", re.I),
    "转录调控(GRAS/TCP/WOX等)": re.compile(
        r"\bGRAS\b|\bTCP\b|\bWOX\b|\bbZIP\b|\bMYB\b|\bWRKY\b|\bNAC\b|\bbHLH\b"
        r"|homeobox|\bMADS\b|transcription factor|转录因子", re.I),
}

NEUTRAL_DESC_OK = re.compile(r"[A-Za-z]{4,}", re.I)  # 同源描述信息量粗判


# ----------------------------------------------------------------------------
# 载入
# ----------------------------------------------------------------------------
def _read(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep=None, engine="python")


def classify(desc: str, pfams: str, hit_desc: str) -> Dict[str, str]:
    """返回 功能类别 / 决定 / 决定依据。"""
    desc_s = desc or ""
    pfam_s = pfams or ""
    hit_s = hit_desc or ""
    text_own = f"{desc_s} {pfam_s}"
    text_all = f"{text_own} {hit_s}"

    # 1. 转座子：自身或同源命中均为元件 → 剔除
    if TE_PAT.search(text_own) and (TE_PAT.search(hit_s) or not NEUTRAL_DESC_OK.search(hit_s)):
        return {"功能类别": "转座子/重复元件", "决定": "drop",
                "决定依据": "注释为转座子/逆转座子元件（TE 结构域或元件命名），直接剔除"}
    if TE_PAT.search(text_own) and NEUTRAL_DESC_OK.search(hit_s) and not TE_PAT.search(hit_s):
        return {"功能类别": "转座子相关(同源提示非TE)", "决定": "neutral",
                "决定依据": "自身注释为 TE 但同源基因有非 TE 明确描述，降级为 neutral 观察"}

    # 2. 未知功能：自身无信息 且 同源无明确描述 → 剔除
    own_unknown = (not desc_s.strip()) or bool(UNKNOWN_PAT.match(desc_s.strip())) or \
                  (desc_s.strip().lower() in {"na", "nan", ".", "-"})
    hit_unknown = (not hit_s.strip()) or bool(UNKNOWN_PAT.match(hit_s.strip())) or \
                  (hit_s.strip().lower() in {"na", "nan", ".", "-"})
    if own_unknown and not (NEUTRAL_DESC_OK.search(pfam_s)):
        if hit_unknown:
            return {"功能类别": "未知功能", "决定": "drop",
                    "决定依据": "自身与同源均无明确功能描述（未知功能蛋白），直接剔除"}
        # 例外：同源有明确描述 → 按同源描述再分类
        for cat, pat in PRIORITY_PATTERNS.items():
            if pat.search(hit_s):
                return {"功能类别": cat, "决定": "priority",
                        "决定依据": f"自身未知但同源基因有明确描述→按同源归入 {cat}"}
        return {"功能类别": "同源提示功能(非优先类)", "决定": "neutral",
                "决定依据": "自身未知，同源基因有明确描述（非优先类别），保留观察"}

    # 3. 优先类别
    for cat, pat in PRIORITY_PATTERNS.items():
        if pat.search(text_all):
            return {"功能类别": cat, "决定": "priority",
                    "决定依据": f"功能注释命中优先类别：{cat}"}

    # 4. 其余
    if own_unknown:
        return {"功能类别": "未知功能", "决定": "drop",
                "决定依据": "无明确功能描述（未知功能蛋白），直接剔除"}
    return {"功能类别": "其他已知功能", "决定": "neutral",
            "决定依据": "有功能注释但不属于优先类别（激素/叶片发育/细胞分裂/转录调控）"}


# ----------------------------------------------------------------------------
# 主流程
# ----------------------------------------------------------------------------
def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="功能注释与优先级排序：同源映射 + 功能过滤 + 证据整合")
    ap.add_argument("--genes-a", required=True,
                    help="方法 A 输出（*.ranked.tsv 或 *.gene_scores.tsv）")
    ap.add_argument("--genes-b", default="",
                    help="方法 B 输出（*.coloc_results.tsv，可选）")
    ap.add_argument("--homologs", required=True,
                    help="同源映射表：gene_id, hit_id, hit_desc, pident, evalue, [source]")
    ap.add_argument("--gene-desc", default="",
                    help="基因功能注释表：gene_id, description, [pfams]（可选）")
    ap.add_argument("--out", required=True, help="输出 TSV")
    ap.add_argument("--config", default="", help="YAML 配置（可选，读取 pp4_threshold）")
    args = ap.parse_args(argv)

    # 1. 方法 A 表
    a = _read(args.genes_a)
    if "gene_id" not in a.columns:
        raise ValueError(f"--genes-a 缺少 gene_id 列（实际: {list(a.columns)}）")

    # 2. 方法 B 证据（可选；区间→基因需列 gene_id 或 interval 映射表合并）
    b = _read(args.genes_b) if args.genes_b else None
    if b is not None and "gene_id" in b.columns:
        keep = [c for c in ("gene_id", "interval", "tool", "colocalized", "evidence",
                            "confidence_label", "data_mode", "degraded_path") if c in b.columns]
        b = b[keep].rename(columns={"evidence": "coloc_evidence",
                                    "confidence_label": "coloc_confidence"})
        a = a.merge(b, on="gene_id", how="left", suffixes=("", "_coloc"))
    elif b is not None:
        print("[WARN] 方法 B 表无 gene_id 列，无法按基因合并；仅按区间字符串匹配 interval 列。",
              file=sys.stderr)
        if "interval" in a.columns:
            a = a.merge(b.rename(columns={"evidence": "coloc_evidence"}),
                        on="interval", how="left", suffixes=("", "_coloc"))

    # 3. 同源映射
    h = _read(args.homologs)
    need_h = {"gene_id", "hit_id", "hit_desc"}
    if not need_h.issubset(h.columns):
        raise ValueError(f"--homologs 缺少列 {need_h - set(h.columns)}")
    h_keep = [c for c in ("gene_id", "hit_id", "hit_desc", "pident", "evalue", "source") if c in h.columns]
    h = h[h_keep].drop_duplicates("gene_id")
    a = a.merge(h, on="gene_id", how="left")

    # 4. 基因自身注释（可选）
    if args.gene_desc:
        d = _read(args.gene_desc)
        d_keep = [c for c in ("gene_id", "description", "pfams") if c in d.columns]
        a = a.merge(d[d_keep].drop_duplicates("gene_id"), on="gene_id", how="left")

    # 5. 分类裁定
    cls_rows = []
    for _, r in a.iterrows():
        cls_rows.append(classify(
            str(r.get("description", "") or ""),
            str(r.get("pfams", "") or ""),
            str(r.get("hit_desc", "") or "")))
    a = pd.concat([a.reset_index(drop=True), pd.DataFrame(cls_rows)], axis=1)

    # 6. 共定位证据列整合（优先 PP.H4/cos，其次联合评分）
    if "coloc_evidence" not in a.columns:
        a["coloc_evidence"] = "NA(未运行方法B)"
    if "joint_z" in a.columns:
        a["joint_score"] = a["joint_z"]
    elif "priority_score" in a.columns:
        a["joint_score"] = a["priority_score"]
    else:
        a["joint_score"] = pd.NA

    # 7. 置信度总标注
    def conf(r) -> str:
        bits = []
        if str(r.get("data_mode", "")) == "truncated" or "截断" in str(r.get("confidence_label", "")):
            bits.append("截断数据模式")
        if str(r.get("coloc_confidence", "")):
            bits.append("coloc:" + str(r.get("coloc_confidence"))[:30])
        if not bits:
            bits.append("完整数据模式")
        return " | ".join(bits)

    a["总置信度标注"] = a.apply(conf, axis=1)

    # 8. 列序 + 输出
    front = ["gene_id", "chr", "start", "end"]
    front = [c for c in front if c in a.columns]
    evid = ["joint_score", "coloc_evidence", "tool", "colocalized", "data_mode"]
    evid = [c for c in evid if c in a.columns]
    annot = ["hit_id", "hit_desc", "pident", "evalue", "source", "description", "pfams"]
    annot = [c for c in annot if c in a.columns]
    verdict = ["功能类别", "决定", "决定依据", "总置信度标注"]
    rest = [c for c in a.columns if c not in front + evid + annot + verdict]
    a = a[front + evid + annot + verdict + rest]

    order = {"priority": 0, "neutral": 1, "drop": 2}
    a = a.sort_values(by="决定", key=lambda s: s.map(order), kind="stable")

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    a.to_csv(args.out, sep="\t", index=False)

    # 汇总
    n = len(a)
    counts = a["决定"].value_counts().to_dict()
    print(f"[OK] 功能优先级排序完成: {n} 基因")
    print(f"     priority={counts.get('priority', 0)}  neutral={counts.get('neutral', 0)}  drop={counts.get('drop', 0)}")
    print(f"     输出: {args.out}")
    if "截断" in " ".join(a["总置信度标注"].astype(str)):
        print("     ⚠ 含截断数据模式结果：优先级排序用途，不可解释为严格统计显著性")
    return 0


if __name__ == "__main__":
    sys.exit(main())
