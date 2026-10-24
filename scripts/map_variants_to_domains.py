#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stage 6 — 变异位点映射到线性结构域图（Nature 风格 PDF）。

用法:
    python map_variants_to_domains.py --map protein_map.tsv --out domain_map.pdf

输入为 Stage 4 输出 TSV（需含 aa_pos / aa_ref / aa_alt / motif 列）。
DELLA/GRAS 内置结构域可在 --domains 覆盖。
"""
import argparse
import csv
import sys

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
except ImportError:
    sys.exit("[ERROR] 缺少 matplotlib，请先: pip install matplotlib")

# 内置 DELLA/GRAS 结构域（紫花苜蓿参考编号；按实际蛋白调整）
DEFAULT_DOMAINS = [
    ("N-var", 1, 55, "#999999"),
    ("DELLA", 56, 70, "#D55E00"),
    ("TVHYNP", 81, 95, "#E69F00"),
    ("LHRI", 200, 260, "#0072B2"),
    ("VHIID", 330, 380, "#009E73"),
    ("LHRII", 420, 460, "#56B4E9"),
    ("PFY", 500, 530, "#CC79A7"),
    ("SAW", 540, 570, "#F0E442"),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--map", required=True, help="Stage 4 输出 TSV")
    ap.add_argument("--out", required=True, help="输出 PDF/PNG")
    ap.add_argument("--domains", default="",
                    help='自定义结构域 "名:起-止:#颜色;..." 覆盖内置')
    ap.add_argument("--protein-len", type=int, default=0,
                    help="蛋白长度；0=由结构域/位点推断")
    args = ap.parse_args()

    domains = DEFAULT_DOMAINS
    if args.domains:
        domains = []
        for chunk in args.domains.split(";"):
            nm, rng, color = chunk.split(":")
            s, e = rng.split("-")
            domains.append((nm, int(s), int(e), color))

    variants = []
    with open(args.map, encoding="utf-8") as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            if rec.get("status") != "ok" or not rec.get("aa_pos"):
                continue
            variants.append((int(rec["aa_pos"]),
                             f"{rec['aa_ref']}{rec['aa_pos']}{rec['aa_alt']}"))

    plen = args.protein_len or max(
        [max(e for _, _, e, _ in domains)] + [p for p, _ in variants] + [1]
    )

    plt.rcParams.update({"font.family": "Arial", "font.size": 8})
    fig, ax = plt.subplots(figsize=(7, 1.6))

    # 主链
    ax.add_patch(Rectangle((1, 0.35), plen - 1, 0.3, facecolor="#EEEEEE",
                           edgecolor="black", linewidth=0.6, zorder=1))
    # 结构域
    for nm, s, e, color in domains:
        ax.add_patch(Rectangle((s, 0.35), e - s + 1, 0.3, facecolor=color,
                               edgecolor="black", linewidth=0.5, zorder=2))
        ax.text((s + e) / 2, 0.80, nm, ha="center", va="bottom", fontsize=7)

    # 变异位点
    for i, (pos, label) in enumerate(sorted(variants)):
        ax.plot([pos, pos], [0.20, 0.65], color="#D55E00", linewidth=1.0, zorder=3)
        ax.text(pos, 0.02, label, ha="center", va="top", fontsize=5.5,
                rotation=90, color="#D55E00")

    ax.set_xlim(0, plen * 1.02)
    ax.set_ylim(-0.8, 1.1)
    ax.set_yticks([])
    ax.set_xlabel("Amino acid position")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.8)
    fig.tight_layout()
    fig.savefig(args.out, dpi=600, bbox_inches="tight")
    print(f"[OK] 结构域图 -> {args.out}  (蛋白长 {plen} aa, 变异 {len(variants)} 个)")


if __name__ == "__main__":
    main()
