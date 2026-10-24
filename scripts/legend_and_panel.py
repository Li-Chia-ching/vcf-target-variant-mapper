#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stage 9 — 保守性色条图例 + 拼版 + 图注草稿。

用法:
    python legend_and_panel.py --png render.png --out figure1_panel.png \
        --caption figure1_caption.txt

左侧放 PyMOL 渲染图，右侧/下方放 bwr_r 色条图例（蓝=保守 → 红=可变）。
"""
import argparse
import os

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
except ImportError:
    raise SystemExit("[ERROR] 缺少 matplotlib")

try:
    from PIL import Image
except ImportError:
    raise SystemExit("[ERROR] 缺少 pillow: pip install pillow")

BWR = LinearSegmentedColormap.from_list("bwr_r_custom", ["#2166AC", "#F7F7F7", "#B2182B"])

CAPTION_TEMPLATE = (
    "Figure 1. Conservation mapping of {gene} variants on the predicted 3D structure. "
    "The structure is colored by Shannon entropy derived from a multiple sequence "
    "alignment of the DELLA/GRAS family (blue, strictly conserved, entropy < 0.2; "
    "white, moderately conserved, entropy < 1.0; red, variable). Missense variants "
    "identified from resequencing VCF are shown as sticks. Solvent accessibility was "
    "computed with the Shrake-Rupley algorithm (buried threshold: ASA < 20 Å²). "
    "Rendering: PyMOL, spectrum b blue_white_red, ray 3000×2250."
)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--png", required=True, help="PyMOL 渲染 PNG")
    ap.add_argument("--out", required=True, help="拼版输出 PNG")
    ap.add_argument("--caption", default="", help="图注输出 txt（可选）")
    ap.add_argument("--gene", default="target gene", help="图注中的基因名")
    args = ap.parse_args()

    img = Image.open(args.png).convert("RGB")
    fig = plt.figure(figsize=(8, 4), dpi=600)
    gs = fig.add_gridspec(1, 2, width_ratios=[4, 0.35], wspace=0.08)

    ax0 = fig.add_subplot(gs[0])
    ax0.imshow(img)
    ax0.axis("off")

    ax1 = fig.add_subplot(gs[1])
    grad = [[i / 255.0] for i in range(256)]
    ax1.imshow(grad, aspect="auto", cmap=BWR, origin="lower")
    ax1.set_xticks([])
    ax1.set_yticks([0, 128, 255])
    ax1.set_yticklabels(["Variable", "Moderate", "Conserved"], fontsize=7)
    for side in ("top", "right", "bottom", "left"):
        ax1.spines[side].set_visible(False)

    fig.savefig(args.out, dpi=600, bbox_inches="tight", facecolor="white")
    print(f"[OK] 拼版图 -> {args.out}")

    if args.caption:
        with open(args.caption, "w", encoding="utf-8") as fh:
            fh.write(CAPTION_TEMPLATE.format(gene=args.gene) + "\n")
        print(f"[OK] 图注草稿 -> {args.caption}")
    if os.path.abspath(args.png) == os.path.abspath(args.out):
        print("[WARN] --out 与 --png 相同，输入图被覆盖")


if __name__ == "__main__":
    main()
