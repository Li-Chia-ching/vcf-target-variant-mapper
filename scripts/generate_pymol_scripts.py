#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stage 8b — 生成 PyMOL 渲染脚本 (.pml)。

用法:
    python generate_pymol_scripts.py --pdb model_entropy.pdb --out render.pml

生成的 .pml 用 `pymol -cq render.pml` 执行：
spectrum b blue_white_red、白色背景、ray 3000x2250、输出 PNG。
"""
import argparse
import os

PML_TEMPLATE = """# 自动生成 — 保守性（熵 B-factor）渲染脚本
reinitialize
load {pdb}, model
hide everything
show cartoon, model
color gray80, model
spectrum b, blue_white_red, model
bg_color white
set cartoon_fancy_helices, 1
set ray_trace_mode, 1
set ray_shadows, 0
# 变异位点高亮（可选：预先在 --mutants 中给出 链+残基号，逗号分隔）
{mutant_block}
orient model
ray 3000, 2250
png {png}, dpi=600
quit
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb", required=True, help="熵 B-factor PDB（Stage 8a 输出）")
    ap.add_argument("--out", required=True, help="输出 .pml 路径")
    ap.add_argument("--png", default="", help="渲染 PNG 路径（默认与 .pml 同名）")
    ap.add_argument("--mutants", default="",
                    help='变异位点高亮，如 "A/56/,A/81/"（PyMOL selection 语法）')
    args = ap.parse_args()

    png = args.png or os.path.splitext(args.out)[0] + ".png"
    if args.mutants:
        sel = args.mutants.rstrip(",")
        mutant_block = (
            f"select mutants, ({sel})\n"
            "show sticks, mutants\n"
            "util.cbag mutants\n"
            "set stick_radius, 0.12, mutants"
        )
    else:
        mutant_block = "# 未指定变异位点高亮"

    pml = PML_TEMPLATE.format(pdb=os.path.basename(args.pdb),
                              png=os.path.basename(png),
                              mutant_block=mutant_block)
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(pml)
    print(f"[OK] PyMOL 脚本 -> {args.out}")
    print(f"     执行: pymol -cq {args.out}   (渲染输出 {os.path.basename(png)})")
    print("[提示] .pml 里的 load/png 用文件名，运行前请将 PDB 与 .pml 放同一目录")


if __name__ == "__main__":
    main()
