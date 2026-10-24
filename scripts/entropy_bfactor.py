#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stage 8a — 将 Shannon 熵写入 PDB 的 B-factor 列（副本操作！）。

用法:
    python entropy_bfactor.py --pdb model.pdb --entropy entropy_asa.tsv \
        --out model_entropy.pdb

⚠️ B-factor 写入会覆盖 AlphaFold 的 pLDDT —— 只对副本操作，原始 PDB 永不改动。
此处按行级文本重写 B-factor（60.46-66.15 列），保留其余原子行格式。
"""
import argparse
import csv
import os
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb", required=True, help="输入 PDB（不会被修改）")
    ap.add_argument("--entropy", required=True, help="Stage 7 输出 TSV（含 resi/entropy）")
    ap.add_argument("--out", required=True, help="输出 PDB（熵值已写入 B-factor）")
    ap.add_argument("--scale", default="auto",
                    help="auto=按熵最大值线性缩放到 0-100；或给定除数，如 2.0")
    args = ap.parse_args()

    if os.path.abspath(args.pdb) == os.path.abspath(args.out):
        sys.exit("[ERROR] --out 不得与 --pdb 相同（会覆盖 pLDDT 原件）")

    ent = {}
    with open(args.entropy, encoding="utf-8") as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            if rec.get("entropy") in ("NA", "", None):
                continue
            ent[int(rec["resi"])] = float(rec["entropy"])

    if not ent:
        sys.exit("[ERROR] 熵表为空")

    if args.scale == "auto":
        max_e = max(ent.values()) or 1.0
        factor = 100.0 / max_e
    else:
        factor = 100.0 / float(args.scale)

    n_wr = 0
    with open(args.pdb, encoding="utf-8") as fin, \
            open(args.out, "w", encoding="utf-8", newline="") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    resi = int(line[22:26])
                except ValueError:
                    fout.write(line)
                    continue
                if resi in ent:
                    b = min(100.0, ent[resi] * factor)
                    line = line[:60] + f"{b:6.2f}" + line[66:]
                    n_wr += 1
            fout.write(line)

    print(f"[OK] {n_wr} 个原子行写入熵 B-factor -> {args.out}")
    print("[提醒] 输出 PDB 的 B-factor 现为熵值（0-100），不再是 pLDDT")


if __name__ == "__main__":
    main()
