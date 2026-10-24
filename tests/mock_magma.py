#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mock_magma.py — 冒烟测试用 MAGMA CLI 仿真器
==========================================

模拟真实 MAGMA 的两步调用协议（供 tests/run_smoke_tests.py 验证
run_magma_multitrait.py 的 subprocess 管线完整性）：

  第 1 步  --annotate --snp-loc=... --gene-loc=... --out=...
           产出注释文件（内容不校验，仅需存在）
  第 2 步  --pval=... N=... --gene-annot=... --out=...
           产出 <out>.genes.out，列为 GENE CHR START STOP NSNPS NPARAM N ZSTAT

ZSTAT 由基因名确定性伪随机生成（种子=基因名哈希），保证可复现。
⚠️ 仅用于管线测试；正式分析必须使用真实 MAGMA（ctg.cncr.nl）。
"""

import hashlib
import os
import sys


def parse_args(argv):
    kv = {}
    flags = set()
    for a in argv:
        if "=" in a and a.startswith("--"):
            k, v = a[2:].split("=", 1)
            kv[k] = v
        elif a.startswith("--"):
            flags.add(a[2:])
    return kv, flags


def fake_z(gene: str) -> float:
    h = hashlib.md5(gene.encode()).hexdigest()
    return round((int(h[:8], 16) % 2000) / 100.0 - 5.0, 4)  # [-5, 15)


def main(argv):
    kv, flags = parse_args(argv)
    if "annotate" in flags:
        out = kv.get("out", "genes.annot")
        gene_loc = kv.get("gene-loc", "")
        with open(out, "w", encoding="utf-8") as fh:
            fh.write("// mock magma annotate output\n")
            if gene_loc and os.path.exists(gene_loc):
                with open(gene_loc, encoding="utf-8") as g:
                    fh.write(g.read())
        print(f"[mock-magma] annotate -> {out}")
        return 0

    if "pval" in kv:
        out = kv.get("out", "trait_genes")
        annot = kv.get("gene-annot", "")
        genes = []
        if annot and os.path.exists(annot):
            with open(annot, encoding="utf-8") as fh:
                for line in fh:
                    if line.startswith("//") or not line.strip():
                        continue
                    parts = line.split()
                    if parts:
                        genes.append(parts[0])
        with open(out + ".genes.out", "w", encoding="utf-8") as fh:
            fh.write("GENE CHR START STOP NSNPS NPARAM N ZSTAT\n")
            for g in genes:
                fh.write(f"{g} 1 0 0 10 1 212 {fake_z(g)}\n")
        print(f"[mock-magma] gene analysis -> {out}.genes.out ({len(genes)} genes)")
        return 0

    print(f"[mock-magma] 未识别调用: {argv}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
