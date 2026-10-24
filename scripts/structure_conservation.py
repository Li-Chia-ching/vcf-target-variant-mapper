#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stage 7 — 计算逐残基溶剂可及性 (ASA, Shrake-Rupley) 与多序列比对 Shannon 熵。

用法:
    python structure_conservation.py --pdb model.pdb \
        --msa family.aa.msa --out entropy_asa.tsv

阈值（先验常规）: ASA < 20 Å² 判埋藏；熵 < 0.2 严格保守，< 1.0 中度保守。
PDB 残基编号与 MSA 序列编号偏移用 --offset 校正。
"""
import argparse
import math
import sys
from collections import Counter

try:
    from Bio.PDB import PDBParser
    from Bio.PDB.SASA import ShrakeRupley
    from Bio import AlignIO
except ImportError:
    sys.exit("[ERROR] 缺少 biopython，请先: pip install biopython")

AA = "ACDEFGHIKLMNPQRSTVWY"


def shannon_entropy(col):
    n = len(col)
    if n == 0:
        return 0.0
    counts = Counter(c for c in col.upper() if c in AA)
    total = sum(counts.values())
    if total == 0:
        return 0.0
    ent = 0.0
    for c in counts.values():
        p = c / total
        ent -= p * math.log2(p)
    return ent


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb", required=True, help="结构 PDB（AlphaFold 等）")
    ap.add_argument("--msa", required=True, help="多序列比对 (fasta/clustal 等)")
    ap.add_argument("--out", required=True, help="输出 TSV")
    ap.add_argument("--msa-format", default="fasta", help="AlignIO 格式，默认 fasta")
    ap.add_argument("--offset", type=int, default=0,
                    help="MSA 位置 = PDB 残基号 + offset")
    ap.add_argument("--asa-buried", type=float, default=20.0, help="埋藏 ASA 阈值")
    ap.add_argument("--ent-strict", type=float, default=0.2)
    ap.add_argument("--ent-moderate", type=float, default=1.0)
    args = ap.parse_args()

    # ASA
    structure = PDBParser(QUIET=True).get_structure("model", args.pdb)
    model = next(structure.get_models())
    sr = ShrakeRupley()
    sr.compute(model, level="R")
    asa = {}
    for chain in model:
        for res in chain:
            if res.id[0] != " ":
                continue
            asa[res.id[1]] = res.sasa

    # 熵
    aln = AlignIO.read(args.msa, args.msa_format)
    msa_len = aln.get_alignment_length()
    entropy = {i + 1: shannon_entropy(aln[:, i]) for i in range(msa_len)}
    print(f"[INFO] MSA: {len(aln)} 条序列 x {msa_len} 位; PDB 残基 {len(asa)} 个")

    n_strict = n_mod = n_buried = 0
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write("resi\taa\tASA_A2\tburied\tentropy\tconservation\n")
        for resi in sorted(asa):
            aa = ""
            for chain in model:
                if resi in chain and chain[resi].id[0] == " ":
                    aa = chain[resi].get_resname()
                    break
            ent = entropy.get(resi + args.offset)
            ent_s = f"{ent:.4f}" if ent is not None else "NA"
            if ent is None:
                cons = "NA"
            elif ent < args.ent_strict:
                cons = "strict"; n_strict += 1
            elif ent < args.ent_moderate:
                cons = "moderate"; n_mod += 1
            else:
                cons = "variable"
            buried = "YES" if asa[resi] < args.asa_buried else "NO"
            if buried == "YES":
                n_buried += 1
            fh.write(f"{resi}\t{aa}\t{asa[resi]:.2f}\t{buried}\t{ent_s}\t{cons}\n")

    print(f"[OK] -> {args.out}")
    print(f"     埋藏(<{args.asa_buried}Å²): {n_buried} | "
          f"严格保守(<{args.ent_strict}): {n_strict} | "
          f"中度保守(<{args.ent_moderate}): {n_mod}")


if __name__ == "__main__":
    main()
