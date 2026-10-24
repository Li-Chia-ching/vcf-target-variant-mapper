#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stage 5 — 参考 vs 突变体蛋白全局比对与逐位差异报告。

用法:
    python align_and_diff.py --ref ref.faa --mut mut.faa --out align_diff.txt

使用 Biopython PairwiseAligner（global / BLOSUM62），报告替换、InDel 与移码判定。
"""
import argparse
import sys

try:
    from Bio import Align
    from Bio.Align import substitution_matrices
except ImportError:
    sys.exit("[ERROR] 缺少 biopython，请先: pip install biopython")


def read_fasta(path):
    seq = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    break
            elif line:
                seq.append(line)
    return "".join(seq)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ref", required=True, help="参考蛋白 FASTA")
    ap.add_argument("--mut", required=True, help="突变体蛋白 FASTA")
    ap.add_argument("--out", required=True, help="输出报告路径")
    args = ap.parse_args()

    ref, mut = read_fasta(args.ref), read_fasta(args.mut)
    if not ref or not mut:
        sys.exit("[ERROR] FASTA 读取为空")

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1

    aln = aligner.align(ref, mut)[0]
    a, b = str(aln[0]), str(aln[1])

    diffs, ref_pos, mut_pos = [], 0, 0
    subs = indels = 0
    for x, y in zip(a, b):
        if x != "-":
            ref_pos += 1
        if y != "-":
            mut_pos += 1
        if x != y:
            if x != "-" and y != "-":
                diffs.append(f"SUB   ref{ref_pos}{x} -> mut{mut_pos}{y}")
                subs += 1
            elif x == "-":
                diffs.append(f"INS   mut{mut_pos}{y} (ref 位点 {ref_pos} 前)")
                indels += 1
            else:
                diffs.append(f"DEL   ref{ref_pos}{x} (mut 位点 {mut_pos} 后)")
                indels += 1

    total_gap = sum(1 for x, y in zip(a, b) if x == "-" or y == "-")
    frameshift = (total_gap % 3 != 0) or any(
        1 for x, y in zip(a, b) if x == "-" or y == "-"
    ) and (len(ref) - len(mut)) % 3 != 0

    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write("# 参考 vs 突变体蛋白差异报告\n")
        fh.write(f"# ref_len={len(ref)}  mut_len={len(mut)}\n")
        fh.write(f"# substitutions={subs}  indels={indels}  aligned_len={len(a)}\n")
        fh.write(f"# non-triplet-length-change={'YES' if (len(ref)-len(mut)) % 3 else 'NO'}\n")
        fh.write(f"# frameshift_suspected={'YES' if frameshift else 'NO'}\n")
        fh.write("## 比对\n")
        fh.write(a + "\n")
        fh.write(b + "\n")
        fh.write("## 逐位差异\n")
        if diffs:
            fh.write("\n".join(diffs) + "\n")
        else:
            fh.write("(无差异)\n")

    print(f"[OK] 替换 {subs} / InDel {indels} -> {args.out}")
    if frameshift:
        print("[WARN] 疑似移码（InDel 非 3 倍数），请人工核查翻译框")


if __name__ == "__main__":
    main()
