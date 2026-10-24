#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stage 4 — 将 HGVS.p 错义变异映射到参考蛋白序列并校验。

用法:
    python map_missense_to_protein.py --csv variants.csv \
        --ref-fasta protein.faa --out protein_map.tsv

HGVS.p 三字母格式正则: p\\.([A-Z][a-z]{2})(\\d+)([A-Z][a-z]{2})
位点必须与参考蛋白 FASTA 交叉验证，不一致记录 rejected_reason。
"""
import argparse
import csv
import re
import sys

HGVS_P_RE = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")

AA_3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}

# 可选先验基序（DELLA/GRAS，紫花苜蓿参考编号；按实际蛋白调整）
MOTIFS = [
    ("DELLA", 56, 70),
    ("TVHYNP", 81, 95),
    ("LHRI", 200, 260),
    ("VHIID", 330, 380),
    ("LHRII", 420, 460),
    ("PFY", 500, 530),
    ("SAW", 540, 570),
]


def read_fasta(path):
    seq, name = [], None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    break
                name = line[1:].split()[0]
            elif line:
                seq.append(line)
    return name, "".join(seq)


def motif_of(pos):
    for name, start, end in MOTIFS:
        if start <= pos <= end:
            return name
    return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", required=True, help="Stage 2 输出的变异 CSV")
    ap.add_argument("--ref-fasta", required=True, help="参考蛋白 FASTA")
    ap.add_argument("--out", required=True, help="输出 TSV")
    ap.add_argument("--motifs", default="",
                    help='可选: 自定义基序 "名称:起-止;名称:起-止"，覆盖内置 MOTIFS')
    args = ap.parse_args()

    motifs = MOTIFS
    if args.motifs:
        motifs = []
        for chunk in args.motifs.split(";"):
            nm, rng = chunk.split(":")
            s, e = rng.split("-")
            motifs.append((nm, int(s), int(e)))

    _, ref_seq = read_fasta(args.ref_fasta)
    print(f"[INFO] 参考蛋白长度: {len(ref_seq)} aa")

    rows, n_ok, n_rej = [], 0, 0
    with open(args.csv, encoding="utf-8-sig") as fh:
        for rec in csv.DictReader(fh):
            hgvs = (rec.get("HGVS_p") or "").strip()
            if not hgvs:
                continue
            m = HGVS_P_RE.match(hgvs)
            if not m:
                rows.append({**rec, "aa_pos": "", "aa_ref": "", "aa_alt": "",
                             "motif": "", "status": "rejected",
                             "rejected_reason": "HGVS.p 非三字母错义格式"})
                n_rej += 1
                continue
            ref3, pos, alt3 = m.group(1), int(m.group(2)), m.group(3)
            aa_ref = AA_3_TO_1.get(ref3, "X")
            aa_alt = AA_3_TO_1.get(alt3, "X")
            reason = ""
            if pos > len(ref_seq):
                reason = "位点超出蛋白长度"
            elif ref_seq[pos - 1] != aa_ref:
                reason = f"参考不一致(FASTA={ref_seq[pos-1]}, HGVS={aa_ref})"
            status = "rejected" if reason else "ok"
            if status == "ok":
                n_ok += 1
            else:
                n_rej += 1
            rows.append({**rec, "aa_pos": pos, "aa_ref": aa_ref, "aa_alt": aa_alt,
                         "motif": motif_of(pos) if status == "ok" else "",
                         "status": status, "rejected_reason": reason})

    fields = list(rows[0].keys()) if rows else \
        ["CHROM", "POS", "HGVS_p", "aa_pos", "aa_ref", "aa_alt",
         "motif", "status", "rejected_reason"]
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"[OK] 映射完成: {n_ok} 条接受 / {n_rej} 条拒绝 -> {args.out}")
    if n_rej:
        print("[提示] 被拒条目请检查 HGVS.p 格式与参考蛋白版本是否一致")


if __name__ == "__main__":
    main()
