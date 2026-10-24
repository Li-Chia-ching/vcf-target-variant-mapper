#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_simulated_data.py — 生成双性状模拟汇总统计（完整 + 截断两种模式）
================================================================

产物（tests/simdata/）：
  complete/ph.complete.tsv, mf.complete.tsv   完整 p 值分布（含 1 个共享强信号区）
  truncated/ph.truncated.tsv, mf.truncated.tsv 仅保留 p<1e-5 行（华智截断风格）
  genes.tsv                                   5 个基因（含信号基因/中性/TE/未知）
  panels/panel_complete.tsv, panel_trunc.tsv  区间面板（snp,pos,af,PH.*,MF.*）
  ld/panel_complete.ld.tsv, panel_trunc.ld.tsv  LD 矩阵（AR(1) r=0.7，正定）
  homologs.tsv, gene_desc.tsv                 同源映射与功能注释（供 M3 测试）

模拟设计：
  chr1:52000-58000 为 PH+MF 共享强信号（p~1e-8），落入 G1（DELLA 同源→priority）
  chr1:250000 附近 仅 PH 信号（非共享），落入 G3（TE→drop）
  其余 p ~ U(0,1)；G4 未知功能（→drop），G2 激酶（→neutral），G5 同源提示 TCP（→priority）
"""

import os

import numpy as np
import pandas as pd

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "simdata")
RNG = np.random.default_rng(20261024)


def make_sumstats(complete: bool) -> None:
    rows = []
    # 背景 SNP：chr1 40-70kb 与 240-270kb，chr2 90-115kb 与 290-320kb，步长 1 kb
    regions = [("chr1", 40000, 70000), ("chr1", 240000, 270000),
               ("chr2", 90000, 115000), ("chr2", 290000, 320000)]
    for chrom, lo, hi in regions:
        for pos in range(lo, hi + 1, 1000):
            rows.append({"chrom": chrom, "pos": pos})
    df = pd.DataFrame(rows)
    df["snp"] = df["chrom"] + "_" + df["pos"].astype(str)
    n = len(df)
    # 基础 p 值
    df["PH_p"] = RNG.uniform(1e-6, 1.0, n)
    df["MF_p"] = RNG.uniform(1e-6, 1.0, n)
    # 植入共享强信号 chr1:55000±2000（G1 区）
    shared = (df["chrom"] == "chr1") & (df["pos"].between(53000, 57000))
    df.loc[shared, "PH_p"] = 10 ** RNG.uniform(-9, -6, shared.sum())
    df.loc[shared, "MF_p"] = 10 ** RNG.uniform(-9, -6, shared.sum())
    # 植入 PH-only 信号 chr1:250000±2000（G3 区，非共享）
    ph_only = (df["chrom"] == "chr1") & (df["pos"].between(248000, 252000))
    df.loc[ph_only, "PH_p"] = 10 ** RNG.uniform(-8, -6, ph_only.sum())
    df.loc[ph_only, "MF_p"] = RNG.uniform(0.1, 1.0, ph_only.sum())

    for t in ("PH", "MF"):
        df[f"{t}_beta"] = RNG.normal(0, 0.05, n)
        df[f"{t}_se"] = np.abs(RNG.normal(0.02, 0.005, n)) + 0.005
        df["af"] = RNG.uniform(0.05, 0.95, n)

    if not complete:  # 华智截断：仅保留任一性状 p<1e-5 的位点
        df = df[(df["PH_p"] < 1e-5) | (df["MF_p"] < 1e-5)]

    sub = "complete" if complete else "truncated"
    d = os.path.join(OUT, sub)
    os.makedirs(d, exist_ok=True)
    out_cols = ["snp", "chrom", "pos", "af"]
    for t in ("PH", "MF"):
        out_cols += [f"{t}_p", f"{t}_beta", f"{t}_se"]
    df[out_cols].to_csv(os.path.join(d, f"ph.{sub}.tsv"), sep="\t", index=False,
                        columns=["snp", "chrom", "pos", "af", "PH_p", "PH_beta", "PH_se"])
    df[out_cols].to_csv(os.path.join(d, f"mf.{sub}.tsv"), sep="\t", index=False,
                        columns=["snp", "chrom", "pos", "af", "MF_p", "MF_beta", "MF_se"])


def make_genes() -> None:
    genes = pd.DataFrame([
        {"gene_id": "G1_DELLA", "chr": "chr1", "start": 50000, "end": 60000, "strand": "+"},
        {"gene_id": "G2_KINASE", "chr": "chr2", "start": 100000, "end": 105000, "strand": "-"},
        {"gene_id": "G3_TE", "chr": "chr1", "start": 245000, "end": 255000, "strand": "+"},
        {"gene_id": "G4_UNKNOWN", "chr": "chr2", "start": 300000, "end": 305000, "strand": "+"},
        {"gene_id": "G5_TCP", "chr": "chr1", "start": 150000, "end": 155000, "strand": "+"},
    ])
    genes.to_csv(os.path.join(OUT, "genes.tsv"), sep="\t", index=False)


def make_panels_and_ld() -> None:
    """区间面板：共享信号区 53000-57000（完整/截断两版）+ AR(1) LD。"""
    for sub, tag in (("complete", "complete"), ("truncated", "trunc")):
        ph = pd.read_csv(os.path.join(OUT, sub, f"ph.{sub}.tsv"), sep="\t")
        mf = pd.read_csv(os.path.join(OUT, sub, f"mf.{sub}.tsv"), sep="\t")
        m = ph.merge(mf, on=["snp", "chrom", "pos", "af"], suffixes=("_PH", "_MF"))
        panel = m[(m["chrom"] == "chr1") & (m["pos"].between(52000, 58000))].copy()
        pdir = os.path.join(OUT, "panels")
        os.makedirs(pdir, exist_ok=True)
        p = pd.DataFrame({
            "snp": panel["snp"], "pos": panel["pos"], "af": panel["af"],
            "PH.beta": panel["PH_beta"], "PH.se": panel["PH_se"], "PH.p": panel["PH_p"],
            "MF.beta": panel["MF_beta"], "MF.se": panel["MF_se"], "MF.p": panel["MF_p"],
        })
        name = "panel_complete" if tag == "complete" else "panel_trunc"
        p.to_csv(os.path.join(pdir, f"{name}.tsv"), sep="\t", index=False)
        # LD：AR(1) r=0.7（to_numpy 防止 insert 按索引对齐写空）
        k = len(p)
        idx = np.arange(k)
        ld = 0.7 ** np.abs(idx[:, None] - idx[None, :])
        ldf = pd.DataFrame(ld, columns=p["snp"].to_numpy())
        ldf.insert(0, "snp", p["snp"].to_numpy())
        ldir = os.path.join(OUT, "ld")
        os.makedirs(ldir, exist_ok=True)
        ldf.to_csv(os.path.join(ldir, f"{name}.ld.tsv"), sep="\t", index=False)


def make_functional_inputs() -> None:
    pd.DataFrame([
        {"gene_id": "G1_DELLA", "hit_id": "LOC_00017224", "hit_desc": "DELLA protein (gibberellin signaling)", "pident": 78.5, "evalue": 1e-120, "source": "M.truncatula"},
        {"gene_id": "G2_KINASE", "hit_id": "LOC_00033110", "hit_desc": "receptor-like protein kinase", "pident": 65.2, "evalue": 1e-80, "source": "M.truncatula"},
        {"gene_id": "G3_TE", "hit_id": "LOC_00090001", "hit_desc": "Ty1-copia retrotransposon gag protein", "pident": 88.0, "evalue": 1e-50, "source": "M.truncatula"},
        {"gene_id": "G4_UNKNOWN", "hit_id": "AT1G00000", "hit_desc": "hypothetical protein", "pident": 30.1, "evalue": 2.0, "source": "A.thaliana"},
        {"gene_id": "G5_TCP", "hit_id": "LOC_00021086", "hit_desc": "TCP transcription factor 4", "pident": 72.3, "evalue": 1e-90, "source": "M.truncatula"},
    ]).to_csv(os.path.join(OUT, "homologs.tsv"), sep="\t", index=False)
    pd.DataFrame([
        {"gene_id": "G1_DELLA", "description": "DELLA protein", "pfams": "GRAS"},
        {"gene_id": "G2_KINASE", "description": "protein kinase", "pfams": "Pkinase"},
        {"gene_id": "G3_TE", "description": "Retrotransposon Ty1-copia", "pfams": "RVT_1 gag_pre-integrs"},
        {"gene_id": "G4_UNKNOWN", "description": "unknown", "pfams": ""},
        {"gene_id": "G5_TCP", "description": "No annotation", "pfams": ""},
    ]).to_csv(os.path.join(OUT, "gene_desc.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    make_sumstats(complete=True)
    make_sumstats(complete=False)
    make_genes()
    make_panels_and_ld()
    make_functional_inputs()
    print(f"[OK] 模拟数据已生成: {OUT}")
