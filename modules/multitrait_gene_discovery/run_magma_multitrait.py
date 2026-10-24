#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_magma_multitrait.py — 方法 A（第一代）：基因水平多性状联合评分
==============================================================

【方法定位】
geneGPS（GPS-N / GPS-M）式"SNP 信号 → 基因窗口聚合"多性状框架 + MAGMA 基因级检验。
⚠️ geneGPS 在公开检索中未找到对应仓库，本模块按其**概念**实现基因级多性状联合评分
（GPS-N：窗口内显著 SNP 计数；GPS-M：窗口内信号幅度），不绑定任何第三方 geneGPS 包；
严格的基因级 p 值检验由成熟引擎 MAGMA（可执行文件路径可在配置中指定）承担。

【截断数据模式切换（核心逻辑）】
  完整汇总统计（p 值覆盖全区间）→ 标准 MAGMA 基因级分析（SNP-wise mean 模型）；
  截断信号表（仅 p < 1e-5，如华智生物交付）→ 自动切换"显著信号基因级富集模式"：
    - 仅对携带截断信号的基因打分排序（GPS-N/GPS-M），不对全基因组做基因级检验；
    - 结果定位为**优先级排序**而非严格基因级显著性；
    - 每行输出强制携带标注："截断数据模式，基因级显著性不可直接解释"。
  检测逻辑：任一性状 max(p) <= p_trunc * 1.5 判为截断；--mode 可强制指定。

【输入】
  - 每个性状一份汇总统计 TSV（自动识别列名：SNP/CHR/BP/P[/BETA/SE]）
  - 基因坐标表（gene_id chr start end [strand]）
  - 配置 YAML（可选，见 example_config.yaml）

【输出】
  - <out-prefix>.gene_scores.tsv   基因级评分总表（含模式与置信度标注）
  - <out-prefix>.ranked.tsv        按联合评分排序的优先级列表
  - <out-prefix>.run_report.tsv    运行模式、截断检测证据、引擎选择记录

【用法】
  python run_magma_multitrait.py \
      --sumstats data/ph.tsv data/mf.tsv data/sh.tsv \
      --trait-names PH MF SH \
      --genes data/genes.tsv \
      --out-prefix output/multitrait/gene_scores \
      [--config example_config.yaml] [--mode auto|complete|truncated]

作者: vcf-target-variant-mapper v2.0 多性状共有基因发现模块
依赖: pandas, numpy, scipy, pyyaml（MAGMA 二进制可选，见 magma.engine 配置）
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

# ----------------------------------------------------------------------------
# 常量
# ----------------------------------------------------------------------------
P_TRUNC_DEFAULT = 1e-5          # 华智生物截断阈值
TRUNC_TOL = 1.5                  # max(p) <= P_TRUNC*TRUNC_TOL 判为截断
WINDOW_UPSTREAM_KB = 2.0         # 基因窗口：上游 2 kb
WINDOW_DOWNSTREAM_KB = 1.0       # 基因窗口：下游 1 kb
TRUNCATED_LABEL = "截断数据模式，基因级显著性不可直接解释"
COMPLETE_LABEL = "完整数据模式，基因级统计可直接解释"

MODE_AUTO = "auto"
MODE_COMPLETE = "complete"
MODE_TRUNCATED = "truncated"


# ----------------------------------------------------------------------------
# 配置
# ----------------------------------------------------------------------------
@dataclass
class Config:
    """运行配置（YAML 可覆盖，命令行优先级最高）。"""
    mode: str = MODE_AUTO
    p_trunc: float = P_TRUNC_DEFAULT
    window_upstream_kb: float = WINDOW_UPSTREAM_KB
    window_downstream_kb: float = WINDOW_DOWNSTREAM_KB
    magma_binary: str = ""            # MAGMA 可执行文件路径（可配置）
    magma_engine: str = "auto"        # auto | magma | emulated
    magma_bfile: str = ""             # 可选参考基因型 --bfile（LD 校正用）
    magma_model: str = ""             # 空 = MAGMA 默认 snp-wise=mean
    min_snps_per_gene: int = 1        # GPS 评分最小 SNP 数
    n_per_trait: Dict[str, int] = field(default_factory=dict)  # 每性状样本量
    warn_truncated: bool = True

    @classmethod
    def from_yaml(cls, path: str) -> "Config":
        if yaml is None:
            raise RuntimeError("pyyaml 未安装，无法读取配置文件")
        with open(path, "r", encoding="utf-8") as fh:
            raw = yaml.safe_load(fh) or {}
        mag = raw.get("magma", {}) or {}
        win = raw.get("gene_window", {}) or {}
        det = raw.get("truncation_detection", {}) or {}
        cfg = cls(
            mode=raw.get("data_mode", MODE_AUTO),
            p_trunc=float(det.get("p_trunc", P_TRUNC_DEFAULT)),
            window_upstream_kb=float(win.get("upstream_kb", WINDOW_UPSTREAM_KB)),
            window_downstream_kb=float(win.get("downstream_kb", WINDOW_DOWNSTREAM_KB)),
            magma_binary=str(mag.get("binary", "")),
            magma_engine=str(mag.get("engine", "auto")),
            magma_bfile=str(mag.get("bfile", "")),
            magma_model=str(mag.get("model", "")),
            min_snps_per_gene=int(raw.get("min_snps_per_gene", 1)),
            n_per_trait={str(k): int(v) for k, v in (raw.get("n_per_trait", {}) or {}).items()},
            warn_truncated=bool(det.get("warn", True)),
        )
        return cfg


# ----------------------------------------------------------------------------
# 数据加载
# ----------------------------------------------------------------------------
SNP_ALIASES = ("SNP", "snp", "RSID", "rsid", "ID", "id", "variant", "VARIANT")
CHR_ALIASES = ("CHR", "chr", "CHROM", "chrom", "#CHROM")
BP_ALIASES = ("BP", "bp", "POS", "pos", "POSITION", "position", "BP_hg38")
P_ALIASES = ("P", "p", "PVAL", "pval", "p_value", "P_VALUE", "PVALUE")
BETA_ALIASES = ("BETA", "beta", "EFFECT", "effect", "b")
SE_ALIASES = ("SE", "se", "STDERR", "stderr", "sebeta", "SEBETA")


def _pick(cols: Sequence[str], aliases: Sequence[str]) -> Optional[str]:
    upper = {c.upper(): c for c in cols}
    for a in aliases:
        if a.upper() in upper:
            return upper[a.upper()]
    return None


def _pick_suffix(cols: Sequence[str], suffix: str) -> Optional[str]:
    """性状前缀列兜底：匹配 *_p / *.p / *-p（如 PH_p、MF.p）。"""
    pat = re.compile(rf"(^|[._-]){re.escape(suffix)}$", re.I)
    for c in cols:
        if pat.search(c):
            return c
    return None


def load_sumstats(path: str) -> pd.DataFrame:
    """读取性状汇总统计，标准化为 snp/chr/bp/p[/beta/se]。"""
    if not os.path.exists(path):
        raise FileNotFoundError(f"汇总统计不存在: {path}")
    sep = "\t" if path.endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep, comment="#")
    if df.empty:
        raise ValueError(f"汇总统计为空: {path}")
    col_snp = _pick(df.columns, SNP_ALIASES) or _pick_suffix(df.columns, "snp")
    col_chr = _pick(df.columns, CHR_ALIASES) or _pick_suffix(df.columns, "chr") or _pick_suffix(df.columns, "chrom")
    col_bp = _pick(df.columns, BP_ALIASES) or _pick_suffix(df.columns, "bp") or _pick_suffix(df.columns, "pos")
    col_p = _pick(df.columns, P_ALIASES) or _pick_suffix(df.columns, "p") or _pick_suffix(df.columns, "pval")
    missing = [n for n, c in [("SNP", col_snp), ("CHR", col_chr), ("BP", col_bp), ("P", col_p)] if c is None]
    if missing:
        raise ValueError(f"{path} 缺少必要列: {missing}（实际列: {list(df.columns)}）")
    out = pd.DataFrame({
        "snp": df[col_snp].astype(str),
        "chr": df[col_chr].astype(str).str.lower(),
        "bp": pd.to_numeric(df[col_bp], errors="coerce"),
        "p": pd.to_numeric(df[col_p], errors="coerce"),
    })
    col_beta = _pick(df.columns, BETA_ALIASES) or _pick_suffix(df.columns, "beta")
    col_se = _pick(df.columns, SE_ALIASES) or _pick_suffix(df.columns, "se") or _pick_suffix(df.columns, "sebeta")
    if col_beta is not None:
        out["beta"] = pd.to_numeric(df[col_beta], errors="coerce")
    if col_se is not None:
        out["se"] = pd.to_numeric(df[col_se], errors="coerce")
    out = out.dropna(subset=["bp", "p"])
    out = out[(out["p"] > 0) & (out["p"] <= 1)]
    if out.empty:
        raise ValueError(f"{path} 清洗后无有效行（检查 P/BP 列）")
    return out.reset_index(drop=True)


def load_genes(path: str) -> pd.DataFrame:
    """读取基因坐标表：gene_id chr start end [strand]。"""
    if not os.path.exists(path):
        raise FileNotFoundError(f"基因坐标表不存在: {path}")
    df = pd.read_csv(path, sep=None, engine="python", comment="#")
    cols = {c.lower(): c for c in df.columns}
    need = ["gene_id", "chr", "start", "end"]
    for n in need:
        if n not in cols:
            raise ValueError(f"基因表缺少列 {n}（实际: {list(df.columns)}）")
    out = pd.DataFrame({
        "gene_id": df[cols["gene_id"]].astype(str),
        "chr": df[cols["chr"]].astype(str).str.lower(),
        "start": pd.to_numeric(df[cols["start"]], errors="coerce"),
        "end": pd.to_numeric(df[cols["end"]], errors="coerce"),
    })
    strand_col = cols.get("strand")
    out["strand"] = df[strand_col].astype(str) if strand_col else "+"
    out = out.dropna(subset=["start", "end"])
    if out.empty:
        raise ValueError(f"基因表清洗后为空: {path}")
    return out.reset_index(drop=True)


# ----------------------------------------------------------------------------
# 截断检测
# ----------------------------------------------------------------------------
def detect_truncation(per_trait: Dict[str, pd.DataFrame], p_trunc: float) -> Tuple[bool, pd.DataFrame]:
    """判断汇总统计是否为截断信号表。

    判据：某性状全部 p <= p_trunc * TRUNC_TOL（即不含非显著位点）→ 该性状截断；
    任一性状截断 → 整体按截断处理（多性状联合检验要求全部性状完整）。
    """
    rows = []
    any_trunc = False
    for name, df in per_trait.items():
        pmax = float(df["p"].max())
        n_sig = int((df["p"] <= p_trunc).sum())
        frac_sig = n_sig / len(df)
        is_trunc = pmax <= p_trunc * TRUNC_TOL
        any_trunc = any_trunc or is_trunc
        rows.append({
            "trait": name, "n_snps": len(df), "p_min": float(df["p"].min()),
            "p_max": pmax, "n_below_p_trunc": n_sig, "frac_below_p_trunc": round(frac_sig, 4),
            "verdict": "TRUNCATED" if is_trunc else "COMPLETE",
        })
    return any_trunc, pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# 基因窗口与 GPS 评分
# ----------------------------------------------------------------------------
def expand_windows(genes: pd.DataFrame, up_kb: float, down_kb: float) -> pd.DataFrame:
    """按链向扩展基因窗口：+ 链 start-up/end+down；- 链 start-down/end+up。"""
    up = int(up_kb * 1000)
    down = int(down_kb * 1000)
    g = genes.copy()
    minus = g["strand"].astype(str).str.strip() == "-"
    g["win_start"] = np.where(minus, g["start"] - down, g["start"] - up)
    g["win_end"] = np.where(minus, g["end"] + up, g["end"] + down)
    g["win_start"] = g["win_start"].clip(lower=1).astype(np.int64)
    g["win_end"] = g["win_end"].astype(np.int64)
    return g


def assign_snps_to_genes(snps: pd.DataFrame, genes_win: pd.DataFrame) -> pd.DataFrame:
    """把 SNP 按坐标落入基因窗口（按染色体分组区间匹配）。"""
    frames = []
    g_by_chr = dict(tuple(genes_win.groupby("chr")))
    for chrom, s in snps.groupby("chr"):
        g = g_by_chr.get(chrom)
        if g is None or g.empty:
            continue
        s_sorted = s.sort_values("bp")
        bp = s_sorted["bp"].to_numpy()
        for _, row in g.iterrows():
            lo = np.searchsorted(bp, row["win_start"], side="left")
            hi = np.searchsorted(bp, row["win_end"], side="right")
            if hi > lo:
                sub = s_sorted.iloc[lo:hi].copy()
                sub["gene_id"] = row["gene_id"]
                frames.append(sub)
    if not frames:
        return pd.DataFrame(columns=["snp", "chr", "bp", "p", "gene_id"])
    return pd.concat(frames, ignore_index=True)


def gps_scores(assigned: pd.DataFrame, p_trunc: float, min_snps: int) -> pd.DataFrame:
    """GPS-N（窗口内显著 SNP 计数）与 GPS-M（窗口内信号幅度）。

    GPS-M = sum(-log10(p)) of SNPs with p <= p_trunc（幅度聚合）；
    截断模式下两值即为优先级排序依据（非统计检验量）。
    """
    if assigned.empty:
        return pd.DataFrame(columns=["gene_id", "gps_n", "gps_m", "n_snps_window", "min_p"])
    rows = []
    for gid, sub in assigned.groupby("gene_id"):
        sig = sub[sub["p"] <= p_trunc]
        rows.append({
            "gene_id": gid,
            "gps_n": int(len(sig)),
            "gps_m": float(np.sum(-np.log10(sig["p"]))) if len(sig) else 0.0,
            "n_snps_window": int(len(sub)),
            "min_p": float(sub["p"].min()),
        })
    out = pd.DataFrame(rows)
    out = out[out["n_snps_window"] >= min_snps]
    return out


# ----------------------------------------------------------------------------
# MAGMA（真实二进制或内置仿真引擎）
# ----------------------------------------------------------------------------
def find_magma(cfg: Config) -> Tuple[Optional[str], str]:
    """返回 (可执行路径或 None, 引擎描述)。engine=emulated 强制内置引擎。"""
    if cfg.magma_engine == "emulated":
        return None, "emulated(强制)"
    cand = cfg.magma_binary or os.environ.get("MAGMA_BIN", "")
    if cand and os.path.exists(cand):
        return cand, f"magma({cand})"
    for name in ("magma", "magma.exe"):
        hit = shutil.which(name)
        if hit:
            return hit, f"magma(PATH:{hit})"
    if cfg.magma_engine == "magma":
        raise FileNotFoundError(
            "配置要求真实 MAGMA（magma.engine=magma），但未找到可执行文件。"
            "请在配置 magma.binary 或环境变量 MAGMA_BIN 中指定路径。")
    return None, "emulated(magma 二进制缺失，内置 SNP-wise mean 近似)"


def run_magma_binary(magma: str, cfg: Config, sumstats: pd.DataFrame,
                     genes_win: pd.DataFrame, n: int, out_prefix: str,
                     trait: str) -> Optional[pd.DataFrame]:
    """调用真实 MAGMA：--annotate + --pval 基因分析，解析 .genes.out 的 ZSTAT。

    可测性：magma 路径以 .py 结尾时用当前 Python 解释器调用（tests/mock_magma.py
    用于冒烟测试 subprocess 管线；生产环境请指向真实 MAGMA 二进制）。
    """
    tmp = tempfile.mkdtemp(prefix=f"magma_{trait}_")
    exe = [sys.executable, magma] if magma.endswith(".py") else [magma]
    try:
        snp_loc = os.path.join(tmp, "snps.loc")
        gene_loc = os.path.join(tmp, "genes.loc")
        sumstats[["snp", "chr", "bp"]].drop_duplicates("snp").to_csv(snp_loc, sep=" ", header=False, index=False)
        genes_win[["gene_id", "chr", "win_start", "win_end"]].to_csv(gene_loc, sep=" ", header=False, index=False)
        annot = os.path.join(tmp, "genes.annot")
        cmd_a = exe + ["--annotate",
                       f"--snp-loc={snp_loc}", f"--gene-loc={gene_loc}",
                       f"--out={annot}"]
        subprocess.run(cmd_a, check=True, capture_output=True, text=True, timeout=1800)
        pval_file = os.path.join(tmp, "trait.pval")
        sumstats[["snp", "p"]].to_csv(pval_file, sep="\t", index=False)
        out = os.path.join(tmp, "trait_genes")
        cmd_g = exe + [f"--pval={pval_file}", f"N={n}",
                       f"--gene-annot={annot}", f"--out={out}"]
        if cfg.magma_bfile:
            cmd_g.append(f"--bfile={cfg.magma_bfile}")
        if cfg.magma_model:
            cmd_g.append(f"--model={cfg.magma_model}")
        subprocess.run(cmd_g, check=True, capture_output=True, text=True, timeout=3600)
        genes_out = out + ".genes.out"
        if not os.path.exists(genes_out):
            raise RuntimeError(f"MAGMA 未产出 {genes_out}")
        res = pd.read_csv(genes_out, sep=r"\s+")
        col_z = _pick(res.columns, ("ZSTAT", "zstat", "Z"))
        col_g = _pick(res.columns, ("GENE", "gene"))
        if col_z is None or col_g is None:
            raise RuntimeError(f"MAGMA .genes.out 列异常: {list(res.columns)}")
        return pd.DataFrame({
            "gene_id": res[col_g].astype(str),
            f"magma_z_{trait}": pd.to_numeric(res[col_z], errors="coerce"),
        })
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, RuntimeError) as exc:
        sys.stderr.write(f"[WARN] MAGMA 真实调用失败（{trait}）: {exc}\n"
                         f"[WARN] 回退内置 emulated 引擎，Z 值为近似（LD 未校正）。\n")
        return None
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def emulated_gene_z(assigned: pd.DataFrame, genes_win: pd.DataFrame, trait: str) -> pd.DataFrame:
    """内置 SNP-wise mean 近似：gene Z = mean(z_snp) * sqrt(k)（LD 未校正，标注 emulated）。

    仅在 MAGMA 二进制不可用时使用；输出列名带 magma_z_ 前缀但 run_report 中
    标注引擎为 emulated，避免与真实 MAGMA 混淆。
    """
    if assigned.empty:
        return pd.DataFrame(columns=["gene_id", f"magma_z_{trait}"])
    a = assigned.copy()
    # 由 p 值换算单侧 z（效应方向未知时取正，仅用于排序）
    from scipy.stats import norm
    a["z"] = norm.isf(a["p"] / 2.0).clip(0, 10)
    rows = []
    for gid, sub in a.groupby("gene_id"):
        k = len(sub)
        rows.append({"gene_id": gid, f"magma_z_{trait}": float(sub["z"].mean() * np.sqrt(k))})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# 多性状联合
# ----------------------------------------------------------------------------
def joint_scores(score_tab: pd.DataFrame, traits: Sequence[str]) -> pd.DataFrame:
    """跨性状联合：z 列 Stouffer 合并；GPS-N/M 取跨性状求和与携带信号性状数。"""
    out = score_tab.copy()
    z_cols = [f"magma_z_{t}" for t in traits if f"magma_z_{t}" in out.columns]
    gpsn_cols = [f"gps_n_{t}" for t in traits if f"gps_n_{t}" in out.columns]
    gpsm_cols = [f"gps_m_{t}" for t in traits if f"gps_m_{t}" in out.columns]
    if z_cols:
        zmat = out[z_cols].fillna(0.0).to_numpy()
        out["joint_z"] = zmat.sum(axis=1) / np.sqrt(len(z_cols))
    else:
        out["joint_z"] = np.nan
    if gpsn_cols:
        nmat = out[gpsn_cols].fillna(0).to_numpy()
        out["gps_n_sum"] = nmat.sum(axis=1)
        out["n_traits_with_signal"] = (nmat > 0).sum(axis=1)
    if gpsm_cols:
        out["gps_m_sum"] = out[gpsm_cols].fillna(0.0).to_numpy().sum(axis=1)
    # 优先级排序分：截断模式用 GPS 复合；完整模式用 joint_z
    if out["joint_z"].notna().any():
        out["priority_score"] = out["joint_z"].fillna(0)
    else:
        out["priority_score"] = out.get("gps_m_sum", pd.Series(0, index=out.index)).fillna(0)
    return out


# ----------------------------------------------------------------------------
# 主流程
# ----------------------------------------------------------------------------
def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="方法 A：geneGPS 概念 + MAGMA 基因级多性状联合评分（含截断模式切换）")
    ap.add_argument("--sumstats", nargs="+", required=True, help="每性状汇总统计 TSV/CSV")
    ap.add_argument("--trait-names", nargs="+", required=True, help="与 --sumstats 一一对应的性状名")
    ap.add_argument("--genes", required=True, help="基因坐标表 gene_id chr start end [strand]")
    ap.add_argument("--out-prefix", required=True, help="输出前缀")
    ap.add_argument("--config", default="", help="YAML 配置（可选）")
    ap.add_argument("--mode", default="", choices=["", MODE_AUTO, MODE_COMPLETE, MODE_TRUNCATED],
                    help="强制数据模式（默认 auto=自动检测）")
    args = ap.parse_args(argv)

    if len(args.sumstats) != len(args.trait_names):
        print("错误: --sumstats 与 --trait-names 数量不一致", file=sys.stderr)
        return 2

    cfg = Config.from_yaml(args.config) if args.config else Config()
    if args.mode:
        cfg.mode = args.mode

    # 1. 载入
    per_trait = {t: load_sumstats(p) for t, p in zip(args.trait_names, args.sumstats)}
    genes = load_genes(args.genes)
    genes_win = expand_windows(genes, cfg.window_upstream_kb, cfg.window_downstream_kb)

    # 2. 截断检测
    any_trunc, trunc_tab = detect_truncation(per_trait, cfg.p_trunc)
    if cfg.mode == MODE_AUTO:
        mode = MODE_TRUNCATED if any_trunc else MODE_COMPLETE
    else:
        mode = cfg.mode
    label = TRUNCATED_LABEL if mode == MODE_TRUNCATED else COMPLETE_LABEL
    if mode == MODE_TRUNCATED and cfg.warn_truncated:
        sys.stderr.write(
            "\n" + "=" * 78 + "\n"
            "⚠️  检测到截断汇总统计（仅含显著信号位点，如 p<1e-5）。\n"
            f"    本模块自动切换为『显著信号基因级富集模式』。\n"
            f"    {TRUNCATED_LABEL}\n"
            "    结果仅用于候选基因优先级排序，不可当作基因级全基因组检验结果引用。\n"
            "    如需严格基因级显著性，请向数据方索取完整（未截断）汇总统计。\n"
            + "=" * 78 + "\n\n")

    # 3. 引擎选择
    magma, engine_desc = find_magma(cfg)

    # 4. 逐性状：SNP→基因、GPS 评分、基因级 Z
    score_frames: List[pd.DataFrame] = []
    engine_per_trait: Dict[str, str] = {}
    base = genes_win[["gene_id", "chr", "win_start", "win_end", "strand"]].copy()
    for trait, ss in per_trait.items():
        assigned = assign_snps_to_genes(ss, genes_win)
        gps = gps_scores(assigned, cfg.p_trunc, cfg.min_snps_per_gene)
        gps = gps.rename(columns={"gps_n": f"gps_n_{trait}", "gps_m": f"gps_m_{trait}",
                                  "n_snps_window": f"n_snps_window_{trait}",
                                  "min_p": f"min_p_{trait}"})
        score_frames.append(gps.set_index("gene_id"))

        if mode == MODE_TRUNCATED:
            # 截断模式：不做全基因组基因级检验（无非显著位点，检验不可解释）
            engine_per_trait[trait] = "skipped(截断模式→显著信号富集)"
        elif magma is not None:
            res = run_magma_binary(magma, cfg, ss, genes_win,
                                   cfg.n_per_trait.get(trait, 212), args.out_prefix, trait)
            if res is not None:
                score_frames.append(res.set_index("gene_id"))
                engine_per_trait[trait] = engine_desc
            else:
                score_frames.append(emulated_gene_z(assigned, genes_win, trait).set_index("gene_id"))
                engine_per_trait[trait] = "emulated(magma 调用失败回退)"
        else:
            score_frames.append(emulated_gene_z(assigned, genes_win, trait).set_index("gene_id"))
            engine_per_trait[trait] = engine_desc

    # 5. 汇总 + 联合
    tab = base.set_index("gene_id")
    for f in score_frames:
        tab = tab.join(f, how="left")
    tab = tab.reset_index()
    tab = joint_scores(tab, args.trait_names)
    tab["data_mode"] = mode
    tab["confidence_label"] = label
    tab["magma_engine"] = "; ".join(f"{t}={e}" for t, e in engine_per_trait.items())

    # 截断模式：仅保留携带信号的基因进入排序表
    ranked = tab.copy()
    if mode == MODE_TRUNCATED:
        sig_cols = [f"gps_n_{t}" for t in args.trait_names if f"gps_n_{t}" in ranked.columns]
        if sig_cols:
            mask = ranked[sig_cols].fillna(0).sum(axis=1) > 0
            ranked = ranked[mask]
    ranked = ranked.sort_values("priority_score", ascending=False)

    # 6. 输出
    out_dir = os.path.dirname(args.out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    tab.to_csv(f"{args.out_prefix}.gene_scores.tsv", sep="\t", index=False)
    ranked.to_csv(f"{args.out_prefix}.ranked.tsv", sep="\t", index=False)

    report = pd.DataFrame({
        "item": ["data_mode", "confidence_label", "magma_engine", "p_trunc",
                 "window_upstream_kb", "window_downstream_kb", "n_genes_input",
                 "n_genes_scored", "n_genes_ranked"] +
                [f"trunc_check.{r.trait}" for r in trunc_tab.itertuples()],
        "value": [mode, label, engine_desc, cfg.p_trunc,
                  cfg.window_upstream_kb, cfg.window_downstream_kb, len(genes),
                  len(tab), len(ranked)] +
                 [f"{r.verdict}: p_max={r.p_max:.3g}, n={r.n_snps}" for r in trunc_tab.itertuples()],
    })
    report.to_csv(f"{args.out_prefix}.run_report.tsv", sep="\t", index=False)
    trunc_tab.to_csv(f"{args.out_prefix}.truncation_check.tsv", sep="\t", index=False)

    print(f"[OK] 方法 A 完成: mode={mode}, engine={engine_desc}")
    print(f"     gene_scores: {args.out_prefix}.gene_scores.tsv ({len(tab)} 基因)")
    print(f"     ranked:      {args.out_prefix}.ranked.tsv ({len(ranked)} 基因)")
    print(f"     ⚠ {label}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
