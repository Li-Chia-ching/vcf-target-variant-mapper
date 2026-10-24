#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_smoke_tests.py — v2.0 回归测试 + 双模式冒烟测试
================================================

A. 回归测试：v1 九脚本 + 模块三脚本 + 测试脚本全部 py_compile 通过；
   R 脚本语法解析通过（plot_effect_distribution.R / run_coloc.R）。
B. 冒烟测试（模拟双性状 summary stats）：
   B1 方法 A 完整模式   → MAGMA(mock 引擎) 正常执行，输出"完整数据模式"标注
   B2 方法 A 截断模式   → 自动切换显著信号富集模式 + 截断标注 + 仅信号基因入排序表
   B3 方法 B 完整区间   → colocboost 正常执行（LD 充分）
   B4 方法 B 截断区间 a → 警告 + 阻断（降级路径 a）
   B5 方法 B 截断区间 b → POEMColoc 降级（降级路径 b）
   B6 方法 B 截断区间 c → 跳过 coloc（降级路径 c）
   B7 功能优先级排序   → priority/drop/neutral 裁定正确

用法（在包根目录或任意位置）：
  python tests/run_smoke_tests.py
"""

import os
import subprocess
import sys

PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TESTS = os.path.join(PKG, "tests")
SIM = os.path.join(TESTS, "simdata")
PY = sys.executable

R_CANDIDATES = [
    os.environ.get("RSCRIPT", ""),
    "Rscript",
    r"C:\Program Files\R\R-4.6.1\bin\Rscript.exe",
]

RESULTS = []


def find_rscript() -> str:
    for c in R_CANDIDATES:
        if not c:
            continue
        try:
            r = subprocess.run([c, "--version"], capture_output=True, text=True, timeout=60)
            if r.returncode == 0:
                return c
        except OSError:
            continue
    return ""


def record(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, ok, detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" — {detail}" if detail and not ok else ""))


def run(cmd, cwd=None, timeout=600):
    return subprocess.run(cmd, cwd=cwd or PKG, capture_output=True, text=True, timeout=timeout)


def check(name: str, cond: bool, detail: str = "") -> None:
    record(name, bool(cond), detail)


# ----------------------------------------------------------------------------
# A. 回归测试
# ----------------------------------------------------------------------------
def test_regression(rscript: str) -> None:
    print("\n=== A. 回归测试（v1 VCF 映射链 + 模块脚本）===")
    py_scripts = [
        "scripts/extract_target_gene_variants.py", "scripts/map_missense_to_protein.py",
        "scripts/align_and_diff.py", "scripts/map_variants_to_domains.py",
        "scripts/structure_conservation.py", "scripts/entropy_bfactor.py",
        "scripts/generate_pymol_scripts.py", "scripts/legend_and_panel.py",
        "modules/multitrait_gene_discovery/run_magma_multitrait.py",
        "modules/multitrait_gene_discovery/functional_prioritization.py",
        "tests/mock_magma.py", "tests/make_simulated_data.py", "tests/run_smoke_tests.py",
    ]
    fails = []
    for s in py_scripts:
        p = os.path.join(PKG, s)
        r = run([PY, "-m", "py_compile", p])
        if r.returncode != 0:
            fails.append(s)
    check("py_compile 全部 Python 脚本", not fails, f"失败: {fails}")

    r = run([PY, os.path.join(TESTS, "make_simulated_data.py")])
    check("生成模拟双性状数据", r.returncode == 0 and os.path.isdir(SIM),
          r.stderr[-300:] if r.returncode else "")

    if rscript:
        r = run([rscript, "-e",
                 "for (f in commandArgs(TRUE)) {tryCatch(parse(file=f), error=function(e) stop(paste(f, conditionMessage(e))))}"] +
                [os.path.join(PKG, "scripts", "plot_effect_distribution.R"),
                 os.path.join(PKG, "modules", "multitrait_gene_discovery", "run_coloc.R")])
        check("R 脚本语法解析（v1 + 模块）", r.returncode == 0, r.stderr[-300:])
    else:
        record("R 脚本语法解析（v1 + 模块）", False, "未找到 Rscript，跳过")


# ----------------------------------------------------------------------------
# B1/B2 方法 A 双模式
# ----------------------------------------------------------------------------
def _write_cfg_a(path: str, mock_magma: bool) -> None:
    mag = (f"  binary: \"{os.path.join(TESTS, 'mock_magma.py').replace(chr(92), '/')}\"\n"
           if mock_magma else "  binary: \"\"\n")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("truncation_detection:\n  p_trunc: 1.0e-5\n  tolerance: 1.5\n"
                 "gene_window:\n  upstream_kb: 2\n  downstream_kb: 1\n"
                 "magma:\n" + mag +
                 "  engine: auto\n  bfile: \"\"\n  model: \"\"\n"
                 "n_per_trait:\n  total: 212\n")


def test_method_a() -> None:
    print("\n=== B1. 方法 A 完整模式（MAGMA mock 引擎 subprocess 管线）===")
    out1 = os.path.join(SIM, "out_a_complete")
    cfg = os.path.join(SIM, "cfg_a_mock.yaml")
    _write_cfg_a(cfg, mock_magma=True)
    r = run([PY, os.path.join(PKG, "modules", "multitrait_gene_discovery", "run_magma_multitrait.py"),
             "--sumstats", os.path.join(SIM, "complete", "ph.complete.tsv"),
             os.path.join(SIM, "complete", "mf.complete.tsv"),
             "--trait-names", "PH", "MF",
             "--genes", os.path.join(SIM, "genes.tsv"),
             "--out-prefix", os.path.join(out1, "gene_scores"),
             "--config", cfg,
             "--mode", "complete"])
    check("B1 方法 A 完整模式运行", r.returncode == 0, r.stderr[-400:])
    gs = os.path.join(out1, "gene_scores.gene_scores.tsv")
    if os.path.exists(gs):
        import pandas as pd
        t = pd.read_csv(gs, sep="\t")
        check("B1 data_mode=complete", (t["data_mode"] == "complete").all())
        check("B1 置信度标注=完整数据模式", t["confidence_label"].str.contains("完整数据模式").all())
        check("B1 产出 magma_z 列（MAGMA 引擎执行）", any(c.startswith("magma_z_") for c in t.columns))
        rep = open(os.path.join(out1, "gene_scores.run_report.tsv"), encoding="utf-8").read()
        check("B1 run_report 记录 magma 引擎", "magma(" in rep or "emulated" in rep,
              "run_report 无引擎记录")
    else:
        check("B1 gene_scores 产物存在", False, "未生成")

    print("\n=== B2. 方法 A 截断模式（p<1e-5 信号表）===")
    out2 = os.path.join(SIM, "out_a_truncated")
    r = run([PY, os.path.join(PKG, "modules", "multitrait_gene_discovery", "run_magma_multitrait.py"),
             "--sumstats", os.path.join(SIM, "truncated", "ph.truncated.tsv"),
             os.path.join(SIM, "truncated", "mf.truncated.tsv"),
             "--trait-names", "PH", "MF",
             "--genes", os.path.join(SIM, "genes.tsv"),
             "--out-prefix", os.path.join(out2, "gene_scores"),
             "--config", cfg,
             "--mode", "truncated"])
    check("B2 方法 A 截断模式运行", r.returncode == 0, r.stderr[-400:])
    check("B2 输出截断警告", "截断数据模式，基因级显著性不可直接解释" in (r.stderr + r.stdout))
    gs = os.path.join(out2, "gene_scores.gene_scores.tsv")
    if os.path.exists(gs):
        import pandas as pd
        t = pd.read_csv(gs, sep="\t")
        check("B2 data_mode=truncated", (t["data_mode"] == "truncated").all())
        check("B2 每行截断标注",
              t["confidence_label"].str.contains("截断数据模式，基因级显著性不可直接解释").all())
        rk = pd.read_csv(os.path.join(out2, "gene_scores.ranked.tsv"), sep="\t")
        genes_with_signal = set(rk["gene_id"])
        check("B2 排序表仅含信号基因", genes_with_signal <= {"G1_DELLA", "G3_TE"},
              f"排序表基因: {sorted(genes_with_signal)}")
    else:
        check("B2 gene_scores 产物存在", False, "未生成")


# ----------------------------------------------------------------------------
# B3–B6 方法 B
# ----------------------------------------------------------------------------
def _write_cfg(path: str, action: str) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"truncation_detection:\n  p_trunc: 1.0e-5\n  tolerance: 1.5\n"
                 f"method_b:\n  primary: colocboost\n  pp4_threshold: 0.8\n"
                 f"  truncation:\n    action: {action}\n"
                 f"n_per_trait:\n  total: 212\n")


def test_method_b(rscript: str) -> None:
    if not rscript:
        for n in ("B3", "B4", "B5", "B6"):
            record(f"{n} 方法 B", False, "未找到 Rscript")
        return
    mod = os.path.join(PKG, "modules", "multitrait_gene_discovery", "run_coloc.R")

    print("\n=== B3/B5. 方法 B：完整区间 colocboost + 截断区间 POEMColoc（auto）===")
    out = os.path.join(SIM, "out_b_auto")
    cfg_auto = os.path.join(SIM, "cfg_auto.yaml")
    _write_cfg(cfg_auto, "auto")
    r = run([rscript, mod,
             "--panels-dir", os.path.join(SIM, "panels"),
             "--ld-dir", os.path.join(SIM, "ld"),
             "--traits", "PH,MF",
             "--out-prefix", os.path.join(out, "coloc"),
             "--config", cfg_auto], timeout=1200)
    check("B3/B5 方法 B 运行（auto）", r.returncode == 0, r.stderr[-500:])
    res = os.path.join(out, "coloc.coloc_results.tsv")
    if os.path.exists(res):
        import pandas as pd
        t = pd.read_csv(res, sep="\t").set_index("interval")
        if "panel_complete" in t.index:
            row = t.loc["panel_complete"]
            check("B3 完整区间走 colocboost", "colocboost" in str(row["tool"]),
                  f"tool={row['tool']}")
            check("B3 完整区间数据模式=complete", row["data_mode"] == "complete")
        else:
            check("B3 完整区间记录存在", False, f"intervals={list(t.index)}")
        if "panel_trunc" in t.index:
            row = t.loc["panel_trunc"]
            check("B5 截断区间降级 POEMColoc", "POEMColoc" in str(row["tool"]),
                  f"tool={row['tool']}")
            check("B5 降级路径标注", "truncated→poemcoloc" in str(row["degraded_path"]),
                  f"degraded_path={row['degraded_path']}")
            check("B5 截断置信度标注", "截断数据模式" in str(row["confidence_label"]))
        else:
            check("B5 截断区间记录存在", False, f"intervals={list(t.index)}")
    else:
        check("B3/B5 coloc_results 产物存在", False, "未生成")

    print("\n=== B4. 方法 B 截断区间：降级路径 (a) 警告 + 阻断 ===")
    out = os.path.join(SIM, "out_b_warn")
    cfg_warn = os.path.join(SIM, "cfg_warn.yaml")
    _write_cfg(cfg_warn, "warn")
    r = run([rscript, mod,
             "--panels-dir", os.path.join(SIM, "panels"),
             "--ld-dir", os.path.join(SIM, "ld"),
             "--traits", "PH,MF",
             "--out-prefix", os.path.join(out, "coloc"),
             "--config", cfg_warn], timeout=1200)
    check("B4 方法 B 运行（warn）", r.returncode == 0, r.stderr[-500:])
    check("B4 输出明确警告", "拒绝静默运行" in (r.stdout + r.stderr) or "降级路径 (a)" in (r.stdout + r.stderr))
    res = os.path.join(out, "coloc.coloc_results.tsv")
    if os.path.exists(res):
        import pandas as pd
        t = pd.read_csv(res, sep="\t").set_index("interval")
        if "panel_trunc" in t.index:
            row = t.loc["panel_trunc"]
            check("B4 截断区间被阻断", "blocked" in str(row["tool"]) or "warn_block" in str(row["degraded_path"]),
                  f"tool={row['tool']}, degraded={row['degraded_path']}")
        else:
            check("B4 截断区间记录存在", False)
    else:
        check("B4 coloc_results 产物存在", False)

    print("\n=== B6. 方法 B 截断区间：降级路径 (c) 跳过 coloc ===")
    out = os.path.join(SIM, "out_b_skip")
    cfg_skip = os.path.join(SIM, "cfg_skip.yaml")
    _write_cfg(cfg_skip, "skip")
    r = run([rscript, mod,
             "--panels-dir", os.path.join(SIM, "panels"),
             "--ld-dir", os.path.join(SIM, "ld"),
             "--traits", "PH,MF",
             "--out-prefix", os.path.join(out, "coloc"),
             "--config", cfg_skip], timeout=1200)
    check("B6 方法 B 运行（skip）", r.returncode == 0, r.stderr[-500:])
    res = os.path.join(out, "coloc.coloc_results.tsv")
    if os.path.exists(res):
        import pandas as pd
        t = pd.read_csv(res, sep="\t").set_index("interval")
        if "panel_trunc" in t.index:
            row = t.loc["panel_trunc"]
            check("B6 截断区间跳过 coloc", "skipped" in str(row["tool"]),
                  f"tool={row['tool']}")
            check("B6 方法 A 为主证据标注", "方法 A" in str(row["evidence"]),
                  f"evidence={row['evidence']}")
        else:
            check("B6 截断区间记录存在", False)
    else:
        check("B6 coloc_results 产物存在", False)


# ----------------------------------------------------------------------------
# B7 功能优先级排序
# ----------------------------------------------------------------------------
def test_functional() -> None:
    print("\n=== B7. 功能优先级排序（同源映射 + 功能过滤）===")
    out = os.path.join(SIM, "out_c")
    gs = os.path.join(SIM, "out_a_truncated", "gene_scores.gene_scores.tsv")
    r = run([PY, os.path.join(PKG, "modules", "multitrait_gene_discovery", "functional_prioritization.py"),
             "--genes-a", gs,
             "--homologs", os.path.join(SIM, "homologs.tsv"),
             "--gene-desc", os.path.join(SIM, "gene_desc.tsv"),
             "--out", os.path.join(out, "prioritized_candidates.tsv")])
    check("B7 功能排序运行", r.returncode == 0, r.stderr[-400:])
    f = os.path.join(out, "prioritized_candidates.tsv")
    if os.path.exists(f):
        import pandas as pd
        t = pd.read_csv(f, sep="\t").set_index("gene_id")
        verdict = t["决定"].to_dict()
        check("B7 G1 DELLA → priority", verdict.get("G1_DELLA") == "priority", f"{verdict.get('G1_DELLA')}")
        check("B7 G3 TE → drop", verdict.get("G3_TE") == "drop", f"{verdict.get('G3_TE')}")
        check("B7 G4 未知功能 → drop", verdict.get("G4_UNKNOWN") == "drop", f"{verdict.get('G4_UNKNOWN')}")
        check("B7 G2 激酶 → neutral", verdict.get("G2_KINASE") == "neutral", f"{verdict.get('G2_KINASE')}")
        check("B7 G5 同源提示 TCP → priority", verdict.get("G5_TCP") == "priority", f"{verdict.get('G5_TCP')}")
        check("B7 含同源 ID 与过滤理由列",
              all(c in t.columns for c in ("hit_id", "hit_desc", "决定依据")))
    else:
        check("B7 输出产物存在", False)


# ----------------------------------------------------------------------------
def main() -> int:
    print("vcf-target-variant-mapper v2.0 — 回归 + 双模式冒烟测试")
    print(f"包目录: {PKG}")
    rscript = find_rscript()
    print(f"Python: {PY}\nRscript: {rscript or '(未找到)'}")

    test_regression(rscript)
    test_method_a()
    test_method_b(rscript)
    test_functional()

    n_pass = sum(1 for _, ok, _ in RESULTS if ok)
    n_fail = len(RESULTS) - n_pass
    print("\n" + "=" * 60)
    print(f"测试汇总: {n_pass} PASS / {n_fail} FAIL / 共 {len(RESULTS)}")
    if n_fail:
        print("失败项:")
        for name, ok, detail in RESULTS:
            if not ok:
                print(f"  ✗ {name} — {detail}")
    print("=" * 60)
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
