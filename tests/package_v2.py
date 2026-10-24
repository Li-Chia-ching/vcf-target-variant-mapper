#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""package_v2.py — 打包 vcf-target-variant-mapper v2.0 zip 并校验关键文件齐备。"""
import os
import zipfile

PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(os.path.dirname(PKG), "vcf-target-variant-mapper.zip")
EXCLUDE_DIR = {"__pycache__", ".git"}

n = 0
with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as z:
    for root, dirs, files in os.walk(PKG):
        dirs[:] = [d for d in dirs if d not in EXCLUDE_DIR]
        for f in sorted(files):
            p = os.path.join(root, f)
            arc = os.path.join("vcf-target-variant-mapper",
                               os.path.relpath(p, PKG).replace("\\", "/"))
            z.write(p, arc)
            n += 1
print(f"zip written: {OUT} ({n} files, {os.path.getsize(OUT)} bytes)")

with zipfile.ZipFile(OUT) as z:
    names = z.namelist()
required = ["SKILL.md", "README.md", "CHANGELOG_v2.0.md", "example_config.yaml",
            "requirements.txt", "DESCRIPTION",
            "run_magma_multitrait.py", "run_coloc.R", "functional_prioritization.py",
            "run_smoke_tests.py", "mock_magma.py", "make_simulated_data.py"]
missing = [r for r in required if not any(x.endswith(r) for x in names)]
if missing:
    raise SystemExit(f"zip 缺少关键文件: {missing}")
print(f"zip validation OK — {len(names)} entries")
