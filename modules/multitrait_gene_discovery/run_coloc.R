#!/usr/bin/env Rscript
# =============================================================================
# run_coloc.R — 方法 B（第二代）：贝叶斯共定位确认 + 截断数据降级策略
# =============================================================================
#
# 【方法定位】
#   优先 colocboost（GitHub: StatFunGen/colocboost）多性状共定位；
#   备选 HyPrColoc（GitHub: jrs95/hyprcoloc，RcppEigen 锁定 0.3.3.9.3）；
#   coloc v5/6 SuSiE（需 LD）仅在用户显式要求时运行。
#
# 【截断数据铁律（不得静默运行）】
#   若仅提供截断信号（p < 1e-5 保留，其余缺失）：
#     先检查输入是否满足 colocboost/HyPrColoc 统计前提（需要完整 p 值分布做
#     先验校准），不满足时按以下降级路径之一处理（truncation.action 配置）：
#       (a) warn    — 明确警告 + 建议获取完整区域统计数据（默认，分析中止）
#       (b) poemcoloc— LD 信息充分时改用 POEMColoc（面向 top-SNP/不完整数据）
#       (c) skip    — 跳过共定位，方法 A 评分作为主要证据
#   绝不在截断数据上静默输出 colocboost/HyPrColoc 的 PP 值。
#
# 【输入（两种格式，auto-detect）】
#   格式 1 "interval 面板"（推荐，02_analysis/08_GWAS/coloc_stage1 风格）：
#     --panels-dir DIR  每个区间一个 TSV：snp, pos, af, <trait>.beta/.se/.p ...
#     --ld-dir DIR      可选：同名 <interval>.ld.tsv（SNP ID 为行列名的 LD 矩阵）
#   格式 2 "逐性状汇总统计"：
#     --sumstats T1.tsv T2.tsv ...  各含 snp/chr/bp/p[/beta/se]，内部聚类为区间
#
# 【输出】
#   <out-prefix>.coloc_results.tsv   区间×性状集 共定位证据（PP.H4 或 cos_index）
#   <out-prefix>.run_report.tsv      引擎/降级路径/警告记录
#
# 【用法】
#   Rscript run_coloc.R --panels-dir panels/ --ld-dir ld/ \
#       --traits PH,MF,SH --out-prefix output/multitrait/coloc \
#       [--config example_config.yaml] [--mode auto|complete|truncated]
# =============================================================================

suppressPackageStartupMessages({
  ok_yaml <- requireNamespace("yaml", quietly = TRUE)
})

# -----------------------------------------------------------------------------
# 参数解析（轻量，无 optparse 依赖）
# -----------------------------------------------------------------------------
parse_args <- function(argv) {
  res <- list(panels_dir = "", ld_dir = "", sumstats = character(0),
              traits = "", out_prefix = "", config = "", mode = "auto",
              tool = "auto")
  i <- 1
  while (i <= length(argv)) {
    key <- argv[i]
    if (key %in% c("--panels-dir", "--ld-dir", "--traits", "--out-prefix",
                   "--config", "--mode", "--tool")) {
      val <- if (i + 1 <= length(argv)) argv[i + 1] else stop(paste("缺少参数值:", key))
      # 先剥前导 --，再把剩余 - 换成 _（否则 --out-prefix 会变成 __out_prefix）
      nm <- gsub("-", "_", sub("^--", "", key))
      res[[nm]] <- val
      i <- i + 2
    } else if (key == "--sumstats") {
      j <- i + 1
      while (j <= length(argv) && !startsWith(argv[j], "--")) j <- j + 1
      res$sumstats <- argv[(i + 1):(j - 1)]
      i <- j
    } else {
      stop(paste("未知参数:", key))
    }
  }
  res
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (!nzchar(args$out_prefix)) stop("必须指定 --out-prefix")

# -----------------------------------------------------------------------------
# 配置
# -----------------------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (is.character(a) && !nzchar(a[1]))) b else a

cfg <- list(
  p_trunc = 1e-5,
  trunc_tol = 1.5,
  truncation_action = "warn",       # warn | poemcoloc | skip | auto
  pp4_pass = 0.8,
  primary_tool = "colocboost",      # colocboost | hyprcoloc
  susie = FALSE,
  n_per_trait = 212,
  cluster_bp = 500000
)
if (nzchar(args$config)) {
  if (!ok_yaml) stop("配置文件需要 R 包 yaml")
  raw <- yaml::read_yaml(args$config)
  det <- raw$truncation_detection %||% list()
  cfg$p_trunc <- as.numeric(det$p_trunc %||% cfg$p_trunc)
  cfg$trunc_tol <- as.numeric(det$tolerance %||% cfg$trunc_tol)
  cfg$truncation_action <- as.character(raw$method_b$truncation$action %||% cfg$truncation_action)
  cfg$pp4_pass <- as.numeric(raw$method_b$pp4_threshold %||% cfg$pp4_pass)
  cfg$primary_tool <- as.character(raw$method_b$primary %||% cfg$primary_tool)
  cfg$susie <- isTRUE(raw$method_b$susie)
  cfg$n_per_trait <- as.numeric(raw$n_per_trait$total %||% cfg$n_per_trait)
}

mode <- args$mode
action <- cfg$truncation_action
traits <- strsplit(args$traits, ",", fixed = TRUE)[[1]]

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
warn_box <- function(lines) {
  cat("\n", paste(rep("=", 78), collapse = ""), "\n", sep = "")
  for (ln in lines) cat("⚠️  ", ln, "\n", sep = "")
  cat(paste(rep("=", 78), collapse = ""), "\n\n", sep = "")
}

# -----------------------------------------------------------------------------
# 工具函数
# -----------------------------------------------------------------------------
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

read_table_auto <- function(path) {
  sep <- if (grepl("\\.(csv)$", path, ignore.case = TRUE)) "," else "\t"
  utils::read.delim(path, sep = sep, check.names = FALSE, stringsAsFactors = FALSE)
}

# 截断检测：性状 p 列全部 <= p_trunc*tol → 截断
check_truncation <- function(pvec, p_trunc, tol) {
  pvec <- pvec[is.finite(pvec)]
  if (!length(pvec)) return(list(truncated = NA, p_max = NA, n = 0))
  pmax <- max(pvec)
  list(truncated = pmax <= p_trunc * tol, p_max = pmax, n = length(pvec))
}

# -----------------------------------------------------------------------------
# 输入装载
# -----------------------------------------------------------------------------
load_panels <- function(panels_dir) {
  files <- list.files(panels_dir, pattern = "\\.(tsv|csv)$", full.names = TRUE)
  if (!length(files)) stop("panels-dir 无面板文件: ", panels_dir)
  panels <- list()
  for (f in files) {
    id <- sub("\\.(tsv|csv)$", "", basename(f))
    panels[[id]] <- read_table_auto(f)
  }
  panels
}

load_ld <- function(ld_dir, interval_id, snps) {
  f <- file.path(ld_dir, paste0(interval_id, ".ld.tsv"))
  if (!file.exists(f)) return(NULL)
  ld <- as.matrix(read_table_auto(f))
  rn <- ld[, 1]
  ld <- ld[, -1, drop = FALSE]
  storage.mode(ld) <- "numeric"
  if (anyNA(rn) || any(!nzchar(rn))) {
    warning(paste0("[", interval_id, "] LD 文件行名缺失/为空，视为 LD 不充分"))
    return(NULL)
  }
  rownames(ld) <- rn
  colnames(ld) <- rn
  if (!all(snps %in% rn)) return(NULL)   # LD 不覆盖全部 SNP → 不充分
  ld[snps, snps, drop = FALSE]
}

# 面板 → 每性状 dataset 列表（beta/se 或 z；完整 p）
panel_trait_data <- function(panel, trait) {
  beta_col <- paste0(trait, ".beta"); se_col <- paste0(trait, ".se"); p_col <- paste0(trait, ".p")
  if (!all(c(p_col) %in% colnames(panel))) return(NULL)
  keep <- is.finite(panel[[p_col]])
  d <- panel[keep, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  out <- list(snp = as.character(d$snp), position = as.numeric(d$pos),
              pvalues = as.numeric(d[[p_col]]), N = cfg$n_per_trait, type = "quant")
  if (all(c(beta_col, se_col) %in% colnames(d))) {
    out$beta <- as.numeric(d[[beta_col]])
    out$varbeta <- as.numeric(d[[se_col]])^2
  }
  if ("af" %in% colnames(d)) out$MAF <- pmin(as.numeric(d$af), 1 - as.numeric(d$af))
  out
}

# -----------------------------------------------------------------------------
# 共定位引擎
# -----------------------------------------------------------------------------
run_colocboost <- function(datasets, ld) {
  if (!has_pkg("colocboost")) stop("colocboost 未安装")
  sumstat <- lapply(datasets, function(d) {
    df <- data.frame(variant = d$snp, n = d$N, stringsAsFactors = FALSE)
    if (!is.null(d$beta)) {
      df$beta <- d$beta; df$sebeta <- sqrt(d$varbeta)
    } else {
      # 由 p 值换算 |z|（方向未知取正，colocboost 对输入 z 的要求）
      df$z <- qnorm(d$pvalues / 2, lower.tail = FALSE)
    }
    df
  })
  res <- colocboost::colocboost(sumstat = sumstat, LD = ld)
  cos <- res$cos_details$cos
  if (is.null(cos) || !length(cos$cos_index)) {
    return(list(ok = TRUE, coloc = FALSE, detail = "no colocalized signal"))
  }
  idx <- cos$cos_index
  list(ok = TRUE, coloc = TRUE,
       detail = paste0("cos_index=", paste(idx, collapse = ","),
                       "; n_cos_snps=", length(idx)))
}

run_hyprcoloc <- function(datasets) {
  if (!has_pkg("hyprcoloc")) stop("hyprcoloc 未安装")
  snps <- Reduce(intersect, lapply(datasets, function(d) d$snp))
  if (length(snps) < 2) stop("交集 SNP < 2")
  betas <- do.call(cbind, lapply(datasets, function(d) {
    if (!is.null(d$beta)) d$beta[match(snps, d$snp)]
    else qnorm(d$pvalues[match(snps, d$snp)] / 2, lower.tail = FALSE)
  }))
  ses <- do.call(cbind, lapply(datasets, function(d) {
    if (!is.null(d$varbeta)) sqrt(d$varbeta[match(snps, d$snp)]) else rep(1, length(snps))
  }))
  colnames(betas) <- names(datasets); colnames(ses) <- names(datasets)
  res <- hyprcoloc::hyprcoloc(betas, ses, trait.names = names(datasets),
                              snp.id = snps, uniform.priors = TRUE)
  r <- res$results
  if (is.null(r) || !nrow(r)) return(list(ok = TRUE, coloc = FALSE, detail = "no cluster"))
  pp <- as.numeric(r$posterior_prob[1])
  traits_hit <- as.character(r$traits[1])
  cand <- as.character(r$candidate_snp[1])
  list(ok = TRUE, coloc = pp >= cfg$pp4_pass,
       detail = paste0("PP=", round(pp, 4), "; traits=", traits_hit, "; candidate=", cand))
}

run_poemcoloc <- function(datasets, ld, chr_hint = "1") {
  # 真实 API（与 02_analysis/08_GWAS/coloc_stage1/run_supplementary.R 跑通版一致）：
  #   POEMColoc(dA, list(B = dB), gds_file = NULL, R2 = r2, MAF = maf, window_size = 1e6)
  #   - top-SNP 数据集：pvalues 为 lead SNP 标量，s = 0.5 放在 dataset 内
  #   - R2 = r^2 矩阵，rownames/colnames = 位点坐标（非 SNP ID）
  #   - MAF 按坐标命名；lead SNP 用 which.min(p)
  #   - 返回列表，po[[k]]$summary[["PP.H4"]]（或 PP.H4.abf）
  if (!has_pkg("POEMColoc")) stop("POEMColoc 未安装")
  if (is.null(ld)) stop("POEMColoc 需要 LD 矩阵")
  t1 <- datasets[[1]]
  snp2pos <- setNames(as.character(t1$position), as.character(t1$snp))
  if (!all(t1$snp %in% rownames(ld))) stop("LD 未覆盖全部 SNP")
  # R2：按坐标命名（LD 行列名为 SNP ID → 映射为坐标）
  r2 <- (ld[t1$snp, t1$snp, drop = FALSE])^2
  pos_chr <- snp2pos[rownames(r2)]
  dimnames(r2) <- list(pos_chr, pos_chr)
  maf <- if (!is.null(t1$MAF)) setNames(pmin(t1$MAF, 1 - t1$MAF), as.character(t1$position)) else NULL

  mk_top <- function(d) {
    i <- which.min(d$pvalues)
    if (!length(i) || !is.finite(d$pvalues[i])) stop("lead SNP 选择失败")
    list(pos = d$position[i], chr = chr_hint, N = d$N,
         type = "quant", pvalues = d$pvalues[i], s = 0.5)
  }

  best <- NA; note <- character(0)
  dA <- mk_top(t1)
  for (k in 2:length(datasets)) {
    tname <- names(datasets)[k]
    dB <- mk_top(datasets[[k]])
    po <- tryCatch(
      POEMColoc::POEMColoc(dA, setNames(list(dB), tname), gds_file = NULL,
                           R2 = r2, MAF = maf, window_size = 1e6),
      error = function(e) { note <<- c(note, paste0(tname, ":err=", conditionMessage(e))); NULL })
    if (!is.null(po) && length(po)) {
      r1 <- po[[1]]
      pp4 <- tryCatch(as.numeric(r1$summary[["PP.H4"]]), error = function(e) NA)
      if (!length(pp4) || is.na(pp4))
        pp4 <- tryCatch(as.numeric(r1$summary[["PP.H4.abf"]]), error = function(e) NA)
      if (!length(pp4)) pp4 <- NA
      note <- c(note, paste0(tname, ":PP4=", round(pp4, 4)))
      if (is.finite(pp4) && (is.na(best) || pp4 > best)) best <- pp4
    }
  }
  list(ok = TRUE, coloc = is.finite(best) && best >= cfg$pp4_pass,
       detail = paste0("POEMColoc(降级b) max PP4=", round(best, 4),
                       "（top-SNP 模式，s=0.5; ", paste(note, collapse = "; "), "）"))
}

# -----------------------------------------------------------------------------
# 主流程
# -----------------------------------------------------------------------------
log_msg("方法 B 启动: primary=", cfg$primary_tool, " action=", action)

if (nzchar(args$panels_dir)) {
  panels <- load_panels(args$panels_dir)
} else if (length(args$sumstats)) {
  stop("格式 2（逐性状汇总统计）请先用 build_coloc_inputs 风格脚本聚类为区间面板；",
       "本模块直接消费面板目录，避免重复实现聚类逻辑。")
} else {
  stop("必须指定 --panels-dir")
}

results <- list()
reports <- list()
push_report <- function(item, value) {
  reports[[length(reports) + 1]] <<- data.frame(item = item, value = as.character(value))
}

for (id in names(panels)) {
  panel <- panels[[id]]
  snps <- as.character(panel$snp)
  datasets <- list()
  for (t in traits) {
    d <- panel_trait_data(panel, t)
    if (!is.null(d)) datasets[[t]] <- d
  }
  if (length(datasets) < 2) {
    push_report(paste0(id, ".skip"), "有效性状 < 2")
    next
  }

  # ---- 截断检测（对参与的每个性状；lapply 保持按名索引）----
  trunc_flags <- lapply(datasets, function(d) check_truncation(d$pvalues, cfg$p_trunc, cfg$trunc_tol))
  any_trunc <- any(sapply(trunc_flags, function(x) isTRUE(x$truncated)))
  mode_here <- if (mode == "auto") if (any_trunc) "truncated" else "complete" else mode

  ld <- if (nzchar(args$ld_dir)) load_ld(args$ld_dir, id, snps) else NULL

  # ---- 截断数据：不得静默运行 ----
  used_tool <- NA; detail <- NA; coloc <- NA; degraded <- ""
  if (mode_here == "truncated") {
    pmax_str <- paste(sapply(names(trunc_flags), function(t)
      paste0(t, "=", format(trunc_flags[[t]]$p_max, digits = 3))), collapse = ", ")
    base_warn <- c(
      paste0("[区间 ", id, "] 检测到截断汇总统计（p_max: ", pmax_str, "）"),
      "仅含 p<阈值 的显著信号位点，不满足 colocboost/HyPrColoc 的统计前提（需完整 p 分布做先验校准）。",
      "本模块拒绝静默运行贝叶斯共定位。")

    if (action == "poemcoloc" || (action == "auto" && !is.null(ld))) {
      if (is.null(ld)) {
        warn_box(c(base_warn, "降级路径 (b) 选中但 LD 信息不充分 → 回退路径 (a)：分析中止。",
                   "建议：获取完整区域统计数据，或提供覆盖全部区间 SNP 的 LD 矩阵后重跑。"))
        push_report(paste0(id, ".degraded"), "poemcoloc 无充分 LD → (a) 中止")
        next
      }
      warn_box(c(base_warn,
                 "降级路径 (b)：LD 信息充分 → 改用 POEMColoc（面向 top-SNP/不完整数据）。",
                 "结果置信度低于完整数据流程，请在报告中显式标注。"))
      chr_hint <- sub(".*?(chr[0-9A-Za-z_]+).*", "\\1", id)
      if (identical(chr_hint, id)) chr_hint <- "1"
      out <- tryCatch(run_poemcoloc(datasets, ld, chr_hint = chr_hint), error = function(e) list(ok = FALSE, detail = conditionMessage(e)))
      detail <- out$detail
      coloc <- if (isTRUE(out$ok)) out$coloc else NA
      used_tool <- "POEMColoc(降级b)"
      degraded <- "truncated→poemcoloc"
      push_report(paste0(id, ".degraded"), "truncated → (b) POEMColoc")
    } else if (action == "skip") {
      warn_box(c(base_warn,
                 "降级路径 (c)：跳过共定位，采用方法 A（geneGPS+MAGMA）评分作为主要证据。",
                 "本区间不产生 coloc 证据列。"))
      used_tool <- "skipped(降级c)"; detail <- "截断数据→方法 A 为主证据"; coloc <- NA
      degraded <- "truncated→skip"
      push_report(paste0(id, ".degraded"), "truncated → (c) skip")
    } else {
      warn_box(c(base_warn,
                 "降级路径 (a)：明确警告并建议获取完整区域统计数据。",
                 "如数据方能提供未截断的区间级汇总统计（含非显著位点），",
                 "则可恢复 colocboost/HyPrColoc 标准流程。本区间暂不产出 coloc 证据。"))
      used_tool <- "blocked(降级a)"; detail <- "截断数据→需完整区域统计"; coloc <- NA
      degraded <- "truncated→warn_block"
      push_report(paste0(id, ".degraded"), "truncated → (a) warn+block")
    }
  } else {
    # ---- 完整数据：colocboost 主 / HyPrColoc 备 ----
    primary <- cfg$primary_tool
    out <- tryCatch({
      if (primary == "colocboost" && has_pkg("colocboost") && !is.null(ld)) {
        used_tool <<- "colocboost"
        run_colocboost(datasets, ld)
      } else if (primary == "colocboost" && has_pkg("colocboost") && is.null(ld)) {
        push_report(paste0(id, ".note"), "colocboost 需 LD 矩阵，本区间无 LD → 切换 HyPrColoc")
        used_tool <<- "hyprcoloc"
        run_hyprcoloc(datasets)
      } else if (has_pkg("hyprcoloc")) {
        used_tool <<- "hyprcoloc"
        run_hyprcoloc(datasets)
      } else {
        stop("colocboost 与 hyprcoloc 均不可用")
      }
    }, error = function(e) {
      if (has_pkg("hyprcoloc") && !identical(used_tool, "hyprcoloc")) {
        push_report(paste0(id, ".fallback"), paste0("主引擎失败(", conditionMessage(e), ") → HyPrColoc"))
        used_tool <<- "hyprcoloc(fallback)"
        tryCatch(run_hyprcoloc(datasets), error = function(e2) list(ok = FALSE, detail = conditionMessage(e2)))
      } else {
        list(ok = FALSE, detail = conditionMessage(e))
      }
    })
    detail <- out$detail
    coloc <- if (isTRUE(out$ok)) out$coloc else NA
  }

  n_traits_hit <- length(datasets)
  results[[length(results) + 1]] <- data.frame(
    interval = id,
    traits = paste(names(datasets), collapse = "+"),
    n_traits = n_traits_hit,
    tool = used_tool,
    colocalized = coloc,
    evidence = detail,
    degraded_path = degraded,
    data_mode = mode_here,
    confidence_label = if (mode_here == "truncated")
      "截断数据模式，共定位证据为降级结果，不可与完整数据 PP 值等同解释"
    else
      "完整数据模式，贝叶斯共定位可直接解释",
    stringsAsFactors = FALSE)
  push_report(paste0(id, ".mode"), mode_here)
}

# -----------------------------------------------------------------------------
# 输出
# -----------------------------------------------------------------------------
out_dir <- dirname(args$out_prefix)
if (nzchar(out_dir) && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

res_tab <- if (length(results)) do.call(rbind, results) else
  data.frame(interval = character(0), traits = character(0), n_traits = integer(0),
             tool = character(0), colocalized = logical(0), evidence = character(0),
             degraded_path = character(0), data_mode = character(0),
             confidence_label = character(0), stringsAsFactors = FALSE)
utils::write.table(res_tab, paste0(args$out_prefix, ".coloc_results.tsv"),
                   sep = "\t", row.names = FALSE, quote = FALSE)

rep_tab <- if (length(reports)) do.call(rbind, reports) else data.frame(item = character(0), value = character(0))
utils::write.table(rep_tab, paste0(args$out_prefix, ".run_report.tsv"),
                   sep = "\t", row.names = FALSE, quote = FALSE)

log_msg("方法 B 完成: ", nrow(res_tab), " 个区间×性状集记录")
log_msg("  coloc_results: ", args$out_prefix, ".coloc_results.tsv")
if (any(res_tab$data_mode == "truncated")) {
  log_msg("  ⚠ 存在截断模式区间——降级结果见 degraded_path 与 confidence_label 列")
}
