#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Stage 3 — 效应类型分布条形图（Nature 风格）
# 用法: Rscript plot_effect_distribution.R <variants.csv> <out.pdf>

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("用法: Rscript plot_effect_distribution.R <variants.csv> <out.pdf>")
in_csv <- args[1]
out_pdf <- args[2]

df <- read.csv(in_csv, fileEncoding = "UTF-8-BOM")
if (!"EFFECT" %in% colnames(df)) stop("输入 CSV 缺少 EFFECT 列")
if (nrow(df) == 0) {
  message("[WARN] CSV 无记录，生成空图")
  df <- data.frame(EFFECT = character(0), n = integer(0))
}

eff <- as.data.frame(table(df$EFFECT), stringsAsFactors = FALSE)
colnames(eff) <- c("EFFECT", "n")
eff <- eff[order(eff$n, decreasing = TRUE), ]
eff$EFFECT <- factor(eff$EFFECT, levels = rev(eff$EFFECT))

p <- ggplot(eff, aes(x = EFFECT, y = n)) +
  geom_col(fill = "#0072B2", width = 0.7) +
  geom_text(aes(label = n), hjust = -0.2, size = 3, family = "Arial") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Number of variants") +
  theme_classic(base_family = "Arial", base_size = 10) +
  theme(
    axis.line = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks = element_line(linewidth = 0.6, colour = "black"),
    panel.grid = element_blank()
  )

ggsave(out_pdf, p, width = 6, height = 2 + 0.4 * nrow(eff), dpi = 600, device = cairo_pdf)
message(paste("[OK] 输出:", out_pdf))
