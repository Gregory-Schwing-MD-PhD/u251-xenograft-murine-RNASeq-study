#!/usr/bin/env Rscript
# ==============================================================================
# IS CONTAMINATION ACTUALLY DISTORTING THE DE? — a direct test
# ==============================================================================
# Contamination (residual rat reads on conserved human genes) scales with the
# host/rat fraction, and Recurrent tumours carry more host than Primary. So IF
# contamination is biasing the Primary-vs-Recurrent DE, the genes that catch rat
# reads (i.e. are "expressed" in the rat-brain Control samples) should show a
# SYSTEMATICALLY different fold-change than clean genes.
#
# This script tests that prediction on the (un-decontaminated) baseline DE:
#   1. Of the baseline significant genes, what fraction are control-detected?
#      (Are your actual hits contaminated?)
#   2. Spearman(mean control CPM, baseline log2FC) across all tested genes.
#      (Is there a global contamination-direction skew?)
#   3. Wilcoxon of log2FC, control-detected vs clean genes. (Distribution shift?)
# Outputs a printed verdict + a 2-panel diagnostic PDF.
#
# Usage:
#   Rscript test_contamination_bias.R <counts.tsv> <metadata.csv> \
#       <baseline_de.tsv> <out_dir>
# ==============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
set.seed(12345)

args <- commandArgs(trailingOnly = TRUE)
COUNTS   <- if (length(args) >= 1) args[1] else
    "ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv"
META     <- if (length(args) >= 2) args[2] else "ANALYSIS/metadata_base.csv"
BASE_DE  <- if (length(args) >= 3) args[3] else
    "ANALYSIS/results_ruvseq/ruvseq_baseline_de.tsv"
OUT_DIR  <- if (length(args) >= 4) args[4] else "ANALYSIS/results_decontamination"
PADJ <- 0.05; LFC <- 1
CPM_CUT <- 1; MIN_CTRL <- 2     # "control-detected" rule (same as the hard filter)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(..., "\n")

# --- load counts, derive control CPM ---
raw <- read.delim(COUNTS, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
ids <- as.character(raw[[1]])
has_name <- ncol(raw) >= 2 && (tolower(colnames(raw)[2]) %in% c("gene_name","gene_symbol","symbol") || !is.numeric(raw[[2]]))
scols <- if (has_name) seq.int(3, ncol(raw)) else seq.int(2, ncol(raw))
mat <- as.matrix(raw[, scols, drop = FALSE]); rownames(mat) <- ids; storage.mode(mat) <- "double"

meta <- read.csv(META, stringsAsFactors = FALSE)
controls <- intersect(meta$sample[meta$Classification == "Control"], colnames(mat))
say("Controls used as contamination proxy:", paste(controls, collapse = ", "))

cpm <- sweep(mat, 2, colSums(mat), "/") * 1e6
ctrl_cpm <- cpm[, controls, drop = FALSE]
gene_tbl <- data.frame(
    gene_id = ids,
    mean_ctrl_cpm = rowMeans(ctrl_cpm),
    ctrl_detected = rowSums(ctrl_cpm >= CPM_CUT) >= MIN_CTRL,
    stringsAsFactors = FALSE)

# --- baseline DE ---
de <- read.delim(BASE_DE, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (!"gene_id" %in% colnames(de)) de$gene_id <- as.character(de[[1]])
de <- de[, c("gene_id", "log2FoldChange", "padj")]
df <- merge(de, gene_tbl, by = "gene_id")
df <- df[is.finite(df$log2FoldChange) & !is.na(df$padj), ]
df$sig <- df$padj < PADJ & abs(df$log2FoldChange) > LFC
say(sprintf("Genes tested: %d | baseline significant: %d", nrow(df), sum(df$sig)))

# ============================ THE THREE CHECKS ============================
# 1. Are the significant hits contamination-prone?
sig <- df[df$sig, ]
pct_sig_contam <- if (nrow(sig)) 100 * mean(sig$ctrl_detected) else NA
pct_all_contam <- 100 * mean(df$ctrl_detected)
say(sprintf("\n[1] Control-detected genes: %.1f%% of ALL tested, %.1f%% of SIGNIFICANT",
            pct_all_contam, pct_sig_contam))
say(sprintf("    (enrichment of contamination among hits: %.2fx)",
            if (is.finite(pct_sig_contam) && pct_all_contam > 0) pct_sig_contam / pct_all_contam else NA))

# 2. Global skew: does contamination level track fold-change?
rho <- suppressWarnings(cor(log10(df$mean_ctrl_cpm + 1), df$log2FoldChange,
                            method = "spearman", use = "complete.obs"))
say(sprintf("\n[2] Spearman(control CPM, log2FC) = %.3f", rho))

# 3. Distribution shift: log2FC of contaminated vs clean genes
wt <- suppressWarnings(wilcox.test(log2FoldChange ~ ctrl_detected, data = df))
med_clean  <- median(df$log2FoldChange[!df$ctrl_detected])
med_contam <- median(df$log2FoldChange[df$ctrl_detected])
say(sprintf("\n[3] median log2FC  clean = %+.3f   contaminated = %+.3f   (Wilcoxon p = %.2g)",
            med_clean, med_contam, wt$p.value))

# ----------------------------- verdict -----------------------------
flags <- c(abs(rho) >= 0.15,
           is.finite(pct_sig_contam) && pct_all_contam > 0 && pct_sig_contam / pct_all_contam >= 1.5,
           abs(med_contam - med_clean) >= 0.25 && wt$p.value < 0.05)
say("\n==============================================================")
if (sum(flags) == 0) {
    say(" VERDICT: little evidence contamination is biasing the DE.")
    say("  -> baseline (post-xengsort) is defensible; gene filtering optional.")
} else if (sum(flags) >= 2) {
    say(" VERDICT: contamination DOES appear to skew the DE (",
        sum(flags), "of 3 flags).")
    say("  -> a gentle (ratio) filter or covariate adjustment is warranted.")
} else {
    say(" VERDICT: weak/mixed evidence (", sum(flags), "of 3 flags) — judgement call.")
}
say("==============================================================")

# ----------------------------- plots -----------------------------
df$group <- ifelse(df$ctrl_detected, "control-detected", "clean")
p1 <- ggplot(df, aes(log10(mean_ctrl_cpm + 1), log2FoldChange, color = sig)) +
    geom_point(alpha = 0.3, size = 0.6) +
    geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#d62728"),
                       labels = c("ns", "sig"), name = NULL) +
    labs(title = sprintf("log2FC vs contamination level (Spearman %.2f)", rho),
         x = "log10(mean control CPM + 1)", y = "baseline log2FC (Rec vs Pri)") +
    theme_bw(base_size = 11)
p2 <- ggplot(df, aes(group, log2FoldChange, fill = group)) +
    geom_violin(alpha = 0.5) + geom_boxplot(width = 0.15, outlier.size = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = c("clean" = "#1f77b4", "control-detected" = "#ff7f0e"), guide = "none") +
    labs(title = sprintf("Wilcoxon p = %.2g", wt$p.value),
         x = NULL, y = "baseline log2FC (Rec vs Pri)") +
    theme_bw(base_size = 11)
pdf_path <- file.path(OUT_DIR, "contamination_bias_test.pdf")
ggsave(pdf_path, p1 + p2 + plot_annotation(
    title = "Is contamination biasing the Primary-vs-Recurrent DE?"),
    width = 11, height = 5)
say("\nWrote", pdf_path)
write.csv(df[order(-df$mean_ctrl_cpm), ],
          file.path(OUT_DIR, "contamination_bias_per_gene.csv"), row.names = FALSE)
say("Wrote", file.path(OUT_DIR, "contamination_bias_per_gene.csv"))
