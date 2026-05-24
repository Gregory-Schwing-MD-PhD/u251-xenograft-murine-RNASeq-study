#!/usr/bin/env Rscript
# ==============================================================================
# RUVSeq CONTAMINATION ADJUSTMENT FOR CROSS-SPECIES (PDX) RNA-SEQ
# ==============================================================================
# Study: U251 human glioblastoma xenograft in RNU/RNU rat brain (orthotopic),
#        +/- MRI-guided LITT. Reads were split into Graft (human) / Host (rat)
#        streams by xengsort, then the Graft stream was quantified against the
#        human transcriptome by Salmon.
#
# PROBLEM
#   xengsort cannot fully remove residual rat reads whose sequence is highly
#   conserved with the human ortholog (ribosomal proteins, translation factors,
#   histones, tubulins, ...). Those reads are misclassified as Graft and then
#   misaligned by Salmon onto the corresponding human gene, producing a
#   per-gene, sample-dependent contamination signal that scales with the rat
#   tissue fraction. Recurrent tumours carry a lower graft fraction than Primary
#   tumours, so the contamination differs systematically between the two groups
#   we compare -- exactly the situation that can manufacture a spurious DE/GSEA
#   signal in conserved gene families.
#
# SOLUTION
#   The three Control samples (rat brain only, ~0.3-4.9% graft) went through the
#   identical library prep, sequencing, xengsort and Salmon pipeline. Any reads
#   they place on human genes ARE contamination, sample- and pipeline-matched.
#   We use RUVSeq::RUVs (Risso et al., Nat Biotechnol 2014) to estimate factors
#   of unwanted variation from these Control replicates, then add those factors
#   as covariates in the DESeq2 model for the Primary-vs-Recurrent contrast.
#
# WORKFLOW
#   1. Load Salmon gene counts + full metadata (9 in-vivo samples)
#   2. Pre-filter low-count genes
#   3. Designate the 3 Control samples as the RUVs replicate set (scIdx)
#   4. Estimate W at k = 1, 2, 3 via RUVs (all genes as controls, cIdx)
#   5. DESeq2 with design ~ W_1 + ... + W_k + Classification on Primary vs
#      Recurrent (Controls used only to estimate W, excluded from the contrast)
#   6. Baseline DESeq2 (~ Classification), same Primary-vs-Recurrent contrast
#   7. Concordance vs baseline: Spearman(log2FC), Jaccard(sig sets), top-20 overlap
#   8. PCA before/after adjustment + Primary-vs-Recurrent silhouette per setting
#   9. Recommend k = argmax silhouette over {1,2,3}
#  10. Export adjusted DE, baseline DE, W factors, concordance, PCA PDF, README
#
# USAGE
#   Rscript ruvseq_contamination_adjustment.R \
#       <gene_counts.tsv> <metadata_full.csv> <out_dir>
#   All three are optional; sensible repo-relative defaults are used otherwise.
# ==============================================================================

suppressPackageStartupMessages({
    library(RUVSeq)          # RUVs(); pulls in EDASeq + Biobase
    library(DESeq2)
    library(edgeR)           # RUVSeq dependency (kept explicit)
    library(ggplot2)
    library(dplyr)
    library(EnsDb.Hsapiens.v86)
    library(cluster)         # silhouette()
    library(patchwork)       # multi-panel PCA grid
})

set.seed(12345)

# ------------------------------------------------------------------------------
# Arguments (with repo-relative defaults)
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
COUNTS_TSV <- if (length(args) >= 1) args[1] else
    "ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv"
META_CSV   <- if (length(args) >= 2) args[2] else "ANALYSIS/metadata_full.csv"
OUT_DIR    <- if (length(args) >= 3) args[3] else "ANALYSIS/results_ruvseq"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Analysis constants
PADJ_CUTOFF   <- 0.05
LOG2FC_CUTOFF <- 1.0
K_VALUES      <- c(1, 2, 3)
GROUP_COLORS  <- c("Control"   = "#1f77b4",
                   "Primary"   = "#ff7f0e",
                   "Recurrent" = "#d62728")

cat("==============================================================\n")
cat(" RUVSeq contamination adjustment\n")
cat("==============================================================\n")
cat("  Counts   :", COUNTS_TSV, "\n")
cat("  Metadata :", META_CSV, "\n")
cat("  Out dir  :", OUT_DIR, "\n\n")

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

# Mirror of map_genes_to_symbols() in create_publication_figure_600_dpi.R so the
# symbol mapping is identical to the existing publication pipeline.
map_genes_to_symbols <- function(gene_ids, db = EnsDb.Hsapiens.v86) {
    clean_ids <- sub("\\..*", "", gene_ids)
    if (mean(grepl("^ENSG", clean_ids)) < 0.1) return(clean_ids)
    symbols <- mapIds(db, keys = clean_ids, column = "SYMBOL",
                      keytype = "GENEID", multiVals = "first")
    ifelse(is.na(symbols), clean_ids, symbols)
}

# Remove the linear effect of the RUV factors W from an expression matrix while
# keeping the grand mean per gene. Equivalent to limma::removeBatchEffect with
# `covariates = W`, implemented in base R to avoid an extra dependency.
#   mat : genes x samples (e.g. assay(vst(...)))
#   W   : samples x k (rows aligned to colnames(mat))
residualize_W <- function(mat, W) {
    W <- as.matrix(W[colnames(mat), , drop = FALSE])
    Y <- t(mat)                                   # samples x genes
    X <- cbind(Intercept = 1, W)                  # samples x (k+1)
    beta  <- solve(crossprod(X), crossprod(X, Y)) # (k+1) x genes
    Wbeta <- W %*% beta[-1, , drop = FALSE]       # samples x genes (W effect only)
    t(Y - Wbeta)                                  # genes x samples
}

# PCA panel + Primary-vs-Recurrent silhouette computed from the SAME projection.
pca_and_sil <- function(mat, meta_df, title) {
    pca <- prcomp(t(mat))
    pv  <- 100 * pca$sdev^2 / sum(pca$sdev^2)
    scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
    colnames(scores) <- c("PC1", "PC2")
    scores$Sample <- rownames(scores)
    scores$Classification <- meta_df$Classification[match(scores$Sample, meta_df$sample)]

    sub <- scores[scores$Classification %in% c("Primary", "Recurrent"), ]
    grp <- factor(sub$Classification, levels = c("Primary", "Recurrent"))
    sil <- NA_real_
    if (nlevels(droplevels(grp)) == 2 && nrow(sub) >= 3) {
        d <- dist(sub[, c("PC1", "PC2")])
        s <- cluster::silhouette(as.integer(grp), d)
        sil <- mean(s[, "sil_width"])
    }

    p <- ggplot(scores, aes(x = PC1, y = PC2, color = Classification)) +
        geom_point(size = 3.5) +
        geom_text(aes(label = Sample), vjust = -0.9, size = 2.8, show.legend = FALSE) +
        scale_color_manual(values = GROUP_COLORS) +
        labs(title = title,
             subtitle = sprintf("PC1 %.1f%% | PC2 %.1f%% | Pri-vs-Rec silhouette = %.3f",
                                pv[1], pv[2], sil),
             x = sprintf("PC1 (%.1f%%)", pv[1]),
             y = sprintf("PC2 (%.1f%%)", pv[2])) +
        theme_bw(base_size = 11) +
        theme(legend.position = "bottom",
              plot.title = element_text(face = "bold"),
              plot.subtitle = element_text(size = 9, color = "grey40"))
    list(plot = p, sil = sil)
}

sig_set <- function(df) {
    df$gene_id[!is.na(df$padj) & df$padj < PADJ_CUTOFF &
               abs(df$log2FoldChange) > LOG2FC_CUTOFF]
}
jaccard <- function(a, b) {
    u <- length(union(a, b))
    if (u == 0) return(NA_real_)
    length(intersect(a, b)) / u
}
top_n_ids <- function(df, n = 20) {
    d <- df[!is.na(df$padj), ]
    if (nrow(d) == 0) return(character(0))
    head(d$gene_id[order(d$padj)], n)
}

# Build a results data.frame from a DESeqDataSet for the Recurrent-vs-Primary
# contrast (log2FC > 0 means up in Recurrent), with gene_id + symbol columns.
extract_results <- function(dds) {
    # independentFiltering = TRUE (default alpha) matches differentialabundance's results()
    res <- results(dds, contrast = c("Classification", "Recurrent", "Primary"),
                   independentFiltering = TRUE)
    df  <- as.data.frame(res)
    df$gene_id <- rownames(df)
    df$symbol  <- map_genes_to_symbols(df$gene_id)
    df[, c("gene_id", "symbol", "baseMean", "log2FoldChange",
           "lfcSE", "stat", "pvalue", "padj")]
}

write_de_tsv <- function(df, path) {
    write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("    wrote", path, "(", nrow(df), "genes )\n")
}

# ==============================================================================
# STEP 1: LOAD DATA
# ==============================================================================
cat("[1/10] Loading counts and metadata...\n")

if (!file.exists(COUNTS_TSV)) stop("Counts file not found: ", COUNTS_TSV)
if (!file.exists(META_CSV))   stop("Metadata file not found: ", META_CSV)

raw <- read.delim(COUNTS_TSV, header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE)

gene_ids <- as.character(raw[[1]])
# Detect an optional gene_name / symbol column in position 2.
has_name_col <- ncol(raw) >= 2 &&
    (tolower(colnames(raw)[2]) %in% c("gene_name", "gene_symbol", "symbol") ||
     !is.numeric(raw[[2]]))
if (has_name_col) {
    gene_name_col <- as.character(raw[[2]])
    count_cols    <- seq.int(3, ncol(raw))
    cat("    Detected gene_name column:", colnames(raw)[2], "\n")
} else {
    gene_name_col <- rep(NA_character_, nrow(raw))
    count_cols    <- seq.int(2, ncol(raw))
}
names(gene_name_col) <- gene_ids

counts <- as.matrix(raw[, count_cols, drop = FALSE])
rownames(counts) <- gene_ids
# Salmon/tximport counts can be fractional after scaling -> round to integers.
counts <- round(counts)
storage.mode(counts) <- "integer"
cat("    Loaded matrix:", nrow(counts), "genes x", ncol(counts), "samples\n")

meta <- read.csv(META_CSV, stringsAsFactors = FALSE)
if (!all(c("sample", "Classification") %in% colnames(meta)))
    stop("Metadata must contain 'sample' and 'Classification' columns.")
rownames(meta) <- meta$sample

# Restrict the count matrix to the samples described in the metadata.
missing <- setdiff(meta$sample, colnames(counts))
if (length(missing) > 0)
    stop("Samples in metadata but absent from counts: ", paste(missing, collapse = ", "))
counts <- counts[, meta$sample, drop = FALSE]

meta$Classification <- factor(meta$Classification,
                              levels = c("Control", "Primary", "Recurrent"))
cat("    Group sizes:\n")
print(table(meta$Classification))

ctrl_samples     <- meta$sample[meta$Classification == "Control"]
contrast_samples <- meta$sample[meta$Classification %in% c("Primary", "Recurrent")]
cat("    Controls (W estimation):", paste(ctrl_samples, collapse = ", "), "\n")
cat("    Contrast (Pri vs Rec) :", paste(contrast_samples, collapse = ", "), "\n\n")
if (length(ctrl_samples) < 2) stop("RUVs needs >= 2 Control replicates.")

# ==============================================================================
# STEP 2: PRE-FILTER LOW-COUNT GENES
# ==============================================================================
cat("[2/10] Pre-filtering...\n")
# Two filters, kept deliberately separate:
#  - counts_ruv: a stable expressed set for RUV W-estimation + PCA
#    (count >= 10 in >= 3 of the 9 samples). `counts` stays the full matrix.
#  - the DE step re-filters per contrast to match differentialabundance exactly
#    (count >= 10 in >= 1 contrast sample); see fit_matched_deseq() below.
counts_ruv <- counts[rowSums(counts >= 10) >= 3, , drop = FALSE]
cat("    Genes for RUV/PCA (>=10 in >=3 of 9):", nrow(counts_ruv), "\n")
cat("    (DE gene set filtered separately to match differentialabundance.)\n\n")

# ==============================================================================
# STEP 3: BUILD SeqExpressionSet + RUVs REPLICATE INDEX (scIdx)
# ==============================================================================
cat("[3/10] Building SeqExpressionSet and Control replicate index...\n")

set <- newSeqExpressionSet(
    counts_ruv,
    phenoData = data.frame(Classification = meta$Classification,
                           row.names = colnames(counts_ruv)))
# Upper-quartile between-lane normalisation (canonical RUVSeq pre-step). Used
# only inside RUVs to estimate W; raw counts are passed to DESeq2 separately.
set <- betweenLaneNormalization(set, which = "upper")

# scIdx: one replicate group = the Control samples (their indices into the
# columns of `set`). A single group needs no -1 padding.
ctrl_idx <- match(ctrl_samples, colnames(counts_ruv))
scIdx <- matrix(ctrl_idx, nrow = 1)
cat("    scIdx (Control sample column indices):", paste(ctrl_idx, collapse = ", "), "\n")
cat("    NOTE: a single 3-replicate group spans <= 2 effective factors;\n")
cat("          k=3 is reported for completeness but may be near-degenerate.\n\n")

# ==============================================================================
# STEP 4-5: RUVs AT EACH k + ADJUSTED DESeq2 ON PRIMARY vs RECURRENT
# ==============================================================================
cat("[4/10 + 5/10] Estimating W and running adjusted DESeq2 per k...\n")

# --- DESeq2 parameterization matched to nf-core/differentialabundance v1.5.0 ---
# Verified against the pipeline's deseq_de.R + therapy_v3_params.yaml so the
# RUVSeq-adjusted DE differs from the differentialabundance baseline ONLY by the
# added W covariate(s):
#   filter   : raw count >= 10 in >= 1 contrast sample (filtering_min_abundance=10,
#              filtering_min_samples=1, applied after subset_to_contrast_samples)
#   DESeq()  : Wald test, parametric fit, ratio size factors, betaPrior=FALSE,
#              minReplicatesForReplace=99 (deseq2_min_replicates_for_replace;
#              moot at n=3 but set for exactness)
#   results(): independentFiltering=TRUE (see extract_results)
DA_MIN_ABUNDANCE <- 10
DA_MIN_SAMPLES   <- 1
fit_matched_deseq <- function(count_mat, coldata, design_form) {
    keep <- rowSums(count_mat >= DA_MIN_ABUNDANCE) >= DA_MIN_SAMPLES
    dds <- DESeqDataSetFromMatrix(countData = count_mat[keep, , drop = FALSE],
                                  colData = coldata, design = design_form)
    DESeq(dds, test = "Wald", fitType = "parametric", sfType = "ratio",
          betaPrior = FALSE, minReplicatesForReplace = 99)
}

W_list   <- list()   # full (9-sample) W per k, for PCA + export
de_list  <- list()   # adjusted DE result data.frame per k

run_adjusted_de <- function(k) {
    cat("  --- k =", k, "---\n")
    ruv <- RUVs(set, cIdx = rownames(set), k = k, scIdx = scIdx)
    Wcols <- grep("^W_", colnames(pData(ruv)), value = TRUE)
    W <- as.matrix(pData(ruv)[, Wcols, drop = FALSE])
    colnames(W) <- paste0("W_", seq_len(ncol(W)))
    rownames(W) <- colnames(counts_ruv)
    W_list[[as.character(k)]] <<- W
    cat("    estimated W with", ncol(W), "factor(s)\n")

    # Subset to the contrast samples; W carried as covariates.
    cd <- data.frame(
        Classification = factor(meta[contrast_samples, "Classification"],
                                levels = c("Primary", "Recurrent")),
        W[contrast_samples, , drop = FALSE],
        row.names = contrast_samples,
        check.names = FALSE)
    design_terms <- c(paste0("W_", seq_len(k)), "Classification")
    design_form  <- as.formula(paste("~", paste(design_terms, collapse = " + ")))
    cat("    design:", deparse(design_form), "\n")

    dds <- fit_matched_deseq(counts[, contrast_samples, drop = FALSE], cd, design_form)
    extract_results(dds)
}

for (k in K_VALUES) {
    res_k <- tryCatch(run_adjusted_de(k), error = function(e) {
        cat("    ERROR at k =", k, ":", conditionMessage(e), "\n")
        cat("    Writing NA-filled stub so downstream files still exist.\n")
        data.frame(gene_id = rownames(counts),
                   symbol  = map_genes_to_symbols(rownames(counts)),
                   baseMean = NA_real_, log2FoldChange = NA_real_,
                   lfcSE = NA_real_, stat = NA_real_,
                   pvalue = NA_real_, padj = NA_real_,
                   stringsAsFactors = FALSE)
    })
    de_list[[as.character(k)]] <- res_k
    write_de_tsv(res_k, file.path(OUT_DIR, sprintf("ruvseq_adjusted_de_k%d.tsv", k)))
}
cat("\n")

# ==============================================================================
# STEP 6: BASELINE DESeq2 (NO RUV ADJUSTMENT)
# ==============================================================================
cat("[6/10] Baseline DESeq2 (~ Classification), Primary vs Recurrent...\n")
cd0 <- data.frame(
    Classification = factor(meta[contrast_samples, "Classification"],
                            levels = c("Primary", "Recurrent")),
    row.names = contrast_samples)
dds0 <- fit_matched_deseq(counts[, contrast_samples, drop = FALSE], cd0, ~ Classification)
res_baseline <- extract_results(dds0)
write_de_tsv(res_baseline, file.path(OUT_DIR, "ruvseq_baseline_de.tsv"))
cat("\n")

# ==============================================================================
# STEP 7: CONCORDANCE METRICS (each k vs baseline)
# ==============================================================================
cat("[7/10] Computing concordance metrics vs baseline...\n")

base_sig   <- sig_set(res_baseline)
base_top20 <- top_n_ids(res_baseline, 20)

concordance_row <- function(label, df) {
    s <- sig_set(df)
    if (label == "baseline") {
        sp  <- 1
        jac <- 1
        ovl <- length(base_top20)
    } else {
        m <- merge(res_baseline[, c("gene_id", "log2FoldChange")],
                   df[, c("gene_id", "log2FoldChange")],
                   by = "gene_id", suffixes = c("_base", "_k"))
        sp  <- suppressWarnings(cor(m$log2FoldChange_base, m$log2FoldChange_k,
                                    method = "spearman", use = "complete.obs"))
        jac <- jaccard(base_sig, s)
        ovl <- length(intersect(base_top20, top_n_ids(df, 20)))
    }
    data.frame(
        setting                  = label,
        n_genes_tested           = sum(!is.na(df$padj)),
        n_sig                    = length(s),
        n_sig_up                 = sum(!is.na(df$padj) & df$padj < PADJ_CUTOFF &
                                       df$log2FoldChange >  LOG2FC_CUTOFF),
        n_sig_down               = sum(!is.na(df$padj) & df$padj < PADJ_CUTOFF &
                                       df$log2FoldChange < -LOG2FC_CUTOFF),
        spearman_lfc_vs_baseline = round(sp, 4),
        jaccard_sig_vs_baseline  = round(jac, 4),
        top20_overlap_vs_baseline = ovl,
        stringsAsFactors = FALSE)
}

concordance <- concordance_row("baseline", res_baseline)
for (k in K_VALUES) {
    concordance <- rbind(concordance,
                         concordance_row(paste0("k", k), de_list[[as.character(k)]]))
}
print(concordance)
cat("\n")

# ==============================================================================
# STEP 8: PCA BEFORE/AFTER + SILHOUETTE-BASED k SELECTION
# ==============================================================================
cat("[8/10] PCA before/after adjustment + silhouettes...\n")

# VST on all 9 samples (blind, unsupervised) for visualisation. vst() alone does
# NOT subtract covariates, so the "after" panels apply residualize_W() to remove
# the estimated W effect -- the scientifically correct way to visualise RUV
# adjustment (cf. limma::removeBatchEffect).
dds_all <- DESeqDataSetFromMatrix(
    countData = counts_ruv,
    colData   = data.frame(Classification = meta$Classification,
                           row.names = colnames(counts_ruv)),
    design    = ~ Classification)
vsd_all <- vst(dds_all, blind = TRUE)
mat_all <- assay(vsd_all)

panel_base <- pca_and_sil(mat_all, meta, "Baseline (no adjustment)")
panels <- list(panel_base$plot)
sil_vec <- c(baseline = panel_base$sil)

for (k in K_VALUES) {
    W <- W_list[[as.character(k)]]
    if (is.null(W)) {            # RUVs failed for this k -> reuse baseline view
        pk <- pca_and_sil(mat_all, meta, sprintf("RUVs k=%d (W unavailable)", k))
    } else {
        mat_adj <- residualize_W(mat_all, W)
        pk <- pca_and_sil(mat_adj, meta, sprintf("RUVs adjusted (k=%d)", k))
    }
    panels[[length(panels) + 1]] <- pk$plot
    sil_vec[paste0("k", k)] <- pk$sil
}

cat("    Primary-vs-Recurrent silhouette by setting:\n")
print(round(sil_vec, 4))

pca_grid <- wrap_plots(panels, ncol = 2) +
    plot_annotation(
        title = "PCA before/after RUVSeq contamination adjustment",
        subtitle = "Coloured by Classification; W effect residualised in adjusted panels")
pca_pdf <- file.path(OUT_DIR, "ruvseq_pca_before_after.pdf")
ggsave(pca_pdf, pca_grid, width = 11, height = 10, device = "pdf")
cat("    wrote", pca_pdf, "\n\n")

# ==============================================================================
# STEP 9: RECOMMENDED k (max silhouette over k = 1,2,3)
# ==============================================================================
cat("[9/10] Selecting recommended k...\n")
k_sils <- sil_vec[paste0("k", K_VALUES)]
k_sils[is.na(k_sils)] <- -Inf            # never recommend a failed/NA setting
recommended_k <- K_VALUES[which.max(k_sils)]
cat(sprintf("    silhouettes: k1=%.3f  k2=%.3f  k3=%.3f  (baseline=%.3f)\n",
            sil_vec["k1"], sil_vec["k2"], sil_vec["k3"], sil_vec["baseline"]))
cat("    Recommended k =", recommended_k, "\n")
if (!is.na(sil_vec["baseline"]) && sil_vec["baseline"] >= max(k_sils))
    cat("    NB: baseline silhouette >= all k; adjustment did not improve\n",
        "       Primary/Recurrent separation. Inspect the PCA before adopting.\n")
writeLines(as.character(recommended_k),
           file.path(OUT_DIR, "ruvseq_recommended_k.txt"))
cat("\n")

# ==============================================================================
# STEP 10: EXPORT W FACTORS + CONCORDANCE SUMMARY
# ==============================================================================
cat("[10/10] Writing W factors and concordance summary...\n")

W_long <- do.call(rbind, lapply(names(W_list), function(kk) {
    W <- W_list[[kk]]
    do.call(rbind, lapply(colnames(W), function(fac) {
        data.frame(sample = rownames(W),
                   Classification = meta[rownames(W), "Classification"],
                   k_setting = paste0("k", kk),
                   factor = fac,
                   value = W[, fac],
                   row.names = NULL,
                   stringsAsFactors = FALSE)
    }))
}))
write.csv(W_long, file.path(OUT_DIR, "ruvseq_estimated_W_factors.csv"),
          row.names = FALSE)
cat("    wrote", file.path(OUT_DIR, "ruvseq_estimated_W_factors.csv"), "\n")

# Append the silhouette to the concordance table for convenience.
concordance$pri_rec_silhouette <- round(sil_vec[concordance$setting], 4)
concordance$recommended <- ifelse(concordance$setting == paste0("k", recommended_k),
                                  "*", "")
write.csv(concordance, file.path(OUT_DIR, "ruvseq_concordance_summary.csv"),
          row.names = FALSE)
cat("    wrote", file.path(OUT_DIR, "ruvseq_concordance_summary.csv"), "\n\n")

# ------------------------------------------------------------------------------
# TOP-LINE SUMMARY (visible in the SLURM log)
# ------------------------------------------------------------------------------
rec_row <- concordance[concordance$setting == paste0("k", recommended_k), ]
cat("==============================================================\n")
cat(" SUMMARY\n")
cat("==============================================================\n")
cat("  Recommended k                         :", recommended_k, "\n")
cat("  Spearman(log2FC) recommended-k vs base :",
    rec_row$spearman_lfc_vs_baseline, "\n")
cat("  Jaccard(sig set) recommended-k vs base :",
    rec_row$jaccard_sig_vs_baseline, "\n")
cat("  Baseline significant genes             :",
    concordance$n_sig[concordance$setting == "baseline"], "\n")
cat("  Recommended-k significant genes        :", rec_row$n_sig, "\n")
cat("  Outputs written to                     :", OUT_DIR, "\n")
cat("==============================================================\n")
cat("DONE.\n")
