#!/usr/bin/env Rscript
# ==============================================================================
# HARD CONTAMINATION FILTER (control-detected genes)
# ==============================================================================
# Cross-species (PDX) RNA-seq: xengsort + Salmon leave residual rat reads on
# conserved human orthologs. The rat-brain Control samples (IL64B, N168B,
# N269B; negligible graft) went through the identical pipeline, so any human-
# gene signal they carry IS contamination.
#
# Rule: drop a gene if its CPM >= CPM_CUTOFF in at least MIN_CONTROLS of the
#       Control samples. Defaults CPM_CUTOFF=1, MIN_CONTROLS=2 (>=2 of 3).
#
# Outputs:
#   - <counts>.decontaminated.tsv  + <lengths>.decontaminated.tsv  (next to the
#     originals; feed these to nf-core/differentialabundance)
#   - results_decontamination/contaminated_genes.tsv   (dropped genes + stats)
#   - results_decontamination/decontamination_sweep.csv (genes removed per cutoff)
#
# The matrices are line-filtered from the originals so their exact values /
# formatting are preserved (only contaminated rows are removed).
#
# Usage:
#   Rscript filter_contaminated_genes.R <counts.tsv> <lengths.tsv> \
#       <metadata.csv> <out_dir> [cpm_cutoff=1] [min_controls=2]
# ==============================================================================

set.seed(12345)
args <- commandArgs(trailingOnly = TRUE)
COUNTS  <- if (length(args) >= 1) args[1] else
    "ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv"
LENGTHS <- if (length(args) >= 2) args[2] else
    "ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_lengths.tsv"
META    <- if (length(args) >= 3) args[3] else "ANALYSIS/metadata_base.csv"
OUTDIR  <- if (length(args) >= 4) args[4] else "ANALYSIS/results_decontamination"
CPM_CUTOFF   <- if (length(args) >= 5) as.numeric(args[5]) else 1
MIN_CONTROLS <- if (length(args) >= 6) as.integer(args[6]) else 2L
SWEEP_CUTOFFS <- c(0.5, 1, 5)

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(..., "\n")

# Read a salmon.merged matrix: gene_id col 1, optional gene_name col 2, then
# per-sample numeric columns. Returns the data.frame + the sample column names.
read_salmon <- function(path) {
    if (!file.exists(path)) stop("File not found: ", path)
    df <- read.delim(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    has_name <- ncol(df) >= 2 &&
        (tolower(colnames(df)[2]) %in% c("gene_name", "gene_symbol", "symbol") ||
         !is.numeric(df[[2]]))
    scols <- if (has_name) seq.int(3, ncol(df)) else seq.int(2, ncol(df))
    list(df = df, ids = as.character(df[[1]]),
         sample_cols = colnames(df)[scols], has_name = has_name)
}

say("==============================================================")
say(" Hard contamination filter")
say("==============================================================")
say("  Counts   :", COUNTS)
say("  Lengths  :", LENGTHS)
say("  Metadata :", META)
say("  Rule     : CPM >=", CPM_CUTOFF, "in >=", MIN_CONTROLS, "controls\n")

ct <- read_salmon(COUNTS)
mat <- as.matrix(ct$df[, ct$sample_cols, drop = FALSE])
rownames(mat) <- ct$ids
storage.mode(mat) <- "double"

meta <- read.csv(META, stringsAsFactors = FALSE)
if (!all(c("sample", "Classification") %in% colnames(meta)))
    stop("Metadata needs 'sample' and 'Classification' columns.")
controls <- intersect(meta$sample[meta$Classification == "Control"], colnames(mat))
tumors   <- intersect(meta$sample[meta$Classification %in% c("Primary", "Recurrent")],
                      colnames(mat))
if (length(controls) < 1) stop("No Control samples present in the matrix.")
say("Controls:", paste(controls, collapse = ", "))
say("Tumors  :", paste(tumors, collapse = ", "), "\n")

# CPM: library size = total counts per sample over all genes in the matrix.
libsize <- colSums(mat)
cpm <- sweep(mat, 2, libsize, "/") * 1e6
ctrl_cpm <- cpm[, controls, drop = FALSE]

# --- sweep: how many genes get flagged at each cutoff / control count ---
sweep_tbl <- do.call(rbind, lapply(SWEEP_CUTOFFS, function(cut) {
    det <- rowSums(ctrl_cpm >= cut)
    data.frame(cpm_cutoff = cut,
               genes_ge1_control = sum(det >= 1),
               genes_ge2_control = sum(det >= 2),
               genes_ge3_control = sum(det >= 3))
}))
say("Removed-gene sweep (genes detected in >= k of", length(controls), "controls):")
print(sweep_tbl)
write.csv(sweep_tbl, file.path(OUTDIR, "decontamination_sweep.csv"), row.names = FALSE)

# --- primary flag ---
det <- rowSums(ctrl_cpm >= CPM_CUTOFF)
contaminated <- det >= MIN_CONTROLS
keep_set <- ct$ids[!contaminated]
say(sprintf("\nPrimary rule -> total %d | drop %d | retain %d genes",
            length(ct$ids), sum(contaminated), length(keep_set)))

# --- dropped-gene report (transparency) ---
report <- data.frame(
    gene_id = ct$ids,
    gene_name = if (ct$has_name) ct$df[[2]] else NA_character_,
    n_controls_detected = det,
    mean_control_cpm = round(rowMeans(ctrl_cpm), 3),
    mean_tumor_cpm = if (length(tumors)) round(rowMeans(cpm[, tumors, drop = FALSE]), 3) else NA_real_,
    stringsAsFactors = FALSE)
for (s in controls) report[[paste0("cpm_", s)]] <- round(ctrl_cpm[, s], 3)
report <- report[contaminated, ]
report <- report[order(-report$mean_control_cpm), ]
write.table(report, file.path(OUTDIR, "contaminated_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
say("Wrote", file.path(OUTDIR, "contaminated_genes.tsv"))

# --- write decontaminated matrices by line-filtering the originals ---
filter_matrix_file <- function(in_path, out_path) {
    lines  <- readLines(in_path)
    header <- lines[1]
    body   <- lines[-1]
    gid    <- sub("\t.*$", "", body)             # gene_id = text before first tab
    kept   <- body[gid %in% keep_set]
    writeLines(c(header, kept), out_path)
    say(sprintf("  wrote %s (%d genes)", out_path, length(kept)))
}
out_counts  <- sub("\\.tsv$", ".decontaminated.tsv", COUNTS)
out_lengths <- sub("\\.tsv$", ".decontaminated.tsv", LENGTHS)
say("\nWriting decontaminated matrices:")
filter_matrix_file(COUNTS, out_counts)
filter_matrix_file(LENGTHS, out_lengths)

say("\nDONE. Point differentialabundance at:")
say("  --matrix                    ", out_counts)
say("  --transcript_length_matrix  ", out_lengths)
