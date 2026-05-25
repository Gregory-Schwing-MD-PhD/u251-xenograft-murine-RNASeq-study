#!/usr/bin/env Rscript

# ==========================================
# MICROENVIRONMENT VISUALIZATION SUITE (RAT - ENSEMBL UPDATE)
# ==========================================
# Features:
# 1. Maps Rat Ensembl IDs -> Gene Symbols (using EnsDb.Rnorvegicus.v79)
# 2. Filters for "Clean" genes (Signal Subtraction)
# 3. Generates Validated Volcano Plots & Heatmaps

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ComplexHeatmap)
    library(circlize)
    library(EnhancedVolcano)
    library(clusterProfiler)
    library(enrichplot)
    library(stringr)
    # library(org.Rn.eg.db) # <-- REMOVED
    library(EnsDb.Rnorvegicus.v79) # <-- ADDED (Native Ensembl Rat DB)
})

# --- ARGUMENTS ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
    stop("Usage: Rscript plot_kitchen_sink_microenvironment.R <DESeq2_Res> <VST> <GMT> <OutPrefix> <Counts> <CleanGenesCSV>")
}

deseq_file   <- args[1]
vst_file     <- args[2]
gmt_file     <- args[3]
out_prefix   <- args[4]
counts_file  <- args[5]
clean_genes_file <- args[6]

# Create output directory
out_dir <- dirname(out_prefix)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ==========================================
# 1. LOAD DATA & MAP IDs
# ==========================================
cat("LOG: Loading Data...\n")
res_df <- fread(deseq_file)
vst_df <- fread(vst_file)
clean_genes <- fread(clean_genes_file)

# --- ID MAPPING FUNCTION (Rat Ensembl) ---
map_rat_ids <- function(ids) {
    # Remove version numbers if present (ENSRNOG...1 -> ENSRNOG...)
    clean_ids <- sub("\\..*", "", ids)
    
    # Map using EnsDb (Keytype is GENEID)
    syms <- suppressWarnings(mapIds(EnsDb.Rnorvegicus.v79, 
                                    keys=clean_ids, 
                                    column="SYMBOL", 
                                    keytype="GENEID", 
                                    multiVals="first"))
                                    
    syms[is.na(syms)] <- clean_ids[is.na(syms)] # Fallback to ID
    return(make.unique(as.character(syms)))
}

# Apply Mapping to Results
# Check if gene_id column exists, otherwise use first column
if (!"gene_id" %in% colnames(res_df)) colnames(res_df)[1] <- "gene_id"
res_df$symbol <- map_rat_ids(res_df$gene_id)

# Apply Mapping to VST
colnames(vst_df)[1] <- "gene_id"
vst_mat <- as.matrix(vst_df[, -1, with=FALSE])
rownames(vst_mat) <- map_rat_ids(vst_df$gene_id)

# ==========================================
# 2. FILTERING (Signal Subtraction)
# ==========================================
# The clean_genes file likely contains IDs. We match on IDs.
valid_ids <- clean_genes$Gene
cat(sprintf("LOG: Filtering to %d validated LITT-response genes.\n", length(valid_ids)))

# Subset VST for Heatmaps (Only validated genes)
# We need to find which VST rows correspond to valid IDs.
# Since we renamed VST rows to Symbols, we need a look-up or re-map.
# Safer strategy: Filter VST by index based on the original ID column
valid_indices <- which(vst_df$gene_id %in% valid_ids)
vst_clean <- vst_mat[valid_indices, , drop=FALSE]

# ==========================================
# 3. ENHANCED VOLCANO (Validated Highlight)
# ==========================================
cat("LOG: Generating Volcano Plot...\n")

# Logic: Highlight Validated genes in Red/Blue, Artifacts in Grey
keyvals <- ifelse(
    res_df$gene_id %in% valid_ids & res_df$log2FoldChange > 0, 'red',
    ifelse(res_df$gene_id %in% valid_ids & res_df$log2FoldChange < 0, 'blue',
    'grey'))

names(keyvals)[keyvals == 'red']  <- 'Validated Upregulated'
names(keyvals)[keyvals == 'blue'] <- 'Validated Downregulated'
names(keyvals)[keyvals == 'grey'] <- 'Artifact/Background'

pdf(paste0(out_prefix, "_Volcano_Validated.pdf"), width=10, height=8)
EnhancedVolcano(res_df,
    lab = res_df$symbol, # Now using Symbols!
    x = 'log2FoldChange',
    y = 'padj',
    title = 'Rat Microenvironment Response',
    subtitle = 'Signal Subtraction Validation',
    pCutoff = 0.05,
    FCcutoff = 1.0,
    colCustom = keyvals,
    legendPosition = 'right',
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    max.overlaps = 20
)
dev.off()

# ==========================================
# 4. HEATMAP (Validated Genes Only)
# ==========================================
cat("LOG: Generating Heatmap...\n")

if (nrow(vst_clean) > 2) {
    mat_scaled <- t(scale(t(vst_clean)))
    col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

    # Dynamic height
    h_calc <- min(30, max(6, nrow(vst_clean)/5))

    pdf(paste0(out_prefix, "_Heatmap_CleanGenes.pdf"), width=10, height=h_calc)
    ht <- Heatmap(mat_scaled,
        name = "Z-score",
        col = col_fun,
        show_row_names = nrow(vst_clean) < 80, # Show names if list is manageable
        show_column_names = TRUE,
        cluster_columns = TRUE,
        column_title = "Validated Rat Response"
    )
    draw(ht)
    dev.off()
} else {
    cat("WARNING: Not enough valid genes for heatmap.\n")
}

# ==========================================
# 5. PATHWAY ANALYSIS (ORA)
# ==========================================
cat("LOG: Running Pathway Analysis (ORA)...\n")

# Use Symbols for Pathway Analysis (easier if GMT is symbols)
# Assuming your Rat GMT uses Symbols. If it uses IDs, change 'universe' to res_df$gene_id
pathways <- read.gmt(gmt_file)
universe <- res_df$symbol

# Get Symbols of Validated Up/Down genes
sig_up_ids <- clean_genes[clean_genes$LFC_LITT > 0]$Gene
sig_down_ids <- clean_genes[clean_genes$LFC_LITT < 0]$Gene

# Map these IDs to Symbols for the enrichment test
sig_up_sym <- res_df$symbol[res_df$gene_id %in% sig_up_ids]
sig_down_sym <- res_df$symbol[res_df$gene_id %in% sig_down_ids]

run_ora <- function(genes, title_suffix) {
    if (length(genes) > 5) { # Need minimal gene set size
        ego <- enricher(genes,
            TERM2GENE = pathways,
            universe = universe,
            pvalueCutoff = 0.1,
            qvalueCutoff = 0.2
        )

        if (!is.null(ego) && nrow(ego@result) > 0) {
            pdf(paste0(out_prefix, "_ORA_Dotplot_", title_suffix, ".pdf"), width=10, height=8)
            print(dotplot(ego, showCategory=20) + ggtitle(paste("Pathways:", title_suffix)))
            dev.off()
            write.csv(ego@result, paste0(out_prefix, "_ORA_Table_", title_suffix, ".csv"))
        }
    }
}

run_ora(sig_up_sym, "Upregulated_Clean")
run_ora(sig_down_sym, "Downregulated_Clean")

cat("LOG: Analysis Complete.\n")
