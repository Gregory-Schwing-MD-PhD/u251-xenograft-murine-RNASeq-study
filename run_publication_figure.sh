#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --job-name=pub_figure
#SBATCH --output=publication_figure_%j.out
#SBATCH --error=publication_figure_%j.err

set -euo pipefail

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║      PUBLICATION FIGURE GENERATOR - 9 PANEL COMPOSITE          ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

# ==============================================================================
# SETUP
# ==============================================================================
export CONDA_PREFIX="${HOME}/mambaforge/envs/nextflow"
export PATH="${CONDA_PREFIX}/bin:$PATH"
unset JAVA_HOME

export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
mkdir -p "$XDG_RUNTIME_DIR" "$NXF_SINGULARITY_CACHEDIR"

export NXF_SINGULARITY_HOME_MOUNT=true
unset LD_LIBRARY_PATH PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE

# ==============================================================================
# PATH CONFIGURATION
# ==============================================================================
GMT_DIR="ANALYSIS/refs/pathways/human"
COMBINED_GMT="${GMT_DIR}/combined_human.gmt"

STRING_DIR="ANALYSIS/refs/pathways/human/human_string"
RESULTS_DIR="ANALYSIS/results_therapy_v3"
NORM_FILE="${RESULTS_DIR}/tables/processed_abundance/all.rlog.tsv"
METADATA_FILE="ANALYSIS/metadata_therapy.csv"
VISUAL_ABSTRACT_SRC="ASSETS/Visual_abstract.png"
VISUAL_ABSTRACT_LINK="Experiment_Visual_Abstract.png"
TARGET_CONTRAST="therapy_impact"

# Optional positional args, for comparing decontamination approaches:
#   $1 = DE results table to use INSTEAD of the differentialabundance output
#        (repo-relative; e.g. ANALYSIS/results_ruvseq/ruvseq_adjusted_de_k2.tsv)
#   $2 = label suffix for the figure output dir (e.g. ruvseq -> publication_figure_ruvseq)
DE_RESULTS_OVERRIDE="${1:-}"
FIG_LABEL="${2:-}"
OUT_DIR="publication_figure${FIG_LABEL:+_$FIG_LABEL}"
DE_ARG=""
[ -n "$DE_RESULTS_OVERRIDE" ] && DE_ARG="/data/$DE_RESULTS_OVERRIDE"

echo "Verifying input files..."
for path in "$NORM_FILE" "$COMBINED_GMT" "$STRING_DIR" "$METADATA_FILE" "$VISUAL_ABSTRACT_SRC"; do
    if [ ! -e "$path" ]; then
        echo "ERROR: Required path not found: $path"
        exit 1
    fi
done

# Symlink visual abstract into cwd under the name the R script expects
ln -sf "$VISUAL_ABSTRACT_SRC" "$VISUAL_ABSTRACT_LINK"

echo "✓ All required files found"
echo ""

mkdir -p "$OUT_DIR"

# ==============================================================================
# CONTAINER SETUP
# ==============================================================================
IMG_PATH="${NXF_SINGULARITY_CACHEDIR}/go2432-bioconductor.sif"

if [[ ! -f "$IMG_PATH" ]]; then
    echo "Pulling Bioconductor container..."
    singularity pull "$IMG_PATH" docker://go2432/bioconductor:latest
fi

echo "✓ Container ready: $IMG_PATH"
echo ""

# ==============================================================================
# RUN PUBLICATION FIGURE GENERATION
# ==============================================================================
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║      GENERATING PUBLICATION FIGURE (Panels A-I)                ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "Configuration:"
echo "  Norm matrix:      $NORM_FILE"
echo "  Metadata File:    $METADATA_FILE"
echo "  Results Dir:      $RESULTS_DIR"
echo "  GMT Dir:          $GMT_DIR"
echo "  Combined GMT:     $COMBINED_GMT"
echo "  STRING Dir:       $STRING_DIR"
echo "  Visual Abstract:  $VISUAL_ABSTRACT_SRC -> $VISUAL_ABSTRACT_LINK"
echo "  Output Dir:       $OUT_DIR"
echo "  Target Contrast:  $TARGET_CONTRAST"
echo "  DE override:      ${DE_RESULTS_OVERRIDE:-<none (use differentialabundance output)>}"
echo ""

singularity exec --bind "$PWD:/data" --pwd /data "$IMG_PATH" \
    Rscript /data/create_publication_figure_600_dpi.R \
    "/data/$NORM_FILE" \
    "/data/$RESULTS_DIR" \
    "/data/$GMT_DIR" \
    "/data/$STRING_DIR" \
    "/data/$OUT_DIR" \
    "$TARGET_CONTRAST" \
    "/data/$METADATA_FILE" \
    "$DE_ARG"

exit_code=$?

if [ $exit_code -ne 0 ]; then
    echo ""
    echo "ERROR: Publication figure generation failed (exit code: $exit_code)"
    exit $exit_code
fi

# ==============================================================================
# VERIFY OUTPUTS
# ==============================================================================
echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║      VERIFYING OUTPUTS                                         ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

required_outputs=(
    "$OUT_DIR/Publication_Figure_9Panel_VOLCANO_COMPLETE.png"
    "$OUT_DIR/Publication_Figure_9Panel_VOLCANO_COMPLETE.pdf"
    "$OUT_DIR/Figure_Caption.txt"
)

all_present=true
for output in "${required_outputs[@]}"; do
    if [ -f "$output" ]; then
        size=$(du -h "$output" | cut -f1)
        echo "✓ $output ($size)"
    else
        echo "✗ MISSING: $output"
        all_present=false
    fi
done

# Optional machine-readable drug-discovery reports (only present when drug
# discovery ran for this contrast). Reported but not required -> no failure.
echo ""
echo "Drug-discovery reports (optional):"
for opt in \
    "$OUT_DIR/${TARGET_CONTRAST}_LLM_Analysis_Report.txt" \
    "$OUT_DIR/${TARGET_CONTRAST}_Drug_Profiles_Comprehensive.csv" \
    "$OUT_DIR/${TARGET_CONTRAST}_Analysis_Report.html"; do
    if [ -f "$opt" ]; then echo "✓ $opt"; else echo "○ not generated: $opt"; fi
done

if [ "$all_present" = false ]; then
    echo ""
    echo "WARNING: Some expected outputs are missing"
    exit 1
fi

# ==============================================================================
# SUMMARY
# ==============================================================================
echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║      PUBLICATION FIGURE GENERATION COMPLETE                    ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "📊 GENERATED OUTPUTS:"
echo ""
echo "Main Figure:"
echo "  • PNG (high-res): $OUT_DIR/Publication_Figure_9Panel_VOLCANO_COMPLETE.png"
echo "  • PDF (vector):   $OUT_DIR/Publication_Figure_9Panel_VOLCANO_COMPLETE.pdf"
echo "  • Captions:       $OUT_DIR/Figure_Caption.txt"
echo "  • LLM report:     $OUT_DIR/${TARGET_CONTRAST}_LLM_Analysis_Report.txt  (all panels A-I)"
echo "  • Drug CSV:       $OUT_DIR/${TARGET_CONTRAST}_Drug_Profiles_Comprehensive.csv"
echo ""
echo "Panel Layout (3×3 grid):"
echo "  A. Experimental Design (Gemini visual abstract)"
echo "  B. Global Structure (PCA biplot + scree plot)"
echo "  C. Subtype Trajectories (with significance markers)"
echo "  D. Semantic Pathway Clustering (tree plot)"
echo "  E. Protein-Protein Interaction Network"
echo "  F. Polypharmacology Network"
echo "  G. Drug BBB Penetration Scores"
echo "  H. Drug-Pathway Gene Overlap Heatmap"
echo "  I. Top 5 Drug Candidates (integrated scoring)"
echo ""
echo "📝 PANEL DESCRIPTIONS:"
echo ""
cat "$OUT_DIR/Figure_Caption.txt"
echo ""
echo "════════════════════════════════════════════════════════════════"
echo "                    ANALYSIS COMPLETE"
echo "════════════════════════════════════════════════════════════════"
