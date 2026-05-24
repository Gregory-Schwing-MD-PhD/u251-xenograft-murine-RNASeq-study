#!/bin/bash
# ==============================================================================
# HARD CONTAMINATION FILTER - SLURM SUBMISSION (Wayne State Warrior HPC)
# ==============================================================================
# Drops human genes detected in the rat-brain Control samples and writes
# decontaminated count + length matrices for nf-core/differentialabundance.
# Fast (matrix filtering); mirrors run_ruvseq.sh's container invocation.
# ==============================================================================
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --job-name=decontam
#SBATCH --mail-user=go2432@wayne.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=slurm_logs/decontam_%j.out
#SBATCH --error=slurm_logs/decontam_%j.err

set -euo pipefail

echo "================================================================"
echo " Hard contamination filter"
echo " Job ID : ${SLURM_JOB_ID:-local}   Node: $(hostname)"
echo " Start  : $(date)"
echo "================================================================"

# Environment (matches run_ruvseq.sh / run_publication_figure.sh)
export CONDA_PREFIX="${HOME}/mambaforge/envs/nextflow"
export PATH="${CONDA_PREFIX}/bin:$PATH"
unset JAVA_HOME
export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
mkdir -p "$XDG_RUNTIME_DIR" "$NXF_SINGULARITY_CACHEDIR"
export NXF_SINGULARITY_HOME_MOUNT=true
unset LD_LIBRARY_PATH PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE

# Paths (repo-root relative, bound to /data inside the container)
R_SCRIPT="ANALYSIS/filter_contaminated_genes.R"
COUNTS="ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv"
LENGTHS="ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_lengths.tsv"
METADATA="ANALYSIS/metadata_base.csv"
OUT_DIR="ANALYSIS/results_decontamination"
CPM_CUTOFF="${1:-1}"        # arg 1: CPM cutoff for absolute mode (default 1)
MIN_CONTROLS="${2:-2}"      # arg 2: min controls for absolute mode (default 2)
MODE="${3:-absolute}"       # arg 3: "absolute" or "ratio" (default absolute)
RATIO_THRESHOLD="${4:-0.5}" # arg 4: control/tumor ratio for ratio mode (default 0.5)

mkdir -p "$OUT_DIR" slurm_logs

echo "Verifying input files..."
for path in "$R_SCRIPT" "$COUNTS" "$LENGTHS" "$METADATA"; do
    if [ ! -e "$path" ]; then echo "ERROR: not found: $path"; exit 1; fi
done
echo "All required files found."

IMG_PATH="${NXF_SINGULARITY_CACHEDIR}/go2432-bioconductor.sif"
if [[ ! -f "$IMG_PATH" ]]; then
    echo "Pulling Bioconductor container..."
    singularity pull "$IMG_PATH" docker://go2432/bioconductor:latest
fi
echo "Container ready: $IMG_PATH"
if [ "$MODE" = "ratio" ]; then
    echo "Rule: RATIO -- control/tumor CPM >= ${RATIO_THRESHOLD}"
else
    echo "Rule: ABSOLUTE -- CPM >= ${CPM_CUTOFF} in >= ${MIN_CONTROLS} controls"
fi

set +e
singularity exec --bind "$PWD:/data" --pwd /data "$IMG_PATH" \
    Rscript "/data/$R_SCRIPT" \
    "/data/$COUNTS" \
    "/data/$LENGTHS" \
    "/data/$METADATA" \
    "/data/$OUT_DIR" \
    "$CPM_CUTOFF" \
    "$MIN_CONTROLS" \
    "$MODE" \
    "$RATIO_THRESHOLD"
exit_code=$?
set -e

echo "================================================================"
echo " End      : $(date)   Exit code: $exit_code"
echo "================================================================"
if [ $exit_code -ne 0 ]; then
    echo "ERROR: decontamination filter failed (exit code: $exit_code)"
    exit $exit_code
fi

echo "Sweep (genes removed per cutoff):"
cat "$OUT_DIR/decontamination_sweep.csv" 2>/dev/null || true
echo ""
echo "Decontaminated matrices written next to the originals in star_salmon/."
echo "Next: sbatch run_de_pdx_v3.sh  (already points at the .decontaminated.tsv files)"
exit 0
