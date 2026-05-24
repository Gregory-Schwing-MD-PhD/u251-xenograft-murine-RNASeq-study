#!/bin/bash
#SBATCH -q secondary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --mail-user=go2432@wayne.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o de_prod_%j.out
#SBATCH -e de_prod_%j.err

set -x

# 1. Environment Setup
source "${HOME}/mambaforge/etc/profile.d/conda.sh"
source activate nextflow
export PATH="${HOME}/mambaforge/envs/nextflow/bin:$PATH"
unset JAVA_HOME

# 2. Path Setup
export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
mkdir -p $XDG_RUNTIME_DIR $NXF_SINGULARITY_CACHEDIR

WORK_DIR="$(pwd)/work"

# 3. Safety locks
export NXF_SINGULARITY_HOME_MOUNT=true
unset LD_LIBRARY_PATH
unset PYTHONPATH
unset R_LIBS
unset R_LIBS_USER
unset R_LIBS_SITE

# 4. Production Command
nextflow run nf-core/differentialabundance \
    -r 1.5.0 \
    -profile singularity \
    --input "$(pwd)/ANALYSIS/metadata_therapy.csv" \
    --contrasts "$(pwd)/ANALYSIS/contrasts_therapy.csv" \
    --matrix "$(pwd)/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.decontaminated.tsv" \
    --transcript_length_matrix "$(pwd)/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_lengths.decontaminated.tsv" \
    --shinyngs_build_app \
    -params-file therapy_v3_params.yaml \
    -w "${WORK_DIR}" \
    -resume \
    -ansi-log false
