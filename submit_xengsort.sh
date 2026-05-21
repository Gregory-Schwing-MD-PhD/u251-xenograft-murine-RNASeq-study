#!/bin/bash
#SBATCH -q secondary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=go2432@wayne.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o xengsort_prod_%j.out
#SBATCH -e xengsort_prod_%j.err

# 1. Environment Setup
source "${HOME}/mambaforge/etc/profile.d/conda.sh"
source activate nextflow
export PATH="${HOME}/mambaforge/envs/nextflow/bin:$PATH"
unset JAVA_HOME

# 2. Path Setup
# Fixes Singularity "Permission Denied" in /run/user/
export XDG_RUNTIME_DIR="${HOME}/xdr"
mkdir -p $XDG_RUNTIME_DIR

# Nextflow work directory (scratch for intermediate files)
export WORK_DIR="$(pwd)/work"
mkdir -p $WORK_DIR

# 3. Singularity Image Prep
# We pre-pull the image into your cache to avoid the "Stream Closed"
# or timeout errors that happen during the Nextflow launch phase.
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
mkdir -p $NXF_SINGULARITY_CACHEDIR
echo "Pre-pulling Singularity image from Docker Hub..."
singularity pull --name ${NXF_SINGULARITY_CACHEDIR}/go2432-xengsort-latest.img docker://go2432/xengsort:latest

# 4. Clean up stale Nextflow locks
find .nextflow/cache -name "LOCK" -delete 2>/dev/null

# 5. Production Run
# Using -profile singularity to trigger the container block in nextflow.config
echo "RUNNING STEP 1: XENGSORT"
nextflow run main.nf -profile singularity \
    --input "ANALYSIS/samplesheet.csv" \
    --host_fasta "ANALYSIS/refs/rat/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz" \
    --graft_fasta "ANALYSIS/refs/human/GRCh38.primary_assembly.genome.fa.gz" \
    --outdir "ANALYSIS" \
    -w "${WORK_DIR}" \
    -resume -ansi-log false
