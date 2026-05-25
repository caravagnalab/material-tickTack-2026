#!/usr/bin/env bash
#SBATCH --job-name=sigtime_part
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --partition=GENOA
#SBATCH --mem=32GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-30%30

source ~/.bashrc
conda activate process

SCRIPT_DIR="/orfeo/cephfs/scratch/cdslab/gsantacatterina/material-tickTack-2026/PCAWG/Scripts/alpha_muts_signature_deconvolution"
cd "${SCRIPT_DIR}"
mkdir -p logs

Rscript 01_partition_and_matrix.R || { echo "Stage 01 (partition) failed for array index ${SLURM_ARRAY_TASK_ID}"; exit 1; }
