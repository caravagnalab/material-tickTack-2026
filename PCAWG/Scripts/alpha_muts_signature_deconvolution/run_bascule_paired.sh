#!/usr/bin/env bash
# ==============================================================================
# run_bascule_paired.sh — submit stage 02 of the time-stratified analysis
# Chain after stage 01:
#   sbatch --dependency=afterok:<stage1_jobid> run_bascule_paired.sh
# ==============================================================================
#SBATCH --job-name=sigtime_basc
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --partition=GENOA
#SBATCH --mem=32GB
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-30%30

source ~/.bashrc
conda activate process

SCRIPT_DIR="/orfeo/cephfs/scratch/cdslab/gsantacatterina/material-tickTack-2026/PCAWG/Scripts/alpha_muts_signature_deconvolution"
cd "${SCRIPT_DIR}"
mkdir -p logs

Rscript 02_bascule_paired.R || { echo "Stage 02 (bascule) failed for array index ${SLURM_ARRAY_TASK_ID}"; exit 1; }
