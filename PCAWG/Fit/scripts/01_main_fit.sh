#!/bin/bash
#SBATCH --job-name=PCAWG_inference_timing_all_cohort
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=70gb
#SBATCH --time=2:00:00
#SBATCH --array=500-1000
#SBATCH --output=log/output_%A_%a.out
#SBATCH --error=log/error_%A_%a.err

module load R/4.4.1

Rscript 01_main_fit.R ${SLURM_ARRAY_TASK_ID} 1170