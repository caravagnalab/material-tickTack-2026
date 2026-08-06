#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=Select_alpha
#SBATCH --ntasks=1                 
#SBATCH --mem=10G                  
#SBATCH --time=01:00:00     
#SBATCH --output=Select_alpha_%A_%a.out
#SBATCH --error=Select_alpha_%A_%a.err

base=/orfeo/scratch/cdslab/antonelloa/material-tickTack-2026/PCAWG/Scripts

module load R
Rscript $base/9_Select_alpha_for_signature_deconvolution.R --sample $SLURM_ARRAY_TASK_ID
