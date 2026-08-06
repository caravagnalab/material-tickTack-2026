#SBATCH --job-name=basculePCAWG
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --partition=GENOA
#SBATCH --mem=32GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-30%30
 
# Conda — `process` provides Rscript and the base reticulate install;
# the R script then switches reticulate to the bascule-env conda env.
source ~/.bashrc
conda activate process
 
# Run from the directory containing the stage scripts so `source("00_setup.R")`
# inside 02_bascule_inference.R resolves correctly.
SCRIPT_DIR="/orfeo/cephfs/scratch/cdslab/gsantacatterina/material-tickTack-2026/PCAWG/Scripts/signature_deconvolution"
cd "${SCRIPT_DIR}"
 
mkdir -p logs
 
Rscript 02_bascule_inference.R \
  || { echo "Stage 02 failed for array index ${SLURM_ARRAY_TASK_ID}"; exit 1; }