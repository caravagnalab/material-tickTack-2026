# ==============================================================================
# 02_bascule_inference.R
# ------------------------------------------------------------------------------
# Pipeline stage 2: BASCULE signature inference + clustering.
#
# Depends on stage 01 having already produced the SBS96 count matrix at:
#   <out_data_tt>/output/SBS/PCAWG.SBS96.all
#
# This script:
#   1. loads the SBS96 matrix and transposes it to samples-as-rows
#   2. fits BASCULE against the COSMIC SBS catalogue, anchoring SBS1/SBS5
#   3. refines de novo signatures and clusters samples (k = 3, then merge)
#   4. plots exposures, joins clusters with sample metadata, and saves
#      RDS + PDF outputs under <out_data_tt>/bascule/
#
# Submitted as a job array, one task per tumour type:
#   sbatch --array=1-30 run_bascule.sh
# ==============================================================================

# Shared setup: defines tt, samples_info, project, out_data_tt.
source("00_setup.R")

# ── Python / BASCULE environment ──────────────────────────────────────────────
# BASCULE has a Pyro/Python backend; it requires a different conda env than
# stage 01. The pybascule import is done once here to fail fast if the env
# is misconfigured (rather than discovering it deep inside bascule::fit).
reticulate::use_condaenv("/orfeo/scratch/cdslab/ggandolfi/miniconda/envs/bascule-env")
py <- reticulate::import("pybascule")

suppressPackageStartupMessages({
  library(bascule)
})

# ── Load the SBS96 matrix produced by stage 01 ────────────────────────────────
data_sbs <- file.path(out_data_tt, "output", "SBS", paste0(project, ".SBS96.all"))

if (!file.exists(data_sbs)) {
  stop("SBS96 matrix not found at ", data_sbs,
       ". Did stage 01_matrix_generation.R run successfully for ", tt, "?")
}

sig_matrix <- read.table(data_sbs, header = TRUE)
rownames(sig_matrix) <- sig_matrix$MutationType

# BASCULE expects samples-as-rows / contexts-as-columns; SigProfiler outputs
# the opposite, so we drop the MutationType column and transpose.
sig_counts <- sig_matrix %>%
  dplyr::select(-MutationType) %>%
  t() %>%
  as.data.frame()

# Drop samples that ended up with zero counts — they break the fit.
sig_counts <- sig_counts[rowSums(sig_counts) > 0, ]

# ── BASCULE fit ───────────────────────────────────────────────────────────────
# We fit against the COSMIC SBS catalogue and force SBS1 and SBS5 to be
# retained (these are near-ubiquitous "clock" signatures; pinning them
# stabilises the de novo discovery for the remaining residual).
cat   <- list(SBS = bascule::COSMIC_sbs %>% as.data.frame())
input <- list(SBS = sig_counts)

x <- bascule::fit(
  counts           = input,
  k_list           = 0,                        # no extra de novo at this stage
  reference_cat    = cat,
  keep_sigs        = c("SBS1", "SBS5"),        # anchored clock signatures
  hyperparameters  = NULL,
  lr               = 0.005,
  optim_gamma      = 0.1,
  n_steps          = 3000,
  py               = NULL,
  enumer           = "parallel",
  nonparametric    = TRUE,
  autoguide        = FALSE,
  filter_dn        = FALSE,
  min_exposure     = 0.1,
  CUDA             = FALSE,
  compile          = FALSE,
  store_parameters = FALSE,
  store_fits       = TRUE,
  seed_list        = 10
)

# ── Refine + cluster ──────────────────────────────────────────────────────────
# refine_denovo_signatures() picks up any unexplained residual as new de novo
# signatures; fit_clustering()/merge_clusters() then groups samples by their
# exposure profiles. k = 3 is the starting point; merge_clusters() collapses
# clusters that are too similar to be statistically distinct.
x_refined         <- refine_denovo_signatures(x)
x_refined_cluster <- fit_clustering(x_refined, cluster = 3)
x_refined_cluster <- merge_clusters(x_refined_cluster)

# ── Diagnostic plots ──────────────────────────────────────────────────────────
plot_exp         <- plot_exposures(x = x,                 sample_name = TRUE)
plot_exp_refined <- plot_exposures(x = x_refined_cluster, sample_name = TRUE)

# ── Join cluster assignments back to sample metadata ──────────────────────────
# Used by 03_downstream_analysis.R for enrichment tests and oncoprints.
sample_info_with_clusters <- x_refined_cluster$clustering$clusters %>%
  dplyr::rename(sampleId = samples) %>%
  dplyr::left_join(samples_info, by = "sampleId")

# ── Persist everything ────────────────────────────────────────────────────────
final_dir <- file.path(out_data_tt, "bascule")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(list(x = x, x_refined_cluster = x_refined_cluster),
        file.path(final_dir, "fits.rds"))

saveRDS(list(plot_exp = plot_exp, plot_exp_refined = plot_exp_refined),
        file.path(final_dir, "plots.rds"))

saveRDS(sample_info_with_clusters,
        file.path(final_dir, "sample_clusters.rds"))

pdf(file.path(final_dir, "exposures.pdf"), width = 12, height = 5)
print(plot_exp)
print(plot_exp_refined)
dev.off()

message("Stage 02 done for tumour type: ", tt)
