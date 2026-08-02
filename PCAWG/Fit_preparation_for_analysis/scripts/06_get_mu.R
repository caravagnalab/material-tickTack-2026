source("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit_preparation_for_analysis/scripts/06_clustering_confidence.R")
#!/usr/bin/env Rscript
# =============================================================================
# pipeline_timing_dtaumin.R
#
# End-to-end per-sample pipeline:
#   1. load a CNAqc object (mutations + CNA + purity)
#   2. mobster fit on the clonal VAF  ->  neutral-tail mutation rate mu
#   3. per timeable segment: mutation burden M, purity, coverage
#   4. tickTack timing (or reuse an existing results_timing)  ->  best_K
#   5. for each K=1 cluster: resolution limit Delta-tau_min  ->  classify the
#      single-cluster call as CERTIFIED-SYNCHRONOUS or CENSORED (underpowered)
#
# Usage:
#   Rscript pipeline_timing_dtaumin.R  <sample.rds>  [out_dir]
#
# Depends: CNAqc, mobster, tickTack, dplyr, tibble.  deltatau_min.R must sit
# next to this file (or edit DTM_PATH).  On the HPC these are already installed.
# =============================================================================

suppressMessages({
  library(dplyr); library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript pipeline_timing_dtaumin.R <sample.rds> [out_dir]")
sample_rds <- args[[1]]
out_dir    <- if (length(args) >= 2) args[[2]] else dirname(sample_rds)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

DTM_PATH <- file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
                      "deltatau_min.R")
if (!file.exists(DTM_PATH)) DTM_PATH <- "deltatau_min.R"
source(DTM_PATH)

# timeable karyotypes tickTack handles + their total copy number
TIMEABLE <- c("2:1" = 3, "2:2" = 4, "2:0" = 2)

# -----------------------------------------------------------------------------
# 1. LOAD
# -----------------------------------------------------------------------------
fit_list <- list.files("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/inference_results_5ncomponents")
sample_rds = paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/inference_results_5ncomponents/",fit_list[7])
message("[1] loading ", sample_rds)
x <- readRDS(sample_rds)
x$results_timing$data$accepted_cna
x$results_timing$best_K
stopifnot(!is.null(x$mutations), !is.null(x$purity))
purity <- as.numeric(x$purity)
muts   <- x$mutations
message("    purity = ", round(purity, 3), " | ", nrow(muts), " mutations")

# -----------------------------------------------------------------------------
# 2. MOBSTER  ->  mutation rate mu from the neutral (1/f) tail
#    mobster is fit on the clonal-diploid VAF spectrum; the Pareto tail gives
#    mu = mutations per effective genome doubling (evolutionary_parameters()).
# -----------------------------------------------------------------------------
mu <- NA_real_; tail_exponent <- NA_real_; mob <- NULL
mob_ok <- requireNamespace("mobster", quietly = TRUE)
if (mob_ok) {
  message("[2] mobster fit for mutation rate ...")
  # fit on 1:1 (diploid) clonal mutations -- the population mobster is designed for
  vaf_tbl <- muts %>%
    filter(!is.na(VAF), VAF > 0.05, VAF < 0.95) %>%
    { if ("karyotype" %in% names(.)) filter(., karyotype == "1:1") else . } %>%
    transmute(VAF = as.numeric(VAF))
  if (nrow(vaf_tbl) >= 30) {
    # NOTE: parallel = FALSE is required -- the parallel scoring backend throws
    # "missing value where TRUE/FALSE needed" on this mobster/conda build.
    mob <- tryCatch(
      mobster::mobster_fit(vaf_tbl, tail = c(TRUE, FALSE), K = 1:2,
                           samples = 3, init = "peaks", silent = TRUE,
                           parallel = FALSE),
      error = function(e) { message("    mobster_fit failed: ", conditionMessage(e)); NULL })
    best_fit = mob$best
    print(best_fit)
    plot(best_fit)
    
    if (!is.null(mob)) {
      ep <- tryCatch(mobster::evolutionary_parameters(mob, ploidy = 2),
                     error = function(e) { message("    no tail / evo params: ", conditionMessage(e)); NULL })
      if (!is.null(ep)) {
        mu <- ep$mu[1]; tail_exponent <- ep$exponent[1]
        message(sprintf("    mu = %.4g (muts / genome doubling) | tail exponent = %.2f", mu, tail_exponent))
      }
    }
  } else message("    too few diploid clonal mutations for mobster (", nrow(vaf_tbl), ")")
} else message("[2] mobster not installed -- skipping mu (Delta-tau_min still computed from M directly)")

# -----------------------------------------------------------------------------
# 3. TIMEABLE SEGMENTS  ->  per-segment M, coverage
# -----------------------------------------------------------------------------
message("[3] enumerating timeable segments")
if (!"karyotype" %in% names(muts) && !is.null(x$cna)) {
  # attach karyotype from CNA if not already on mutations
  message("    (mutations lack karyotype column; relying on x$cna)")
}


accepted_cna <- x$results_timing$data$accepted_cna
data_df <- data.frame(
  Segment = x$results_timing$data$input_data$seg_assignment,
  DP = x$results_timing$data$input_data$DP
)
summary_table <- data_df %>%
  group_by(Segment) %>%
  summarise(
    M = n(),
    coverage = round(mean(DP))
  ) %>% 
  rename(segment_id = Segment) %>% 
  mutate(segment_id = as.double(segment_id))

seg_tbl <- accepted_cna %>% dplyr::inner_join(summary_table, by="segment_id") %>% 
  rowwise() %>%
  mutate(CN_tot = sum( as.integer(strsplit(karyotype, split = ":")[[1]][1]), 
                       as.integer(strsplit(karyotype, split = ":")[[1]][2]))) %>%
  rename(segment_id_original = segment_id) %>% 
  rename(segment_id = segment_name)


# -----------------------------------------------------------------------------
# 4. TIMING (best_K).  Reuse x$results_timing if present, else run tickTack.
# -----------------------------------------------------------------------------
best_K <- NA_integer_
if (!is.null(x$results_timing) && !is.null(x$results_timing$best_K)) {
  best_K <- as.integer(x$results_timing$best_K)
  message("[4] reusing stored results_timing: best_K = ", best_K)
} else if (requireNamespace("tickTack", quietly = TRUE)) {
  message("[4] running tickTack timing ...")
  fit <- tryCatch({
    y <- tickTack::fit(x, max_attempts = 2, INIT = TRUE, tolerance = 0.001)
    y
  }, error = function(e) { message("    tickTack::fit failed: ", conditionMessage(e)); NULL })
  if (!is.null(fit)) {
    best_K <- tryCatch(as.integer(fit$results_timing$best_K), error = function(e) NA_integer_)
    x$results_timing <- fit$results_timing
    message("    best_K = ", best_K)
  }
} else message("[4] no stored timing and tickTack not installed -- best_K unknown")

# -----------------------------------------------------------------------------
# 5. RESOLUTION LIMIT + CLASSIFICATION
#    For a single-clock (K=1) call across S_timeable segments, the shared clock
#    is supported by all timeable segments (pooling), while each segment's own
#    tau precision is set by its M. We report:
#      - per-segment Delta-tau_min (S=1, the segment alone)
#      - cluster Delta-tau_min using pooled S = S_timeable and the median M
#    and classify against an equivalence margin.
# -----------------------------------------------------------------------------
message("[5] resolution limit")

EQ_MARGIN <- 0.10   # biological equivalence margin ("together" = gap < 0.10). ADJUST in Methods.
tau_eval  <- 0.5    # weakest-information point; conservative

seg_res <- seg_tbl %>%
  rowwise() %>%
  mutate(
    dtau_min80_seg = dtau_min(M, purity, coverage, tau = tau_eval, S = 1, target = "80"),
    dtau_min50_seg = dtau_min(M, purity, coverage, tau = tau_eval, S = 1, target = "50")
  ) %>% ungroup()

S_timeable = seg_tbl %>% nrow()
M_med   <- stats::median(seg_tbl$M)
cov_med <- stats::median(seg_tbl$coverage)
dtau_cluster80 <- dtau_min(M_med, purity, round(cov_med), tau = tau_eval, S = S_timeable, target = "80")
dtau_cluster50 <- dtau_min(M_med, purity, round(cov_med), tau = tau_eval, S = S_timeable, target = "50")
# cap at 1: a value >=1 means the whole [0,1] timeline is unresolvable
dtau_cluster80c <- min(dtau_cluster80, 1)

classification <- if (is.na(best_K)) {
  "UNKNOWN (no timing result)"
} else if (best_K >= 2) {
  "RESOLVED (>=2 clocks) -- asynchrony trustworthy (splitting essentially never spurious)"
} else if (dtau_cluster80 <= EQ_MARGIN) {
  sprintf("CERTIFIED-SYNCHRONOUS: had power to resolve gaps > %.2f, saw one clock -> synchronous to within %.2f", dtau_cluster80, dtau_cluster80)
} else {
  sprintf("CENSORED / UNDERPOWERED: could only resolve gaps > %.2f (> margin %.2f) -> consistent with synchrony but NOT certified", dtau_cluster80c, EQ_MARGIN)
}

# -----------------------------------------------------------------------------
# OUTPUT
# -----------------------------------------------------------------------------
sample_id <- if ("sample" %in% names(muts)) muts$sample[1] else basename(sample_rds)
summary_tbl <- tibble(
  sample = sample_id, purity = purity, mu_neutral = mu, tail_exponent = tail_exponent,
  n_timeable_segments = S_timeable, M_median = M_med, coverage_median = cov_med,
  best_K = best_K, eq_margin = EQ_MARGIN,
  dtau_min80_cluster = dtau_cluster80c, dtau_min50_cluster = min(dtau_cluster50, 1),
  classification = classification
)

message("\n================= SAMPLE SUMMARY =================")
print(as.data.frame(summary_tbl))
message("\n--- per-segment resolution ---")
print(as.data.frame(seg_res))
message("\nCLASSIFICATION: ", classification)

out_csv <- file.path(out_dir, paste0(sample_id, "__dtaumin_summary.csv"))
seg_csv <- file.path(out_dir, paste0(sample_id, "__dtaumin_persegment.csv"))
write.csv(summary_tbl, out_csv, row.names = FALSE)
write.csv(seg_res,     seg_csv, row.names = FALSE)
message("\nwrote: ", out_csv, "\n       ", seg_csv)
