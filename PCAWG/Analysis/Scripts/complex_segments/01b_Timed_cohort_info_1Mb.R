.libPaths("~/R/orfeo_R_4.4/")
# -----------------------------------------------------------------
Drivers <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Drivers.rds")
Samples <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Samples.rds")
Segments <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Segments.rds")

# ===================================================================
# Karyotype landscape analysis across samples and cancer types
# ===================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(stringr)
library(forcats)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(grid)
library(gridExtra)
library(scales)
library(png)
library(data.table)

karyotypes <- c("1:0", "1:1", "2:0", "2:1", "2:2", "3:0", "3:1", "3:2", "3:3",
                "4:0", "4:1", "4:2", "4:3", "4:4", "5:0", "5:1", "5:2", "5:3", "5:4",
                "6:0", "6:1", "6:2", "6:3", "7:0", "7:1", "7:2", "8:0", "8:1", "9:0")

karyotype_levels <- karyotypes

chr_order <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

get_mode <- function(x) {
  if (length(x) == 0) return(NA)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# --- Helper: render a pheatmap object to a raster grob via a temp PNG ---
# This sidesteps pheatmap's internal fixed viewport names (e.g. "matrix"),
# which otherwise collide across multiple captures on the same page/device
# and cause every subsequent heatmap to silently render as the first one.
pheatmap_to_raster <- function(pheatmap_obj, width_px = 1600, height_px = 1200, res = 150) {
  tmp_png <- tempfile(fileext = ".png")
  png(tmp_png, width = width_px, height = height_px, res = res)
  grid.draw(pheatmap_obj$gtable)
  dev.off()
  
  img <- readPNG(tmp_png)
  file.remove(tmp_png)
  
  rasterGrob(img, interpolate = TRUE)
}

summarise_karyotypes <- function(df, karyotypes) {
  
  df <- df %>%
    mutate(karyotype = paste(Major, minor, sep = ":"))
  
  df_filtered <- df %>%
    filter(karyotype %in% karyotypes)
  
  # Helper to safely compute min/max, returning NA instead of Inf/-Inf
  # when a group has no non-missing values
  safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
  safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  
  karyotype_summary <- df_filtered %>%
    group_by(karyotype) %>%
    summarise(
      n_segments    = n(),
      
      length_min    = safe_min(length),
      length_max    = safe_max(length),
      length_median = median(length, na.rm = TRUE),
      
      n_min    = safe_min(n),
      n_max    = safe_max(n),
      n_median = median(n, na.rm = TRUE),
      
      chr_mode = get_mode(chr),
      
      .groups = "drop"
    ) %>%
    right_join(tibble(karyotype = karyotypes), by = "karyotype") %>%
    mutate(n_segments = ifelse(is.na(n_segments), 0, n_segments)) %>%
    arrange(factor(karyotype, levels = karyotypes))
  
  chr_distribution <- df_filtered %>%
    count(karyotype, chr, name = "n_segments_chr") %>%
    arrange(factor(karyotype, levels = karyotypes), desc(n_segments_chr))
  
  list(summary = karyotype_summary, chr_distribution = chr_distribution)
}

# --- Helper: count mutations falling within each CNA segment ---
add_mutation_counts <- function(cna, mutations) {
  
  if (nrow(mutations) == 0 || nrow(cna) == 0) {
    cna$n <- 0L
    return(as_tibble(cna))
  }
  
  cna_dt <- as.data.table(cna)
  mut_dt <- as.data.table(mutations)
  
  cna_dt[, seg_id := .I]
  
  # Ensure numeric types line up (readRDS sometimes leaves these as dbl, which is fine for foverlaps)
  mut_dt[, `:=`(mut_start = from, mut_end = to)]
  
  setkey(cna_dt, chr, from, to)
  
  overlaps <- foverlaps(
    mut_dt, cna_dt,
    by.x = c("chr", "mut_start", "mut_end"),
    by.y = c("chr", "from", "to"),
    type = "within",
    nomatch = 0L
  )
  
  n_per_segment <- overlaps[, .N, by = seg_id]
  
  cna_dt <- merge(cna_dt, n_per_segment, by = "seg_id", all.x = TRUE)
  cna_dt[is.na(N), N := 0L]
  setnames(cna_dt, "N", "n")
  cna_dt[, seg_id := NULL]
  
  as_tibble(cna_dt)
}

run_multi_sample_summary <- function(sample_list, karyotypes) {
  
  results <- imap(sample_list, function(df, sample_name) {
    res <- summarise_karyotypes(df, karyotypes)
    res$summary <- res$summary %>% mutate(sample = sample_name, .before = 1)
    res$chr_distribution <- res$chr_distribution %>% mutate(sample = sample_name, .before = 1)
    res
  })
  
  list(
    summary_all = bind_rows(map(results, "summary")),
    chr_distribution_all = bind_rows(map(results, "chr_distribution"))
  )
}

base_dir <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/inference_results_5ncomponents/"
sample_list <- list()
  
for (idx in seq_len(nrow(Samples))) {
  # for (idx in 1:10) {
    
    sample_id <- Samples$sample[idx]
    
    fit_path <- paste0(base_dir, "/",sample_id, ".rds")
    
    if (!file.exists(fit_path)) {
      message("Missing: ", fit_path)
      next
    }
    
    fit <- readRDS(fit_path)
    
    cna_df <- fit$cna %>% filter(length >= 1000000)
    
    if (!"length" %in% colnames(cna_df)) {
      cna_df <- cna_df %>% mutate(length = to - from)
    }
    if (!"n" %in% colnames(cna_df)) {
      cna_df <- add_mutation_counts(cna_df, fit$mutations)
    }
    
    df <- cna_df
    df <- df %>% mutate(sample = sample_id, IntoGen_cancer_type = Samples$IntoGen_cancer_type[idx])
    
    sample_list[[sample_id]] <- df
  }

saveRDS(sample_list, "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Analysis_results/Cohort_stat/complex_segments/sample_list_1Mb.RDS")
multi_results <- run_multi_sample_summary(sample_list, karyotypes)
saveRDS(multi_results, "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Analysis_results/Cohort_stat/complex_segments/multi_results_1Mb.RDS")
multi_results <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Analysis_results/Cohort_stat/complex_segments/multi_results_1Mb.RDS")


summary_all <- multi_results$summary_all
chr_distribution_all <- multi_results$chr_distribution_all

# Raw (segment-level) table with karyotype + sample + cancer_type tagged,
# needed for the boxplots and chromosome composition plots
raw_all <- imap_dfr(sample_list, function(df, sample_name) {
  df %>%
    mutate(karyotype = paste(Major, minor, sep = ":")) %>%
    filter(karyotype %in% karyotypes)
})

# Bring cancer_type into summary_all if not already present
if (!"IntoGen_cancer_type" %in% colnames(summary_all)) {
  sample_to_cancer <- raw_all %>% distinct(sample, IntoGen_cancer_type)
  summary_all <- summary_all %>% left_join(sample_to_cancer, by = "sample")
}

count_matrix <- summary_all %>%
  select(sample, karyotype, n_segments) %>%
  mutate(karyotype = factor(karyotype, levels = karyotype_levels)) %>%
  pivot_wider(names_from = karyotype, values_from = n_segments, values_fill = 0) %>%
  column_to_rownames("sample") %>%
  as.matrix()

count_matrix <- count_matrix[, karyotype_levels[karyotype_levels %in% colnames(count_matrix)]]

prop_matrix <- count_matrix / rowSums(count_matrix)
prop_matrix[is.nan(prop_matrix)] <- 0

pheatmap(
  prop_matrix,
  color = colorRampPalette(c("white", "#3B4CC0", "#B40426"))(100),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "ward.D2",
  cluster_rows              = TRUE,
  cluster_cols              = FALSE,
  fontsize_row               = 8,
  fontsize_col               = 9,
  angle_col                  = 90,
  main                        = "Karyotype proportion profile per sample (clustered)",
  border_color                = NA
)

# -----------------------------------------------------------------
# Combined single figure (all samples pooled)
# -----------------------------------------------------------------

pdf("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Analysis_results/Cohort_stat/complex_segments/Timed_cohort_karyotype_all_cancer_types_1MB.pdf", width = 16, height = 12)
# pdf("~/old_home/re_gecip/cancer_pan/scocomello/material-tickTack-2026/Results/Cohort_stat/complex_segments/test_1.pdf", width = 16, height = 12)

chr_levels_present <- chr_order[chr_order %in% unique(raw_all$chr)]

p_chr <- raw_all %>%
  mutate(
    karyotype = factor(karyotype, levels = karyotype_levels),
    chr       = factor(chr, levels = chr_levels_present)
  ) %>%
  count(karyotype, chr) %>%
  group_by(karyotype) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = karyotype, y = prop, fill = chr)) +
  geom_col() +
  labs(title = "Chromosome composition per karyotype (proportion, all samples)",
       x = "Karyotype", y = "Proportion of segments", fill = "Chr") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 3))

p_length <- raw_all %>%
  mutate(karyotype = factor(karyotype, levels = karyotype_levels)) %>%
  ggplot(aes(x = karyotype, y = length)) +
  geom_boxplot(outlier.size = 0.5, fill = "steelblue", alpha = 0.6) +
  scale_y_log10(labels = scales::comma) +
  labs(title = "Segment length distribution per karyotype (all samples)",
       x = "Karyotype", y = "Length (log10 scale)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p_n <- raw_all %>%
  mutate(karyotype = factor(karyotype, levels = karyotype_levels)) %>%
  ggplot(aes(x = karyotype, y = n)) +
  geom_boxplot(outlier.size = 0.5, fill = "darkorange", alpha = 0.6) +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = round(after_stat(y), 1)),
    vjust = -0.6,
    size = 3,
    color = "black"
  ) +
  scale_y_log10(labels = scales::comma) +
  labs(title = "n distribution per karyotype (all samples)",
       x = "Karyotype", y = "n (log10 scale)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# --- FIX: render pheatmap to a raster grob (via temp PNG) instead of
#     capturing the gtable directly, avoiding pheatmap's internal fixed
#     viewport names colliding across multiple captures ---
heatmap_pheatmap_obj <- pheatmap(
  prop_matrix,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "ward.D2",
  cluster_rows              = TRUE,
  cluster_cols              = FALSE,
  fontsize_row               = 7,
  fontsize_col               = 8,
  angle_col                  = 90,
  main                        = "Karyotype profile (clustered samples)",
  border_color                = NA,
  silent                       = TRUE
)

heatmap_grob <- pheatmap_to_raster(heatmap_pheatmap_obj)
p_heatmap_wrapped <- wrap_elements(full = heatmap_grob)

combined_fig <- (p_heatmap_wrapped | p_chr) /
  (p_length | p_n)

combined_fig <- combined_fig +
  plot_annotation(title = "Karyotype landscape across all samples",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))

print(combined_fig)

dev.off()

# -----------------------------------------------------------------
# Per-cancer-type figures, all saved as pages in one PDF
# -----------------------------------------------------------------

cancer_types_present <- unique(raw_all$IntoGen_cancer_type)

pdf(paste0("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Analysis_results/Cohort_stat/complex_segments/", "karyotype_by_cancer_type_1Mb.pdf"), width = 16, height = 12)

for (ct in cancer_types_present) {
  
  raw_ct     <- raw_all     %>% filter(IntoGen_cancer_type == ct)
  summary_ct <- summary_all %>% filter(IntoGen_cancer_type == ct)
  
  if (nrow(raw_ct) == 0) next
  
  count_matrix_ct <- summary_ct %>%
    select(sample, karyotype, n_segments) %>%
    mutate(karyotype = factor(karyotype, levels = karyotype_levels)) %>%
    pivot_wider(names_from = karyotype, values_from = n_segments, values_fill = 0) %>%
    column_to_rownames("sample") %>%
    as.matrix()
  
  count_matrix_ct <- count_matrix_ct[, karyotype_levels[karyotype_levels %in% colnames(count_matrix_ct)], drop = FALSE]
  
  prop_matrix_ct <- count_matrix_ct / rowSums(count_matrix_ct)
  prop_matrix_ct[is.nan(prop_matrix_ct)] <- 0
  
  cluster_ok <- nrow(prop_matrix_ct) >= 2
  
  # --- FIX: same raster-based capture here ---
  pheatmap_obj_ct <- pheatmap(
    prop_matrix_ct,
    cluster_rows              = cluster_ok,
    cluster_cols               = FALSE,
    clustering_distance_rows  = "euclidean",
    clustering_method         = "ward.D2",
    fontsize_row               = 7,
    fontsize_col               = 8,
    angle_col                  = 90,
    main                        = paste0(ct, ": karyotype profile"),
    border_color                = NA,
    silent                       = TRUE
  )
  
  heatmap_grob_ct <- pheatmap_to_raster(pheatmap_obj_ct)
  p_heatmap_ct_wrapped <- wrap_elements(full = heatmap_grob_ct)
  
  chr_levels_present_ct <- chr_order[chr_order %in% unique(raw_ct$chr)]
  
  p_chr_ct <- raw_ct %>%
    mutate(
      karyotype = factor(karyotype, levels = karyotype_levels),
      chr       = factor(chr, levels = chr_levels_present_ct)
    ) %>%
    count(karyotype, chr) %>%
    group_by(karyotype) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = karyotype, y = prop, fill = chr)) +
    geom_col() +
    labs(title = paste0(ct, ": chromosome composition"),
         x = "Karyotype", y = "Proportion", fill = "Chr") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 3))
  
  p_length_ct <- raw_ct %>%
    mutate(karyotype = factor(karyotype, levels = karyotype_levels)) %>%
    ggplot(aes(x = karyotype, y = length)) +
    geom_boxplot(outlier.size = 0.5, fill = "steelblue", alpha = 0.6) +
    scale_y_log10(labels = scales::comma) +
    labs(title = paste0(ct, ": segment length"),
         x = "Karyotype", y = "Length (log10)") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  p_n_ct <- raw_ct %>%
    mutate(karyotype = factor(karyotype, levels = karyotype_levels)) %>%
    ggplot(aes(x = karyotype, y = n)) +
    scale_y_log10(labels = scales::comma) +
    geom_boxplot(outlier.size = 0.5, fill = "darkorange", alpha = 0.6) +
    stat_summary(
      fun = median,
      geom = "text",
      aes(label = round(after_stat(y), 1)),
      vjust = -0.6,
      size = 3,
      color = "black"
    ) +
    labs(title = paste0(ct, ": n distribution"),
         x = "Karyotype", y = "n (log10 scale)") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  combined_fig_ct <- (p_heatmap_ct_wrapped | p_chr_ct) /
    (p_length_ct | p_n_ct) +
    plot_annotation(
      title = paste0("Karyotype landscape: ", ct),
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  print(combined_fig_ct)
}

dev.off()
