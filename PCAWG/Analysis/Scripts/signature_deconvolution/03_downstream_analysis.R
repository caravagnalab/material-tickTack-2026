# ==============================================================================
# 03_downstream_analysis.R
# ------------------------------------------------------------------------------
# Pipeline stage 3: downstream interpretation of BASCULE clusters.
#
# Unlike stages 01 and 02, this script is NOT a job-array task. It iterates
# over every tumour type for which stage 02 has produced cluster assignments,
# and for each one runs:
#
#   * a chi-squared homogeneity test of `class` vs cluster, with post-hoc
#     standardized residuals (FDR-corrected) to find enriched/depleted cells
#   * extraction of BRCA1 SNVs from CNAqc objects (hard-coded coordinates,
#     GRCh37)
#   * parsing of per-sample VEP annotations, with Ensembl-ID → symbol
#     mapping, filtered to BRCA1
#   * a ComplexHeatmap oncoprint of driver mutations, split by class & cluster
#
# Intended to be run interactively or as a single non-array job. It does not
# use the SLURM/SGE array index.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
  library(ComplexHeatmap)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# ── Paths (kept in sync with 00_setup.R) ──────────────────────────────────────
# This script does NOT source 00_setup.R because that script reads a job-array
# index from the environment and stops if it isn't set. We redefine the small
# set of paths we need here instead.
samples_rds <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/material-tickTack-2026/PCAWG/Data/Samples.rds"
out_data    <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/material-tickTack-2026/PCAWG/Analysis_results/Signatures/"
fit_dir     <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/inference_results_5ncomponents/"
drivers_rds <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/material-tickTack-2026/PCAWG/Data/Drivers.rds"
vep_root    <- "/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/selection/clonal_analysis_PCAWG/"

samples_info <- readRDS(samples_rds)
ttypes       <- unique(samples_info$ttype)


# ==============================================================================
# Helper functions
# ==============================================================================

# ── test_cluster_enrichment ────────────────────────────────────────────────────
# Tests whether `class` is distributed homogeneously across `cluster`, and if
# not, reports which (cluster, class) cells are over- or under-represented.
#
# Method:
#   1. Chi-squared test of independence on the contingency table.
#   2. If global p < alpha, compute standardized Pearson residuals (stdres).
#      |stdres| > 2 is nominally significant; > 3.29 ~ Bonferroni-corrected.
#   3. Convert each residual to a two-sided p-value via the standard normal,
#      then FDR-correct (Benjamini–Hochberg) across all cells.
#
# Returns the chi² result and the per-cell residual table (invisibly).
test_cluster_enrichment <- function(df,
                                    cluster_col = "cluster",
                                    class_col   = "class",
                                    alpha       = 0.05) {

  contingency <- table(df[[cluster_col]], df[[class_col]])
  chi_result  <- chisq.test(contingency)

  # Step 1: global homogeneity.
  cat(sprintf("Chi-squared test of homogeneity: X^2=%.3f, dof=%d, p=%.4e\n",
              chi_result$statistic, chi_result$parameter, chi_result$p.value))
  if (chi_result$p.value >= alpha) {
    cat("  -> No evidence that class incidence varies across clusters.\n")
    return(invisible(list(chi2_test = chi_result)))
  }
  cat("  -> Class incidence is NOT constant across clusters. Running post-hoc analysis...\n\n")

  # Step 2: standardized residuals -> per-cell tests.
  residuals <- chi_result$stdres

  results <- as.data.frame(as.table(residuals)) %>%
    setNames(c("cluster", "class", "std_residual")) %>%
    dplyr::mutate(
      observed        = as.vector(contingency),
      expected        = as.vector(chi_result$expected),
      fold_enrichment = observed / expected,
      enriched        = std_residual > 0,
      # Two-sided p from std residual under N(0,1).
      p_value         = 2 * pnorm(abs(std_residual), lower.tail = FALSE),
      p_adj           = p.adjust(p_value, method = "BH"),
      significant     = p_adj < alpha
    ) %>%
    dplyr::arrange(desc(std_residual))

  cat("Significant enrichments/depletions (FDR corrected):\n")
  results %>%
    dplyr::filter(significant) %>%
    dplyr::mutate(direction = ifelse(enriched, "enriched", "depleted")) %>%
    dplyr::select(cluster, class, direction, observed, expected,
                  fold_enrichment, std_residual, p_adj) %>%
    print()

  invisible(list(chi2_test = chi_result, residuals = results))
}


# ── plot_driver_oncoprint ─────────────────────────────────────────────────────
# Builds a ComplexHeatmap oncoprint of driver alterations for one tumour type,
# with columns split by (class | cluster) so visual blocks line up with the
# enrichment test above.
#
# Alteration encoding (from Drivers$mutation_status):
#   M          -> point mutation             (red, full bar)
#   CI_M       -> clonal/subclonal mutation  (salmon, half-height bar)
#   CNA_driver -> copy-number driver         (blue, thin horizontal line)
#
# Multiple alterations on the same gene/sample are joined by ";", which is
# the format oncoPrint() expects.
plot_driver_oncoprint <- function(clusters, Drivers,
                                  sample_col  = "sample",
                                  cluster_col = "clusters",
                                  class_col   = "class") {

  # 1. Restrict drivers to non-WT calls and attach cluster/class info.
  drivers_annotated <- Drivers %>%
    dplyr::filter(mutation_status != "WT") %>%
    inner_join(clusters %>% dplyr::select(all_of(c(sample_col, cluster_col, class_col))),
               by = c("sample_id" = sample_col))

  samples <- unique(clusters[[sample_col]])
  genes   <- unique(drivers_annotated$gene)

  # 2. Build gene x sample alteration matrix. Cells are ""/single alt/multi-alt
  #    strings joined by ";", as required by oncoPrint().
  mat <- matrix("", nrow = length(genes), ncol = length(samples),
                dimnames = list(genes, samples))

  for (i in seq_len(nrow(drivers_annotated))) {
    row <- drivers_annotated[i, ]
    g   <- row$gene
    s   <- row$sample_id
    alt <- row$mutation_status
    if (mat[g, s] == "") {
      mat[g, s] <- alt
    } else {
      mat[g, s] <- paste(mat[g, s], alt, sep = ";")
    }
  }

  # 3. Order samples by class then cluster so the column split is contiguous.
  sample_meta <- clusters %>%
    dplyr::select(all_of(c(sample_col, cluster_col, class_col))) %>%
    dplyr::distinct() %>%
    dplyr::arrange(.data[[class_col]], .data[[cluster_col]])

  sample_order <- sample_meta[[sample_col]]
  sample_order <- sample_order[sample_order %in% colnames(mat)]
  mat <- mat[, sample_order, drop = FALSE]

  # Keep genes with at least one hit, then sort by mutation frequency.
  mat <- mat[rowSums(mat != "") > 0, , drop = FALSE]
  gene_order <- order(rowSums(mat != ""), decreasing = TRUE)
  mat <- mat[gene_order, , drop = FALSE]

  # 4. Alteration colours and the per-alteration draw functions.
  col <- c(
    "M"          = "#E64B35",
    "CI_M"       = "#F39B7F",
    "CNA_driver" = "#4DBBD5"
  )

  alter_fun <- list(
    background = function(x, y, w, h)
      grid.rect(x, y, w * 0.9, h * 0.9,
                gp = gpar(fill = "#F5F5F5", col = NA)),
    M = function(x, y, w, h)
      grid.rect(x, y, w * 0.9, h * 0.9,
                gp = gpar(fill = col["M"], col = NA)),
    CI_M = function(x, y, w, h)
      grid.rect(x, y, w * 0.9, h * 0.45,
                gp = gpar(fill = col["CI_M"], col = NA)),
    CNA_driver = function(x, y, w, h)
      grid.segments(x - w * 0.4, y, x + w * 0.4, y,
                    gp = gpar(col = col["CNA_driver"], lwd = 3))
  )

  # 5. Top annotation bars: class and cluster, with distinct palettes.
  ann_df <- sample_meta %>%
    filter(.data[[sample_col]] %in% colnames(mat)) %>%
    arrange(match(.data[[sample_col]], colnames(mat)))

  n_clusters <- length(unique(ann_df[[cluster_col]]))
  n_classes  <- length(unique(ann_df[[class_col]]))

  cluster_cols <- setNames(
    colorspace::qualitative_hcl(n_clusters, palette = "Dark 3"),
    unique(ann_df[[cluster_col]])
  )
  class_cols <- setNames(
    colorspace::qualitative_hcl(n_classes, palette = "Set 2"),
    unique(ann_df[[class_col]])
  )

  top_ann <- HeatmapAnnotation(
    class   = ann_df[[class_col]],
    cluster = ann_df[[cluster_col]],
    col     = list(class   = class_cols,
                   cluster = cluster_cols),
    annotation_name_side = "left",
    gap = unit(1, "mm")
  )

  # 6. Column splits: one block per (class, cluster) combination.
  col_split <- factor(
    paste(ann_df[[class_col]], ann_df[[cluster_col]], sep = " | "),
    levels = unique(paste(ann_df[[class_col]], ann_df[[cluster_col]], sep = " | "))
  )

  # 7. Draw and return invisibly so the caller can post-process if needed.
  op <- oncoPrint(
    mat,
    alter_fun       = alter_fun,
    col             = col,
    top_annotation  = top_ann,
    column_split    = col_split,
    column_title_rot = 45,
    column_title_gp  = gpar(fontsize = 9),
    row_names_gp     = gpar(fontsize = 9),
    pct_gp           = gpar(fontsize = 8),
    remove_empty_columns = FALSE,
    remove_empty_rows    = TRUE,
    show_column_names    = FALSE,
    heatmap_legend_param = list(title = "Alteration")
  )

  draw(op, merge_legend = TRUE)
  invisible(op)
}


# ── extract_brca1_snvs ────────────────────────────────────────────────────────
# Pulls SNVs falling in the BRCA1 locus (chr17:41196312–41277500, GRCh37)
# from every sample's CNAqc object and returns one combined data frame.
#
# Note: coordinates are GRCh37-specific. For GRCh38 use chr17:43044295–43170245
# instead.
extract_brca1_snvs <- function(sample_names, fit_dir) {
  brca1_start <- 41196312
  brca1_end   <- 41277500

  lapply(sample_names, function(sample_name) {
    cnaqc <- readRDS(file.path(fit_dir, paste0(sample_name, ".rds")))
    cnaqc$mutations %>%
      dplyr::filter(chr == "chr17",
                    from >= brca1_start,
                    to   <= brca1_end) %>%
      dplyr::select(sample, chr, from, to, ref, alt, consequence_type)
  }) %>%
    dplyr::bind_rows()
}


# ── load_brca1_vep ────────────────────────────────────────────────────────────
# Reads per-sample VEP annotations, maps Ensembl gene IDs to HGNC symbols
# via org.Hs.eg.db, and keeps only the rows annotated as BRCA1.
load_brca1_vep <- function(sample_names, vep_root) {
  lapply(sample_names, function(sample_name) {
    vep_path <- file.path(vep_root, sample_name, "VEP_output.tsv")
    if (!file.exists(vep_path)) {
      message("  VEP file missing for ", sample_name, " — skipping.")
      return(NULL)
    }

    vep_annotations <- read_tsv(vep_path, comment = "##", show_col_types = FALSE)

    ensembl_ids <- unique(vep_annotations$Gene)
    symbols <- mapIds(org.Hs.eg.db,
                      keys      = ensembl_ids,
                      column    = "SYMBOL",
                      keytype   = "ENSEMBL",
                      multiVals = "first")
    symbols <- data.frame(Gene = names(symbols), Gene_symbol = symbols)

    vep_annotations <- left_join(vep_annotations, symbols, by = "Gene")
    vep_annotations %>%
      dplyr::filter(Gene_symbol == "BRCA1") %>%
      dplyr::mutate(sample = sample_name)
  }) %>%
    dplyr::bind_rows()
}


# ==============================================================================
# Main loop: one tumour type at a time
# ==============================================================================

Drivers <- readRDS(drivers_rds)

# Accumulates per-tumour-type enrichment results for any later aggregation.
# Currently left empty for the caller to populate as needed.
cluster_analysis <- dplyr::tibble()

for (tt in ttypes) {
  message("==== Tumour type: ", tt, " ====")

  out_data_tt <- file.path(out_data, tt)
  final_dir   <- file.path(out_data_tt, "bascule")
  cluster_rds <- file.path(final_dir, "sample_clusters.rds")

  if (!file.exists(cluster_rds)) {
    message("  No BASCULE output for ", tt, " — skipping.")
    next
  }

  clusters <- readRDS(cluster_rds)

  # 1. Cluster vs class enrichment (only meaningful with >1 level each).
  if (length(unique(clusters$clusters)) > 1 &&
      length(unique(clusters$class))    > 1) {
    cluster_enrichment_result <- test_cluster_enrichment(
      clusters,
      cluster_col = "clusters",
      class_col   = "class"
    )
  } else {
    message("  Only one cluster or one class — skipping enrichment test.")
  }

  sample_names <- clusters$sample

  # 2. BRCA1 SNVs straight from CNAqc.
  brca1_muts_df <- extract_brca1_snvs(sample_names, fit_dir)

  # 3. VEP-annotated BRCA1 variants joined to cluster assignments.
  vep_annotations_all <- load_brca1_vep(sample_names, vep_root)
  if (nrow(vep_annotations_all) > 0) {
    vep_with_clusters <- vep_annotations_all %>%
      dplyr::left_join(clusters, by = "sample")
  }

  # 4. Driver oncoprint.
  plot_driver_oncoprint(clusters, Drivers)
}
