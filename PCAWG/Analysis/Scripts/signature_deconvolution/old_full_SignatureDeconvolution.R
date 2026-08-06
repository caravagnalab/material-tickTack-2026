
library(dplyr)
library(CNAqc)
library(tidyverse)

# ── Read arguments from the HPC job array ─────────────────────────────────────
# ttype_idx:    job array index — read from SLURM_ARRAY_TASK_ID (or SGE_TASK_ID)
# analysis_idx: first positional argument passed to Rscript
#
# Example submission:
#   sbatch --array=1-30 run.sh 0     # matrix generation for all tumour types
#   sbatch --array=1-30 run.sh 1     # BASCULE inference for all tumour types
#
# Local testing:
#   Rscript script.R 1              (ttype_idx taken from env, analysis_idx = 1)

ttype_idx <- as.integer(
  Sys.getenv("SLURM_ARRAY_TASK_ID",    # SLURM
             unset = Sys.getenv("SGE_TASK_ID",    # SGE / PBS fallback
                                unset = NA_character_))
)

if (is.na(ttype_idx)) {
  stop("Could not read job array index: ",
       "SLURM_ARRAY_TASK_ID and SGE_TASK_ID are both unset. ",
       "Submit this script as a job array.")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Missing required argument: analysis_idx. ",
       "Usage: Rscript script.R <analysis_idx>  (0 = matrix generation, 1 = BASCULE)")
}
analysis_idx <- as.integer(args[1])

if (!analysis_idx %in% c(0L, 1L)) {
  stop("analysis_idx must be 0 (matrix generation) or 1 (BASCULE inference), got: ", analysis_idx)
}

message("ttype_idx: ", ttype_idx, " | analysis_idx: ", analysis_idx)

# ── Setup ──────────────────────────────────────────────────────────────────────
samples_info <- readRDS("/orfeo/cephfs/scratch/cdslab/gsantacatterina/material-tickTack-2026/PCAWG/Data/Samples.rds")
out_data     <- "/orfeo/cephfs/scratch/cdslab/gsantacatterina/material-tickTack-2026/PCAWG/Analysis_results/Signatures/"
dir.create(out_data, recursive = TRUE)

ttypes  <- unique(samples_info$ttype)
tt      <- ttypes[ttype_idx]
project <- "PCAWG"

message("Tumour type: ", tt)

samples_info <- samples_info %>%
  dplyr::filter(ttype == tt) %>%
  dplyr::select(sample, ttype, class) %>%
  dplyr::group_by(class) %>%
  dplyr::mutate(sampleId = paste0("S_", class, dplyr::row_number())) %>%
  dplyr::ungroup()

# ── analysis_idx == 0 : SigProfiler matrix generation ─────────────────────────
if (analysis_idx == 0L) {
  reticulate::use_python("/orfeo/scratch/area/lvaleriani/myconda/bin")
  reticulate::py_config()
  library(SigProfilerMatrixGeneratorR)
  
  get_sig_profiler_table <- function(i) {
    print(i)
    sample_name <- samples_info$sample[i]
    sample_id   <- samples_info$sampleId[i]
    
    cnaqc  <- readRDS(paste0(
      "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/",
      "inference_results_5ncomponents/", sample_name, ".rds"
    ))
    genome <- cnaqc$reference_genome
    cnaqc$mutations$sample <- sample_id
    
    cnaqc$mutations %>%
      dplyr::select(chr, from, to, ref, alt, sample) %>%
      dplyr::distinct() %>%
      tidyr::separate(col = chr, sep = "chr", into = c("tmp", "chr")) %>%
      dplyr::select(-tmp) %>%
      dplyr::mutate(
        Project  = project,
        Genome   = genome,
        mut_type = "SNP",
        Type     = "SOMATIC",
        ID       = sample,
        Sample   = sample
      ) %>%
      dplyr::rename(chrom = chr, pos_start = from, pos_end = to) %>%
      dplyr::select(Project, Sample, ID, Genome, mut_type,
                    chrom, pos_start, pos_end, ref, alt, Type) %>%
      dplyr::filter(ref != alt) %>%
      dplyr::distinct()
  }
  
  sigprofiler_table <- lapply(seq_len(nrow(samples_info)), get_sig_profiler_table) %>%
    dplyr::bind_rows()
  
  sigprofiler_table %>% dplyr::filter(Sample == "S_Classic127") %>% view()
  
  dir.create(out_data_tt, recursive = TRUE)
  
  write.table(sigprofiler_table,
              file  = file.path(out_data_tt, paste0(project, ".txt")),
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  SigProfilerMatrixGeneratorR(
    project     = project,
    genome      = "GRCh37",
    matrix_path = out_data_tt,
    plot        = FALSE
  )
}

# ── analysis_idx == 1 : BASCULE inference ─────────────────────────────────────
if (analysis_idx == 1L) {
  reticulate::use_condaenv("/orfeo/scratch/cdslab/ggandolfi/miniconda/envs/bascule-env")
  py <- reticulate::import("pybascule")
  library(bascule)
  
  out_data_tt <- file.path(out_data, tt)
  data_sbs    <- file.path(out_data_tt, "output", "SBS", paste0(project, ".SBS96.all"))
  
  sig_matrix <- read.table(data_sbs, header = TRUE)
  rownames(sig_matrix) <- sig_matrix$MutationType
  sig_counts <- sig_matrix %>%
    dplyr::select(-MutationType) %>%
    t() %>%
    as.data.frame()
  
  sig_counts <- sig_counts[rowSums(sig_counts) > 0, ]
  
  cat   <- list(SBS = bascule::COSMIC_sbs %>% as.data.frame())
  input <- list(SBS = sig_counts)
  
  x <- bascule::fit(
    counts           = input,
    k_list           = 0,
    reference_cat    = cat,
    keep_sigs        = c("SBS1", "SBS5"),
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
  
  x_refined         <- refine_denovo_signatures(x)
  x_refined_cluster <- fit_clustering(x_refined, cluster = 3)
  x_refined_cluster <- merge_clusters(x_refined_cluster)
  
  plot_exp         <- plot_exposures(x = x,                 sample_name = TRUE)
  plot_exp_refined <- plot_exposures(x = x_refined_cluster, sample_name = TRUE)
  
  sample_info_with_clusters <- x_refined_cluster$clustering$clusters %>%
    dplyr::rename(sampleId = samples) %>%
    dplyr::left_join(samples_info, by = "sampleId")
  
  final_dir <- file.path(out_data_tt, "bascule")
  dir.create(final_dir, recursive = TRUE)
  
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
}




if (analysis_idx == 2) {
  
  library(tidyverse)
  
  test_cluster_enrichment <- function(df, cluster_col = "cluster", class_col = "class", alpha = 0.05) {
    
    contingency <- table(df[[cluster_col]], df[[class_col]])
    chi_result  <- chisq.test(contingency)
    
    # --- Step 1: Global homogeneity test ---
    cat(sprintf("Chi-squared test of homogeneity: X²=%.3f, dof=%d, p=%.4e\n",
                chi_result$statistic, chi_result$parameter, chi_result$p.value))
    if (chi_result$p.value >= alpha) {
      cat("  → No evidence that class incidence varies across clusters.\n")
      return(invisible(list(chi2_test = chi_result)))
    }
    cat("  → Class incidence is NOT constant across clusters. Running post-hoc analysis...\n\n")
    
    # --- Step 2: Standardized Pearson residuals ---
    # residual = (observed - expected) / sqrt(expected)
    # |residual| > 2 ~ nominally significant, > 3.29 ~ Bonferroni-corrected at 0.05
    residuals <- chi_result$stdres  # stdres are preferred over simple Pearson residuals
    
    results <- as.data.frame(as.table(residuals)) %>%
      setNames(c("cluster", "class", "std_residual")) %>%
      dplyr::mutate(
        observed        = as.vector(contingency),
        expected        = as.vector(chi_result$expected),
        fold_enrichment = observed / expected,
        enriched        = std_residual > 0,
        # Two-sided p-value from std residual (approx. standard normal)
        p_value         = 2 * pnorm(abs(std_residual), lower.tail = FALSE),
        p_adj           = p.adjust(p_value, method = "BH"),
        significant     = p_adj < alpha
      ) %>%
      dplyr::arrange(desc(std_residual))
    
    # --- Print enriched cells ---
    cat("Significant enrichments/depletions (FDR corrected):\n")
    results %>%
      dplyr::filter(significant) %>%
      dplyr::mutate(direction = ifelse(enriched, "enriched", "depleted")) %>%
      dplyr::select(cluster, class, direction, observed, expected, fold_enrichment, std_residual, p_adj) %>%
      print()
    
    invisible(list(chi2_test = chi_result, residuals = results))
  }
  
  library(ComplexHeatmap)
  
  plot_driver_oncoprint <- function(clusters, Drivers, sample_col = "sample", 
                                    cluster_col = "clusters", class_col = "class") {
    
    # --- 1. Join clusters/class info onto drivers ---
    drivers_annotated <- Drivers %>%
      dplyr::filter(mutation_status != "WT") %>%
      inner_join(clusters %>% dplyr::select(all_of(c(sample_col, cluster_col, class_col))),
                 by = c("sample_id" = sample_col))
    
    samples  <- unique(clusters[[sample_col]])
    genes    <- unique(drivers_annotated$gene)
    
    # --- 2. Build the alteration matrix (genes x samples) ---
    # Each cell can have multiple alterations — ComplexHeatmap expects
    # a character matrix where multiple hits are separated by ";"
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
        mat[g, s] <- paste(mat[g, s], alt, sep = ";")  # multi-hit
      }
    }
    
    # --- 3. Sort samples by class then cluster ---
    sample_meta <- clusters %>%
      dplyr::select(all_of(c(sample_col, cluster_col, class_col))) %>%
      dplyr::distinct() %>%
      dplyr::arrange(.data[[class_col]], .data[[cluster_col]])
    
    sample_order <- sample_meta[[sample_col]]
    sample_order <- sample_order[sample_order %in% colnames(mat)]
    mat <- mat[, sample_order, drop = FALSE]
    
    # Keep only genes mutated in at least 1 sample
    mat <- mat[rowSums(mat != "") > 0, , drop = FALSE]
    
    # Sort genes by mutation frequency
    gene_order <- order(rowSums(mat != ""), decreasing = TRUE)
    mat <- mat[gene_order, , drop = FALSE]
    
    # --- 4. Alteration colours ---
    col <- c(
      "M"          = "#E64B35",  # red     - point mutation
      "CI_M"       = "#F39B7F",  # salmon  - clonal/subclonal mutation
      "CNA_driver" = "#4DBBD5"   # blue    - copy number
    )
    
    alter_fun <- list(
      background = function(x, y, w, h)
        grid.rect(x, y, w * 0.9, h * 0.9,
                  gp = gpar(fill = "#F5F5F5", col = NA)),
      M = function(x, y, w, h)
        grid.rect(x, y, w * 0.9, h * 0.9,
                  gp = gpar(fill = col["M"], col = NA)),
      CI_M = function(x, y, w, h)
        grid.rect(x, y, w * 0.9, h * 0.45,   # half-height bar
                  gp = gpar(fill = col["CI_M"], col = NA)),
      CNA_driver = function(x, y, w, h)
        grid.segments(x - w * 0.4, y, x + w * 0.4, y,
                      gp = gpar(col = col["CNA_driver"], lwd = 3))
    )
    
    # --- 5. Top annotations: cluster & class ---
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
    
    # --- 6. Column splits: divide by class then cluster ---
    col_split <- factor(
      paste(ann_df[[class_col]], ann_df[[cluster_col]], sep = " | "),
      levels = unique(paste(ann_df[[class_col]], ann_df[[cluster_col]], sep = " | "))
    )
    
    # --- 7. Draw oncoprint ---
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
  
  
  # tt = "PRAD"
  # tt = "PACA"
  cluster_analysis = dplyr::tibble()  
  for (tt in ttypes) {
    print(tt)
    out_data_tt <- file.path(out_data, tt)
    final_dir <- file.path(out_data_tt, "bascule")  
    
    clusters = readRDS(file.path(final_dir, "sample_clusters.rds"))
    if (length(unique(clusters$clusters)) > 1 & length(unique(clusters$class)) > 1) {
      cluster_enrichment_result = test_cluster_enrichment(clusters, cluster_col = "clusters", class_col = "class")
    }
    clusters %>% dplyr::filter(clusters == "G1")
    
    sample_names = clusters$sample
    sample_name = sample_names[2]

    Drivers <- readRDS("/orfeo/cephfs/scratch/cdslab/gsantacatterina/material-tickTack-2026/PCAWG/Data/Drivers.rds")
    unique(Drivers$mutation_status)

    library(org.Hs.eg.db)
    library(AnnotationDbi)

    BRCA1_muts_df = lapply(sample_names, function(sample_name) {
      print(sample_name)
      cnaqc  <- readRDS(paste0(
        "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/",
        "inference_results_5ncomponents/", sample_name, ".rds"
      ))

      cnaqc$mutations %>% dplyr::filter(chr == "chr17", from >= 41196312, to <= 41277500) %>%
        dplyr::select(sample, chr, from, to, ref, alt, consequence_type)
    }) %>% do.call("bind_rows", .)

    sample_name

    vep_annotations_all = lapply(sample_names, function(sample_name) {
      vep_annotations = read_tsv(file.path(
        "/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack/PCAWG/selection/clonal_analysis_PCAWG/",
        sample_name,
        "/VEP_output.tsv"
      ), comment = "##")

      ensembl_ids <- vep_annotations$Gene %>% unique()

      symbols <- mapIds(org.Hs.eg.db,
                        keys    = ensembl_ids,
                        column  = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")

      symbols = data.frame(Gene = names(symbols), Gene_symbol=symbols)
      vep_annotations = left_join(vep_annotations, symbols, by="Gene")
      vep_annotations %>% dplyr::filter(Gene_symbol == "BRCA1") %>% dplyr::mutate(sample = sample_name)
    }) %>% do.call("bind_rows", .)

    vep_annotations_all %>% dplyr::left_join(clusters) %>% view()



    BRCA1_muts_df %>% dplyr::left_join(clusters)


    plot_driver_oncoprint(clusters, Drivers)
    
    
  }  
}
