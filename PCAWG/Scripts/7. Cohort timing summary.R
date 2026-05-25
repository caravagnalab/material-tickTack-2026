library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

# ══════════════════════════════════════════════════════════════════════════════
# COLOUR PALETTES
# ══════════════════════════════════════════════════════════════════════════════

karyotype_palette <- c(
  "2:1"   = "#4e9af1",
  "2:0"   = "#1a5fa8",
  "2:2"   = "#a8d8f0",
  "1:1"   = "#f7c948",
  "1:0"   = "#e07b39",
  "other" = "#cccccc"
)

mutation_palette <- c(
  "CI_M"       = "#d9473f",
  "M"          = "#ff9c57",
  "CNA_driver" = "#7b2d8b",
  "WT"         = "#eeeeee"
)

class_palette <- c(
  "WGD"     = "#7b2d8b",
  "Classic" = "#31a354",
  "HM"      = "#e6550d",
  "other"   = "#bdbdbd"
)

arm_palette <- c(
  "focal"  = "#fc8d59",
  "q_arm"  = "#91bfdb",
  "p_arm"  = "#4575b4",
  "other"  = "#cccccc"
)

clock_rank_palette <- c(
  "1"  = "#2c7bb6",
  "2"  = "#fdae61",
  "3+" = "#d7191c"
)

# ══════════════════════════════════════════════════════════════════════════════
# LOAD DATA
# ══════════════════════════════════════════════════════════════════════════════

Drivers  <- readRDS("~/Documents/material-tickTack-2026/PCAWG/Data/Drivers.rds")
Segments <- readRDS("~/Documents/material-tickTack-2026/PCAWG/Data/Segments.rds")

out_dir            <- "~/Documents/material-tickTack-2026/PCAWG/Plot/"
out_dir_pseudotime <- file.path(out_dir, "drivers+mutation+pseudotime")
out_dir_rank       <- file.path(out_dir, "drivers+mutation+clockrank")
dir.create(out_dir_pseudotime, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_rank,       showWarnings = FALSE, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# PARAMETERS
# ══════════════════════════════════════════════════════════════════════════════

MIN_GENE_PROP <- 0.10

# ══════════════════════════════════════════════════════════════════════════════
# WITHIN-SAMPLE CLOCK RANK  —  derived from Segments, not from Drivers
# ══════════════════════════════════════════════════════════════════════════════

rank_level_order <- c("1", "2", "3+")

segment_event_ranks <- Segments %>%
  distinct(sample, clock_mean) %>%
  filter(!is.na(clock_mean)) %>%
  group_by(sample) %>%
  arrange(clock_mean, .by_group = TRUE) %>%
  mutate(
    clock_rank = if_else(row_number() >= 3, "3+", as.character(row_number()))
  ) %>%
  ungroup() %>%
  select(sample, clock_mean, clock_rank) %>%
  mutate(clock_rank = factor(clock_rank, levels = rank_level_order))

max_clock_rank_per_sample <- segment_event_ranks %>%
  group_by(sample) %>%
  summarise(
    max_clock_rank = rank_level_order[
      max(as.integer(clock_rank), na.rm = TRUE)
    ],
    .groups = "drop"
  ) %>%
  mutate(max_clock_rank = factor(max_clock_rank, levels = rank_level_order))

Drivers <- Drivers %>%
  left_join(segment_event_ranks, by = c("sample", "clock_mean"))

# ══════════════════════════════════════════════════════════════════════════════
# LEGEND LAYOUT HELPER
# ══════════════════════════════════════════════════════════════════════════════

two_row_legend_theme <- theme(
  legend.position    = "bottom",
  legend.box         = "horizontal",
  legend.box.just    = "left",
  legend.spacing.x   = unit(0.3, "cm"),
  legend.spacing.y   = unit(0.15, "cm"),
  legend.key.size    = unit(0.35, "cm"),
  legend.title       = element_text(size = 7, face = "bold"),
  legend.text        = element_text(size = 6),
  legend.margin      = margin(4, 4, 4, 4)
)

LEGEND_NROW <- 2

# ══════════════════════════════════════════════════════════════════════════════
# DEDUPLICATION HELPER
#
# Returns ONE row per gene × sample, keeping the highest-priority record:
#   Priority: CNA_driver > CI_M > M > WT; within ties, timed before untimed.
#
# This is used as the source for ALL marginal plots (karyotype bar, timing
# violin, clock-rank bar, max-rank bar, N-samples bar) so that:
#   - if a sample has 2 mutations on the same gene × segment, only one dot
#     appears in the violin, one bar slice in the karyotype, etc.
#   - the choice of which row to keep is consistent across every panel.
# ══════════════════════════════════════════════════════════════════════════════

dedup_gene_sample <- function(df) {
  df %>%
    group_by(gene, sample) %>%
    arrange(
      match(as.character(mutation_status), c("CNA_driver", "CI_M", "M", "WT")),
      !timed,
      mult_estimate,
      .by_group = TRUE
    ) %>%
    slice(1) %>%
    ungroup()
}

# ══════════════════════════════════════════════════════════════════════════════
# HEATMAP BUILDER
# ══════════════════════════════════════════════════════════════════════════════

build_mutation_heatmap <- function(heatmap_df, heatmap_theme) {
  
  timed_df <- heatmap_df %>% filter(!is.na(timed) & timed == TRUE)
  mult1_df  <- heatmap_df %>% filter(!is.na(mult_class) & mult_class == "mult_1")
  mult2_df  <- heatmap_df %>% filter(!is.na(mult_class) & mult_class == "mult_2")
  
  ggplot(heatmap_df, aes(x = sample, y = gene, fill = mutation_status)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    {if (nrow(timed_df) > 0)
      geom_tile(data = timed_df, aes(x = sample, y = gene),
                fill = NA, colour = "black", linewidth = 0.9,
                inherit.aes = FALSE)} +
    {if (nrow(mult1_df) > 0)
      geom_point(data = mult1_df, aes(x = sample, y = gene),
                 shape = 4, size = 1.2, stroke = 0.5, colour = "white",
                 inherit.aes = FALSE)} +
    {if (nrow(mult2_df) > 0)
      geom_point(data = mult2_df, aes(x = sample, y = gene),
                 shape = 8, size = 1.2, stroke = 0.5, colour = "white",
                 inherit.aes = FALSE)} +
    scale_fill_manual(values = mutation_palette, name = "Mutation",
                      na.value = "grey95", drop = FALSE,
                      guide = guide_legend(nrow = LEGEND_NROW, title.position = "top")) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL, title = "Mutation status") +
    heatmap_theme
}

# ══════════════════════════════════════════════════════════════════════════════
# HEATMAP_DF BUILDER
# Input is already dedup'd; this adds factor levels and annotation columns.
# ══════════════════════════════════════════════════════════════════════════════

build_heatmap_df <- function(drivers_dedup, gene_order, sample_order) {
  drivers_dedup %>%
    mutate(
      gene   = factor(gene,   levels = gene_order),
      sample = factor(sample, levels = sample_order),
      karyotype_class = case_when(
        karyotype == "2:1" ~ "2:1", karyotype == "2:0" ~ "2:0",
        karyotype == "2:2" ~ "2:2", karyotype == "1:1" ~ "1:1",
        karyotype == "1:0" ~ "1:0", TRUE ~ "other"
      ),
      karyotype_class = factor(karyotype_class, levels = names(karyotype_palette)),
      mutation_status = factor(mutation_status, levels = names(mutation_palette)),
      class           = factor(class,           levels = names(class_palette)),
      arm_level_class = factor(
        ifelse(is.na(arm_level_class), "other", arm_level_class),
        levels = names(arm_palette)
      ),
      mult_class = case_when(
        !is.na(mult_estimate) & round(mult_estimate) == 1 ~ "mult_1",
        !is.na(mult_estimate) & round(mult_estimate) >= 2 ~ "mult_2",
        TRUE ~ NA_character_
      )
    )
  # Note: no slice(1) needed here — input is already deduplicated
}

# ══════════════════════════════════════════════════════════════════════════════
# SHARED HELPERS
# ══════════════════════════════════════════════════════════════════════════════

make_meta_strip <- function(data, fill_var, y_label, scale_layer, meta_theme_arg) {
  ggplot(data, aes(x = sample, y = y_label, fill = .data[[fill_var]])) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_layer + scale_x_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL) + meta_theme_arg
}

save_plot <- function(plot, path, n_samples, n_genes, dpi = 180) {
  fig_w     <- max(10, n_samples * 0.28 + 14)
  fig_h     <- max(7,  n_genes   * 0.40 + 8)
  base_font <- max(9, min(16, round(sqrt(fig_w * fig_h) * 0.55)))
  text_override <- theme(
    text             = element_text(size = base_font),
    axis.text.x      = element_blank(), axis.ticks.x = element_blank(),
    axis.text.x.top  = element_text(size = base_font * 0.85),
    axis.text.y      = element_text(size = base_font * 0.85),
    axis.text.y.left = element_text(size = base_font * 0.85, hjust = 1, face = "italic"),
    axis.title       = element_text(size = base_font * 0.90),
    plot.title       = element_text(size = base_font * 1.30, face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = base_font * 0.80, colour = "grey40"),
    legend.text      = element_text(size = base_font * 0.75),
    legend.title     = element_text(size = base_font * 0.80, face = "bold"),
    strip.text       = element_text(size = base_font * 0.85)
  )
  ggsave(plot & text_override, filename = path,
         width = fig_w, height = fig_h, dpi = dpi, limitsize = FALSE)
  message("  Saved -> ", path)
}

# ══════════════════════════════════════════════════════════════════════════════
# SHARED CANCER PALETTE (computed once)
# ══════════════════════════════════════════════════════════════════════════════

cancer_type_list <- unique(na.omit(Drivers$IntoGen_cancer_type))
all_cancer_types <- sort(unique(na.omit(Drivers$IntoGen_cancer_type)))
cancer_palette   <- setNames(
  colorRampPalette(c(
    "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99",
    "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a",
    "#ffff99","#b15928","#8dd3c7","#ffffb3","#bebada"
  ))(length(all_cancer_types)),
  all_cancer_types
)

# ══════════════════════════════════════════════════════════════════════════════
# SHARED PER-LOOP SETUP HELPER
# Filters, deduplicates, builds gene/sample order and shared data frames.
# Returns a named list so both loops stay DRY.
# ══════════════════════════════════════════════════════════════════════════════

prepare_ct <- function(drivers_ct) {
  
  n_total_samples <- n_distinct(drivers_ct$sample)
  
  gene_props <- drivers_ct %>%
    group_by(gene) %>%
    summarise(prop = n_distinct(sample) / n_total_samples, .groups = "drop")
  genes_keep <- gene_props %>% filter(prop >= MIN_GENE_PROP) %>% pull(gene)
  if (length(genes_keep) == 0) return(NULL)
  
  drivers_ct <- drivers_ct %>% filter(gene %in% genes_keep)
  
  # ── One row per gene × sample — used for ALL marginal plots ──────────────
  # This is the single source of truth. Karyotype, timing, clock-rank,
  # and N-samples panels all read from here, so no information is
  # amplified if a gene has multiple mutations or CN events per segment.
  drivers_dedup <- dedup_gene_sample(drivers_ct)
  
  # Gene order: by median clock_mean across deduplicated rows
  gene_order <- drivers_dedup %>%
    filter(!is.na(clock_mean)) %>%
    group_by(gene) %>%
    summarise(median_clock = median(clock_mean, na.rm = TRUE), .groups = "drop") %>%
    arrange(median_clock) %>%
    pull(gene)
  gene_order <- c(setdiff(unique(drivers_dedup$gene), gene_order), gene_order)
  
  sample_order <- drivers_dedup %>%
    distinct(sample, class, wgd_status) %>%
    arrange(class, wgd_status) %>%
    pull(sample) %>%
    unique()
  
  heatmap_df <- build_heatmap_df(drivers_dedup, gene_order, sample_order)
  
  sample_meta <- drivers_dedup %>%
    distinct(sample, class, wgd_status, n_cna, purity,
             ncomponents, best_K, IntoGen_cancer_type, ploidy) %>%
    mutate(sample = factor(sample, levels = sample_order),
           class  = factor(class,  levels = names(class_palette)))
  
  list(
    drivers_ct      = drivers_ct,
    drivers_dedup   = drivers_dedup,
    gene_order      = gene_order,
    sample_order    = sample_order,
    heatmap_df      = heatmap_df,
    sample_meta     = sample_meta,
    n_genes         = length(gene_order),
    n_samples       = length(sample_order),
    n_total_samples = n_total_samples
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# SHARED PANEL BUILDERS (used by both loops)
# ══════════════════════════════════════════════════════════════════════════════

build_kar_bar <- function(heatmap_df, gene_order) {
  # Source: heatmap_df — already one row per gene × sample
  kar_bar_df <- heatmap_df %>%
    filter(!is.na(karyotype_class)) %>%
    count(gene, karyotype_class) %>%
    group_by(gene) %>% mutate(pct = n / sum(n)) %>% ungroup() %>%
    mutate(gene = factor(gene, levels = gene_order))
  
  ggplot(kar_bar_df, aes(x = pct, y = gene, fill = karyotype_class)) +
    geom_col(width = 0.7, position = "stack") +
    scale_fill_manual(values = karyotype_palette, name = "Karyotype", drop = FALSE,
                      guide = guide_legend(nrow = LEGEND_NROW, title.position = "top")) +
    scale_x_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "% samples", y = NULL, title = "Karyotype\ndistribution") +
    theme_minimal(base_size = 9) + two_row_legend_theme +
    theme(axis.text.y  = element_text(size = 8, hjust = 1, face = "italic"),
          axis.text.x  = element_text(size = 7),
          axis.ticks.x = element_line(linewidth = 0.3),
          panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          plot.title   = element_text(size = 8, face = "bold", hjust = 0.5),
          plot.margin  = margin(4, 4, 4, 4))
}

build_timing_violin <- function(drivers_dedup, gene_order) {
  # Source: drivers_dedup — one row per gene × sample.
  # Points coloured by WGD class (Classic / WGD / HM) using class_palette.
  timing_df <- drivers_dedup %>%
    filter(!is.na(clock_mean)) %>%
    mutate(gene  = factor(gene,  levels = gene_order),
           class = factor(class, levels = names(class_palette)))
  
  timing_summary <- timing_df %>%
    group_by(gene) %>%
    summarise(med = median(clock_mean),
              lo  = median(clock_low,  na.rm = TRUE),
              hi  = median(clock_high, na.rm = TRUE), .groups = "drop")
  
  ggplot(timing_df, aes(y = gene, x = clock_mean)) +
    geom_violin(fill = "#b2d8e8", colour = NA, alpha = 0.6,
                scale = "width", trim = TRUE) +
    # Each dot = one sample; colour = WGD class
    geom_point(aes(colour = class), size = 1.2, alpha = 0.8,
               position = position_jitter(height = 0.15, seed = 42)) +
    geom_errorbarh(data = timing_summary,
                   aes(y = gene, xmin = lo, xmax = hi, x = NULL),
                   height = 0.25, colour = "grey30", linewidth = 0.5) +
    geom_point(data = timing_summary, aes(y = gene, x = med),
               shape = 21, fill = "white", colour = "grey20",
               size = 2.5, stroke = 0.8) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                       labels = c("0", "0.5", "1")) +
    scale_colour_manual(values = class_palette, name = "Class",
                        na.value = "grey80", drop = FALSE,
                        guide = guide_legend(nrow = LEGEND_NROW, title.position = "top")) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "Clock mean", y = NULL, title = "Timing\ndistribution") +
    theme_minimal(base_size = 9) + two_row_legend_theme +
    theme(axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x  = element_text(size = 7),
          axis.ticks.x = element_line(linewidth = 0.3),
          panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
          panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
          panel.grid.minor   = element_blank(),
          plot.title  = element_text(size = 8, face = "bold", hjust = 0.5),
          plot.margin = margin(4, 4, 4, 2))
}

build_n_samples_bar <- function(drivers_dedup, gene_order, n_total_samples, ct) {
  # Source: drivers_dedup — one row per gene × sample, so count() = n unique samples
  sample_count_df <- drivers_dedup %>%
    count(gene, name = "n_samples") %>%
    mutate(gene = factor(gene, levels = gene_order), label = as.character(n_samples))
  
  ggplot(sample_count_df, aes(x = n_samples, y = gene)) +
    geom_col(width = 0.7, fill = cancer_palette[ct], colour = NA) +
    geom_text(aes(label = label), hjust = -0.1, size = 2.2, colour = "grey30") +
    scale_x_continuous(name = "N samples", expand = expansion(mult = c(0, 0.35))) +
    scale_y_discrete(drop = FALSE) +
    labs(y = NULL, title = paste0("N samples\n(total = ", n_total_samples, ")")) +
    theme_minimal(base_size = 9) +
    theme(axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x  = element_text(size = 7),
          axis.ticks.x = element_line(linewidth = 0.3),
          panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          plot.title  = element_text(size = 8, face = "bold", hjust = 0.5),
          plot.margin = margin(4, 8, 4, 2))
}

build_meta_stack <- function(sample_meta, meta_theme) {
  make_strip <- function(fill_var, y_label, scale_layer)
    make_meta_strip(sample_meta, fill_var, y_label, scale_layer, meta_theme)
  
  p_class  <- make_strip("class", "Class",
                         scale_fill_manual(values = class_palette, name = "Class", na.value = "grey95",
                                           drop = FALSE, guide = guide_legend(nrow = LEGEND_NROW, title.position = "top")))
  p_ncna   <- make_strip("n_cna", "N CNA",
                         scale_fill_gradientn(colours = c("#fee8c8","#e34a33"), name = "N CNA", na.value = "grey95",
                                              guide = guide_colourbar(title.position = "top", barwidth = unit(1.5,"cm"), barheight = unit(0.3,"cm"))))
  p_purity <- make_strip("purity", "Purity",
                         scale_fill_gradientn(colours = c("#edf8e9","#238b45"), name = "Purity", na.value = "grey95",
                                              limits = c(0,1),
                                              guide = guide_colourbar(title.position = "top", barwidth = unit(1.5,"cm"), barheight = unit(0.3,"cm"))))
  p_ploidy <- make_strip("ploidy", "Ploidy",
                         scale_fill_gradientn(colours = c("#f0f0f0","#636363"), name = "Ploidy", na.value = "grey95",
                                              guide = guide_colourbar(title.position = "top", barwidth = unit(1.5,"cm"), barheight = unit(0.3,"cm"))))
  p_ncomp  <- make_strip("ncomponents", "N comp.",
                         scale_fill_gradientn(colours = c("#ffffd4","#993404"), name = "N comp.", na.value = "grey95",
                                              guide = guide_colourbar(title.position = "top", barwidth = unit(1.5,"cm"), barheight = unit(0.3,"cm"))))
  p_bestk  <- make_strip("best_K", "Best K",
                         scale_fill_gradientn(colours = c("#f7fbff","#2171b5"), name = "Best K", na.value = "grey95",
                                              guide = guide_colourbar(title.position = "top", barwidth = unit(1.5,"cm"), barheight = unit(0.3,"cm"))))
  
  wrap_plots(p_class, p_ncna, p_purity, p_ploidy, p_ncomp, p_bestk, ncol = 1) +
    plot_layout(heights = rep(1, 6))
}

# ══════════════════════════════════════════════════════════════════════════════
# LOOP 1 — PSEUDOTIME VIOLIN PLOTS
# ══════════════════════════════════════════════════════════════════════════════

for (ct in cancer_type_list) {
  
  message("\n── Processing (pseudotime): ", ct, " ──")
  
  prep <- prepare_ct(Drivers %>% filter(IntoGen_cancer_type == ct))
  if (is.null(prep)) { message("  No genes pass filter, skipping."); next }
  
  message("  Genes: ", prep$n_genes, "  |  Samples: ", prep$n_samples)
  
  heatmap_theme <- theme_minimal(base_size = 9) + two_row_legend_theme +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid  = element_blank(),
          plot.title  = element_text(size = 8, face = "bold", hjust = 0.5),
          plot.margin = margin(4, 2, 4, 2))
  
  meta_theme <- theme_minimal(base_size = 9) + two_row_legend_theme +
    theme(axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y  = element_text(size = 7, hjust = 1),
          axis.ticks.y = element_blank(), panel.grid = element_blank(),
          plot.margin  = margin(1, 2, 1, 2))
  
  p_kar_bar   <- build_kar_bar(prep$heatmap_df, prep$gene_order)
  p_timing    <- build_timing_violin(prep$drivers_dedup, prep$gene_order)
  p_mut_heat  <- build_mutation_heatmap(prep$heatmap_df, heatmap_theme)
  p_n_samples <- build_n_samples_bar(prep$drivers_dedup, prep$gene_order,
                                     prep$n_total_samples, ct)
  meta_stack  <- build_meta_stack(prep$sample_meta, meta_theme)
  
  layout_design <- "
  ABCD
  ##E#
"
  final_plot <- (p_kar_bar + p_timing + p_mut_heat + p_n_samples + meta_stack) +
    plot_layout(design = layout_design, widths = c(0.1, 0.1, 0.3, 0.1),
                heights = c(4, 1), guides = "collect") +
    plot_annotation(
      title    = paste0("Driver CNA annotation summary — ", ct),
      subtitle = paste0(
        prep$n_genes, " genes (\u2265", MIN_GENE_PROP * 100, "% of ",
        prep$n_samples, " samples)  |  ",
        "Each dot = one sample, coloured by WGD class  |  ",
        "Heatmap: black border = timed | X = mult 1 | * = mult 2"
      ),
      theme = theme(plot.title    = element_text(size = 13, face = "bold"),
                    plot.subtitle = element_text(size = 8,  colour = "grey40"))
    ) & two_row_legend_theme
  
  save_plot(final_plot,
            file.path(out_dir_pseudotime, paste0("drivers_summary_", ct, ".png")),
            prep$n_samples, prep$n_genes)
}

# ══════════════════════════════════════════════════════════════════════════════
# LOOP 2 — CLOCK RANK PLOTS
# ══════════════════════════════════════════════════════════════════════════════

for (ct in cancer_type_list) {
  
  message("\n── Processing (clock rank): ", ct, " ──")
  
  prep <- prepare_ct(Drivers %>% filter(IntoGen_cancer_type == ct))
  if (is.null(prep)) { message("  No genes pass filter, skipping."); next }
  
  heatmap_theme <- theme_minimal(base_size = 9) + two_row_legend_theme +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid  = element_blank(),
          plot.title  = element_text(size = 8, face = "bold", hjust = 0.5),
          plot.margin = margin(4, 2, 4, 2))
  
  meta_theme <- theme_minimal(base_size = 9) + two_row_legend_theme +
    theme(axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y  = element_text(size = 7, hjust = 1),
          axis.ticks.y = element_blank(), panel.grid = element_blank(),
          plot.margin  = margin(1, 2, 1, 2))
  
  # Panel A: clock_rank distribution per gene
  # Source: drivers_dedup — one row per gene × sample
  clock_rank_bar_df <- prep$drivers_dedup %>%
    filter(!is.na(clock_rank)) %>%
    mutate(gene = factor(gene, levels = prep$gene_order)) %>%
    count(gene, clock_rank) %>%
    group_by(gene) %>% mutate(pct = n / sum(n)) %>% ungroup()
  
  p_clock_rank_bar <- ggplot(clock_rank_bar_df,
                             aes(x = pct, y = gene, fill = clock_rank)) +
    geom_col(width = 0.7, position = "stack") +
    scale_fill_manual(values = clock_rank_palette, name = "Clock rank", drop = FALSE,
                      guide = guide_legend(nrow = LEGEND_NROW, title.position = "top")) +
    scale_x_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "% samples", y = NULL, title = "Within-sample\nclock rank") +
    theme_minimal(base_size = 9) + two_row_legend_theme +
    theme(axis.text.y  = element_text(size = 8, hjust = 1, face = "italic"),
          axis.text.x  = element_text(size = 7),
          axis.ticks.x = element_line(linewidth = 0.3),
          panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          plot.title  = element_text(size = 8, face = "bold", hjust = 0.5),
          plot.margin = margin(4, 4, 4, 4))
  
  # Panel B: max clock_rank of host sample per gene
  # Source: drivers_dedup joined to max_clock_rank_per_sample
  max_rank_bar_df <- prep$drivers_dedup %>%
    filter(!is.na(clock_rank)) %>%
    select(gene, sample) %>%
    left_join(max_clock_rank_per_sample, by = "sample") %>%
    filter(!is.na(max_clock_rank)) %>%
    mutate(gene = factor(gene, levels = prep$gene_order)) %>%
    count(gene, max_clock_rank) %>%
    group_by(gene) %>% mutate(pct = n / sum(n)) %>% ungroup()
  
  p_max_rank_bar <- ggplot(max_rank_bar_df,
                           aes(x = pct, y = gene, fill = max_clock_rank)) +
    geom_col(width = 0.7, position = "stack") +
    scale_fill_manual(values = clock_rank_palette, name = "Max rank\nin sample",
                      drop = FALSE,
                      guide = guide_legend(nrow = LEGEND_NROW, title.position = "top")) +
    scale_x_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "% samples", y = NULL, title = "Max clock rank\nof host sample") +
    theme_minimal(base_size = 9) + two_row_legend_theme +
    theme(axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x  = element_text(size = 7),
          axis.ticks.x = element_line(linewidth = 0.3),
          panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          plot.title  = element_text(size = 8, face = "bold", hjust = 0.5),
          plot.margin = margin(4, 4, 4, 2))
  
  p_mut_heat  <- build_mutation_heatmap(prep$heatmap_df, heatmap_theme)
  p_n_samples <- build_n_samples_bar(prep$drivers_dedup, prep$gene_order,
                                     prep$n_total_samples, ct)
  meta_stack  <- build_meta_stack(prep$sample_meta, meta_theme)
  
  layout_design_rank <- "
  ABCD
  ##E#
"
  final_plot_rank <- (
    p_clock_rank_bar + p_max_rank_bar + p_mut_heat + p_n_samples + meta_stack
  ) +
    plot_layout(design = layout_design_rank,
                widths = c(0.13, 0.10, 0.30, 0.10),
                heights = c(4, 1), guides = "collect") +
    plot_annotation(
      title    = paste0("Driver CNA — clock rank summary — ", ct),
      subtitle = paste0(
        prep$n_genes, " genes (\u2265", MIN_GENE_PROP * 100, "% of ",
        prep$n_samples, " samples)  |  ",
        "Heatmap: black border = timed | X = mult 1 | * = mult 2  |  ",
        "Clock rank derived from Segments.rds"
      ),
      theme = theme(plot.title    = element_text(size = 13, face = "bold"),
                    plot.subtitle = element_text(size = 8,  colour = "grey40"))
    ) & two_row_legend_theme
  
  save_plot(final_plot_rank,
            file.path(out_dir_rank, paste0("drivers_clockrank_", ct, ".png")),
            prep$n_samples, prep$n_genes)
}

message("\n✔  All done.")
