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

out_dir      <- "~/Documents/material-tickTack-2026/PCAWG/Plot/"
out_dir_pseudotime <- file.path(out_dir, "drivers+mutation+pseudotime")
out_dir_rank <- file.path(out_dir, "drivers+mutation+clockrank")
dir.create(out_dir_rank, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_pseudotime, showWarnings = FALSE, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# PARAMETERS
# ══════════════════════════════════════════════════════════════════════════════

MIN_GENE_PROP <- 0.30

# ══════════════════════════════════════════════════════════════════════════════
# WITHIN-SAMPLE CLOCK RANK  —  derived from Segments, not from Drivers
#
# Segments contains all timed CNA events for every sample. We rank the
# distinct clock_mean values within each sample (earliest = 1) and
# collapse everything at rank ≥ 3 into "3+".
# The resulting lookup (sample × clock_mean → clock_rank) is then joined to
# Drivers so that each driver row inherits the rank of its host event.
#
# max_clock_rank_per_sample is computed BEFORE the join so it reflects the
# full CNA complexity of the sample, not just the driver-bearing events.
# ══════════════════════════════════════════════════════════════════════════════

rank_level_order <- c("1", "2", "3+")

# Step 1 ── rank events within each sample using Segments.
# event_rank is the position within the sample (1, 2, 3 …); cap at "3+".
# No case_when needed: row_number() already gives the numeric rank directly.
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

# Step 2 ── max rank per sample (based on ALL segments, before driver join).
# as.integer() on an ordered factor gives position in rank_level_order,
# so max() correctly identifies the highest rank reached by each sample.
max_clock_rank_per_sample <- segment_event_ranks %>%
  group_by(sample) %>%
  summarise(
    max_clock_rank = rank_level_order[
      max(as.integer(clock_rank), na.rm = TRUE)
    ],
    .groups = "drop"
  ) %>%
  mutate(max_clock_rank = factor(max_clock_rank, levels = rank_level_order))

# Step 3 ── join clock_rank into Drivers
# Only samples present in Segments and with a matching clock_mean get a rank.
# Drivers rows with clock_mean = NA (untimed) keep clock_rank = NA.
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
# HEATMAP BUILDER  —  shared between both loops
#
# The mutation heatmap now carries two overlaid annotation layers:
#
#  1. BORDER COLOUR  — whether the CNA carrying the mutation was timed
#     timed == TRUE   →  solid coloured border  (black)
#     timed == FALSE  →  no border / white border
#
#  2. MULTIPLICITY PATTERN  — mult_estimate (if not NA)
#     mult == 1  →  diagonal hatching (simulated with small rotated squares)
#     mult == 2  →  cross-hatching    (simulated with two sets of squares)
#     NA         →  no pattern
#
# ggplot2 does not natively support fill patterns, so the pattern is rendered
# as a second geom_point layer using unicode characters or shape 4/8 symbols
# drawn on top of the tile.  This keeps the code base-R friendly with no
# extra packages.
# ══════════════════════════════════════════════════════════════════════════════

build_mutation_heatmap <- function(heatmap_df, heatmap_theme) {
  
  # Separate layers by annotation type
  timed_df <- heatmap_df %>% filter(!is.na(timed) & timed == TRUE)
  
  mult1_df  <- heatmap_df %>%
    filter(!is.na(mult_class) & mult_class == "mult_1")
  mult2_df  <- heatmap_df %>%
    filter(!is.na(mult_class) & mult_class == "mult_2")
  
  p <- ggplot(heatmap_df,
              aes(x = sample, y = gene, fill = mutation_status)) +
    
    # Base tile (fill = mutation status)
    geom_tile(colour = "white", linewidth = 0.2) +
    
    # Timed border: thick coloured outline drawn as an additional tile layer
    # with fill = NA so only the border shows
    {if (nrow(timed_df) > 0)
      geom_tile(data = timed_df,
                aes(x = sample, y = gene),
                fill      = NA,
                colour    = "black",
                linewidth = 0.9,
                inherit.aes = FALSE)
    } +
    
    # Multiplicity = 1 → single diagonal cross (shape 4, X mark)
    {if (nrow(mult1_df) > 0)
      geom_point(data = mult1_df,
                 aes(x = sample, y = gene),
                 shape       = 4,          # "X"
                 size        = 1.2,
                 stroke      = 0.5,
                 colour      = "white",
                 inherit.aes = FALSE)
    } +
    
    # Multiplicity = 2 → asterisk / star (shape 8)
    {if (nrow(mult2_df) > 0)
      geom_point(data = mult2_df,
                 aes(x = sample, y = gene),
                 shape       = 8,          # "*"
                 size        = 1.2,
                 stroke      = 0.5,
                 colour      = "white",
                 inherit.aes = FALSE)
    } +
    
    scale_fill_manual(values = mutation_palette, name = "Mutation",
                      na.value = "grey95", drop = FALSE,
                      guide = guide_legend(nrow = LEGEND_NROW,
                                           title.position = "top")) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL, title = "Mutation status") +
    heatmap_theme
  p
}

# ══════════════════════════════════════════════════════════════════════════════
# SHARED HEATMAP_DF BUILDER
# Prepares heatmap_df with timed and mult_class columns needed by the builder
# ══════════════════════════════════════════════════════════════════════════════

build_heatmap_df <- function(drivers_ct, gene_order, sample_order) {
  drivers_ct %>%
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
      # mult_class: classify multiplicity for pattern overlay
      mult_class = case_when(
        !is.na(mult_estimate) & round(mult_estimate) == 1 ~ "mult_1",
        !is.na(mult_estimate) & round(mult_estimate) >= 2 ~ "mult_2",
        TRUE ~ NA_character_
      )
    ) %>%
    group_by(gene, sample) %>%
    # Priority: CNA_driver > CI_M > M > WT; within ties keep timed=TRUE first
    arrange(
      match(as.character(mutation_status), c("CNA_driver","CI_M","M","WT")),
      !timed,            # TRUE (timed) before FALSE
      mult_estimate,     # lower mult first (will be overridden by user if needed)
      .by_group = TRUE
    ) %>%
    slice(1) %>%
    ungroup()
}

# ══════════════════════════════════════════════════════════════════════════════
# SHARED HELPERS
# ══════════════════════════════════════════════════════════════════════════════

make_meta_strip <- function(data, fill_var, y_label, scale_layer,
                            meta_theme_arg) {
  ggplot(data, aes(x = sample, y = y_label, fill = .data[[fill_var]])) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_layer + scale_x_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL) + meta_theme_arg
}

save_plot <- function(plot, path, n_samples, n_genes, dpi = 180) {
  fig_w <- max(10, n_samples * 0.28 + 14)
  fig_h <- max(7,  n_genes   * 0.40 + 8)
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
# CANCER TYPE LOOP — ORIGINAL PLOTS  (timing violin)
# ══════════════════════════════════════════════════════════════════════════════

cancer_type_list <- unique(na.omit(Drivers$IntoGen_cancer_type))

# Shared cancer palette (computed once)
all_cancer_types <- sort(unique(na.omit(Drivers$IntoGen_cancer_type)))
cancer_palette   <- setNames(
  colorRampPalette(c(
    "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99",
    "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a",
    "#ffff99","#b15928","#8dd3c7","#ffffb3","#bebada"
  ))(length(all_cancer_types)),
  all_cancer_types
)

for (ct in cancer_type_list) {
  
  message("\n── Processing (original): ", ct, " ──")
  
  drivers_ct <- Drivers %>% filter(IntoGen_cancer_type == ct)
  if (nrow(drivers_ct) == 0) { message("  No data, skipping."); next }
  
  n_total_samples <- n_distinct(drivers_ct$sample)
  gene_props <- drivers_ct %>%
    group_by(gene) %>%
    summarise(prop = n_distinct(sample) / n_total_samples, .groups = "drop")
  genes_keep <- gene_props %>% filter(prop >= MIN_GENE_PROP) %>% pull(gene)
  if (length(genes_keep) == 0) { message("  No genes pass filter, skipping."); next }
  
  drivers_ct <- drivers_ct %>% filter(gene %in% genes_keep)
  message("  Genes retained: ", length(genes_keep))
  
  gene_order <- drivers_ct %>%
    filter(!is.na(clock_mean)) %>%
    group_by(gene) %>%
    summarise(median_clock = median(clock_mean, na.rm = TRUE), .groups = "drop") %>%
    arrange(median_clock) %>% pull(gene)
  gene_order <- c(setdiff(unique(drivers_ct$gene), gene_order), gene_order)
  
  sample_order <- drivers_ct %>%
    distinct(sample, class, wgd_status) %>%
    arrange(class, wgd_status) %>% pull(sample) %>% unique()
  
  n_genes   <- length(gene_order)
  n_samples <- length(sample_order)
  
  heatmap_df <- build_heatmap_df(drivers_ct, gene_order, sample_order)
  
  sample_meta <- drivers_ct %>%
    distinct(sample, class, wgd_status, n_cna, purity,
             ncomponents, best_K, IntoGen_cancer_type, ploidy) %>%
    mutate(sample = factor(sample, levels = sample_order),
           class  = factor(class,  levels = names(class_palette)))
  
  heatmap_theme <- theme_minimal(base_size = 9) + two_row_legend_theme +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid  = element_blank(),
          plot.title  = element_text(size = 8, face = "bold", hjust = 0.5),
          plot.margin = margin(4, 2, 4, 2))
  
  meta_theme <- theme_minimal(base_size = 9) + two_row_legend_theme +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 7, hjust = 1),
          axis.ticks.y = element_blank(), panel.grid = element_blank(),
          plot.margin = margin(1, 2, 1, 2))
  
  # Karyotype bar
  kar_bar_df <- heatmap_df %>%
    filter(!is.na(karyotype_class)) %>%
    count(gene, karyotype_class) %>%
    group_by(gene) %>% mutate(pct = n / sum(n)) %>% ungroup() %>%
    mutate(gene = factor(gene, levels = gene_order))
  
  p_kar_bar <- ggplot(kar_bar_df, aes(x = pct, y = gene, fill = karyotype_class)) +
    geom_col(width = 0.7, position = "stack") +
    scale_fill_manual(values = karyotype_palette, name = "Karyotype", drop = FALSE,
                      guide = guide_legend(nrow = LEGEND_NROW, title.position = "top")) +
    scale_x_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "% samples", y = NULL, title = "Karyotype\ndistribution") +
    theme_minimal(base_size = 9) + two_row_legend_theme +
    theme(axis.text.y = element_text(size = 8, hjust = 1, face = "italic"),
          axis.text.x = element_text(size = 7),
          axis.ticks.x = element_line(linewidth = 0.3),
          panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
          plot.margin = margin(4, 4, 4, 4))
  
  # Timing violin
  timing_df <- drivers_ct %>%
    filter(!is.na(clock_mean)) %>%
    mutate(gene = factor(gene, levels = gene_order))
  timing_summary <- timing_df %>%
    group_by(gene) %>%
    summarise(med = median(clock_mean),
              lo  = median(clock_low,  na.rm = TRUE),
              hi  = median(clock_high, na.rm = TRUE), .groups = "drop")
  
  p_timing <- ggplot(timing_df, aes(y = gene, x = clock_mean)) +
    geom_violin(fill = "#b2d8e8", colour = NA, alpha = 0.6,
                scale = "width", trim = TRUE) +
    geom_point(aes(colour = clock_mean), size = 1.2, alpha = 0.7,
               position = position_jitter(height = 0.15, seed = 42)) +
    geom_errorbarh(data = timing_summary,
                   aes(y = gene, xmin = lo, xmax = hi, x = NULL),
                   height = 0.25, colour = "grey30", linewidth = 0.5) +
    geom_point(data = timing_summary, aes(y = gene, x = med),
               shape = 21, fill = "white", colour = "grey20", size = 2.5, stroke = 0.8) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                       labels = c("0", "0.5", "1")) +
    scale_colour_gradientn(colours = c("#ffffcc","#41b6c4","#253494"),
                           limits = c(0, 1), guide = "none") +
    scale_y_discrete(drop = FALSE) +
    labs(x = "Clock mean", y = NULL, title = "Timing\ndistribution") +
    theme_minimal(base_size = 9) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 7),
          axis.ticks.x = element_line(linewidth = 0.3),
          panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
          panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
          plot.margin = margin(4, 4, 4, 2))
  
  # N samples bar
  sample_count_df <- drivers_ct %>%
    distinct(gene, sample) %>%
    count(gene, name = "n_samples") %>%
    mutate(gene = factor(gene, levels = gene_order), label = paste0(n_samples))
  
  p_cancer_bar <- ggplot(sample_count_df, aes(x = n_samples, y = gene)) +
    geom_col(width = 0.7, fill = cancer_palette[ct], colour = NA) +
    geom_text(aes(label = label), hjust = -0.1, size = 2.2, colour = "grey30") +
    scale_x_continuous(name = "N samples", expand = expansion(mult = c(0, 0.35))) +
    scale_y_discrete(drop = FALSE) +
    labs(y = NULL, title = paste0("N samples\n(total = ", n_total_samples, ")")) +
    theme_minimal(base_size = 9) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 7),
          axis.ticks.x = element_line(linewidth = 0.3),
          panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
          plot.margin = margin(4, 8, 4, 2))
  
  # Mutation heatmap (with timed border + multiplicity pattern)
  p_mut_heat <- build_mutation_heatmap(heatmap_df, heatmap_theme)
  
  # Metadata strips
  p_meta_class <- make_meta_strip(
    sample_meta, "class", "Class",
    scale_fill_manual(values = class_palette, name = "Class", na.value = "grey95",
                      drop = FALSE, guide = guide_legend(nrow = LEGEND_NROW, title.position = "top")),
    meta_theme)
  p_meta_ncna <- make_meta_strip(
    sample_meta, "n_cna", "N CNA",
    scale_fill_gradientn(colours = c("#fee8c8","#e34a33"), name = "N CNA", na.value = "grey95",
                         guide = guide_colourbar(title.position = "top",
                                                 barwidth = unit(1.5,"cm"), barheight = unit(0.3,"cm"))),
    meta_theme)
  p_meta_purity <- make_meta_strip(
    sample_meta, "purity", "Purity",
    scale_fill_gradientn(colours = c("#edf8e9","#238b45"), name = "Purity",
                         na.value = "grey95", limits = c(0,1),
                         guide = guide_colourbar(title.position = "top",
                                                 barwidth = unit(1.5,"cm"), barheight = unit(0.3,"cm"))),
    meta_theme)
  p_meta_ploidy <- make_meta_strip(
    sample_meta, "ploidy", "Ploidy",
    scale_fill_gradientn(colours = c("#f0f0f0","#636363"), name = "Ploidy", na.value = "grey95",
                         guide = guide_colourbar(title.position = "top",
                                                 barwidth = unit(1.5,"cm"), barheight = unit(0.3,"cm"))),
    meta_theme)
  p_meta_ncomp <- make_meta_strip(
    sample_meta, "ncomponents", "N comp.",
    scale_fill_gradientn(colours = c("#ffffd4","#993404"), name = "N comp.", na.value = "grey95",
                         guide = guide_colourbar(title.position = "top",
                                                 barwidth = unit(1.5,"cm"), barheight = unit(0.3,"cm"))),
    meta_theme)
  p_meta_bestk <- make_meta_strip(
    sample_meta, "best_K", "Best K",
    scale_fill_gradientn(colours = c("#f7fbff","#2171b5"), name = "Best K", na.value = "grey95",
                         guide = guide_colourbar(title.position = "top",
                                                 barwidth = unit(1.5,"cm"), barheight = unit(0.3,"cm"))),
    meta_theme)
  
  meta_stack <- wrap_plots(
    p_meta_class, p_meta_ncna, p_meta_purity,
    p_meta_ploidy, p_meta_ncomp, p_meta_bestk,
    ncol = 1
  ) + plot_layout(heights = rep(1, 6))
  
  layout_design <- "
  ABCD
  ##E#
"
  final_plot <- (p_kar_bar + p_timing + p_mut_heat + p_cancer_bar + meta_stack) +
    plot_layout(design = layout_design, widths = c(0.1, 0.1, 0.3, 0.1),
                heights = c(4, 1), guides = "collect") +
    plot_annotation(
      title    = paste0("Driver CNA annotation summary — ", ct),
      subtitle = paste0(
        n_genes, " genes (\u2265", MIN_GENE_PROP * 100, "% of ", n_samples, " samples)  |  ",
        "Heatmap: black border = timed segment | X = mult 1 | * = mult 2"
      ),
      theme = theme(plot.title    = element_text(size = 13, face = "bold"),
                    plot.subtitle = element_text(size = 8,  colour = "grey40"))
    ) & two_row_legend_theme
  
  save_plot(final_plot, file.path(out_dir_pseudotime, paste0("drivers_summary_", ct, ".png")),
            n_samples, n_genes)
}


# ══════════════════════════════════════════════════════════════════════════════
# CANCER TYPE LOOP — CLOCK RANK PLOTS
# Replaces timing violin with:
#   Panel A: clock_rank stacked bar per gene
#   Panel B: max_clock_rank of host sample per gene
# Saved to: out_dir/drivers+mutation+clockrank/
# ══════════════════════════════════════════════════════════════════════════════

for (ct in cancer_type_list) {
  
  message("\n── Processing (clock rank): ", ct, " ──")
  
  drivers_ct <- Drivers %>% filter(IntoGen_cancer_type == ct)
  if (nrow(drivers_ct) == 0) { message("  No data, skipping."); next }
  
  n_total_samples <- n_distinct(drivers_ct$sample)
  gene_props <- drivers_ct %>%
    group_by(gene) %>%
    summarise(prop = n_distinct(sample) / n_total_samples, .groups = "drop")
  genes_keep <- gene_props %>% filter(prop >= MIN_GENE_PROP) %>% pull(gene)
  if (length(genes_keep) == 0) { message("  No genes pass filter, skipping."); next }
  
  drivers_ct <- drivers_ct %>% filter(gene %in% genes_keep)
  
  gene_order <- drivers_ct %>%
    filter(!is.na(clock_mean)) %>%
    group_by(gene) %>%
    summarise(median_clock = median(clock_mean, na.rm = TRUE), .groups = "drop") %>%
    arrange(median_clock) %>% pull(gene)
  gene_order <- c(setdiff(unique(drivers_ct$gene), gene_order), gene_order)
  
  sample_order <- drivers_ct %>%
    distinct(sample, class, wgd_status) %>%
    arrange(class, wgd_status) %>% pull(sample) %>% unique()
  
  n_genes   <- length(gene_order)
  n_samples <- length(sample_order)
  
  heatmap_df <- build_heatmap_df(drivers_ct, gene_order, sample_order)
  
  sample_meta <- drivers_ct %>%
    distinct(sample, class, wgd_status, n_cna, purity,
             ncomponents, best_K, IntoGen_cancer_type, ploidy) %>%
    mutate(sample = factor(sample, levels = sample_order),
           class  = factor(class,  levels = names(class_palette)))
  
  heatmap_theme <- theme_minimal(base_size = 9) + two_row_legend_theme +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid  = element_blank(),
          plot.title  = element_text(size = 8, face = "bold", hjust = 0.5),
          plot.margin = margin(4, 2, 4, 2))
  
  meta_theme <- theme_minimal(base_size = 9) + two_row_legend_theme +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 7, hjust = 1),
          axis.ticks.y = element_blank(), panel.grid = element_blank(),
          plot.margin = margin(1, 2, 1, 2))
  
  # ── Panel A: clock_rank distribution per gene ─────────────────────────────
  clock_rank_bar_df <- drivers_ct %>%
    filter(!is.na(clock_rank)) %>%
    mutate(gene = factor(gene, levels = gene_order)) %>%
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
    labs(x = "% segments", y = NULL, title = "Within-sample\nclock rank") +
    theme_minimal(base_size = 9) + two_row_legend_theme +
    theme(axis.text.y = element_text(size = 8, hjust = 1, face = "italic"),
          axis.text.x = element_text(size = 7),
          axis.ticks.x = element_line(linewidth = 0.3),
          panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
          plot.margin = margin(4, 4, 4, 4))
  
  # ── Panel B: max_clock_rank of host sample per gene ───────────────────────
  # Uses max_clock_rank_per_sample computed from Segments BEFORE driver join
  gene_sample_max_rank <- drivers_ct %>%
    filter(!is.na(clock_rank)) %>%
    distinct(gene, sample) %>%
    left_join(max_clock_rank_per_sample, by = "sample") %>%
    filter(!is.na(max_clock_rank)) %>%
    mutate(gene = factor(gene, levels = gene_order))
  
  max_rank_bar_df <- gene_sample_max_rank %>%
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
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 7),
          axis.ticks.x = element_line(linewidth = 0.3),
          panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
          plot.margin = margin(4, 4, 4, 2))
  
  # Mutation heatmap (with timed border + multiplicity pattern)
  p_mut_heat <- build_mutation_heatmap(heatmap_df, heatmap_theme)
  
  # N samples bar
  sample_count_df <- drivers_ct %>%
    distinct(gene, sample) %>%
    count(gene, name = "n_samples") %>%
    mutate(gene = factor(gene, levels = gene_order), label = paste0(n_samples))
  
  p_cancer_bar <- ggplot(sample_count_df, aes(x = n_samples, y = gene)) +
    geom_col(width = 0.7, fill = cancer_palette[ct], colour = NA) +
    geom_text(aes(label = label), hjust = -0.1, size = 2.2, colour = "grey30") +
    scale_x_continuous(name = "N samples", expand = expansion(mult = c(0, 0.35))) +
    scale_y_discrete(drop = FALSE) +
    labs(y = NULL, title = paste0("N samples\n(total = ", n_total_samples, ")")) +
    theme_minimal(base_size = 9) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 7),
          axis.ticks.x = element_line(linewidth = 0.3),
          panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
          plot.margin = margin(4, 8, 4, 2))
  
  # Metadata strips
  p_meta_class <- make_meta_strip(
    sample_meta, "class", "Class",
    scale_fill_manual(values = class_palette, name = "Class", na.value = "grey95",
                      drop = FALSE, guide = guide_legend(nrow = LEGEND_NROW, title.position = "top")),
    meta_theme)
  p_meta_ncna <- make_meta_strip(
    sample_meta, "n_cna", "N CNA",
    scale_fill_gradientn(colours = c("#fee8c8","#e34a33"), name = "N CNA", na.value = "grey95",
                         guide = guide_colourbar(title.position = "top",
                                                 barwidth = unit(1.5,"cm"), barheight = unit(0.3,"cm"))),
    meta_theme)
  p_meta_purity <- make_meta_strip(
    sample_meta, "purity", "Purity",
    scale_fill_gradientn(colours = c("#edf8e9","#238b45"), name = "Purity",
                         na.value = "grey95", limits = c(0,1),
                         guide = guide_colourbar(title.position = "top",
                                                 barwidth = unit(1.5,"cm"), barheight = unit(0.3,"cm"))),
    meta_theme)
  p_meta_ploidy <- make_meta_strip(
    sample_meta, "ploidy", "Ploidy",
    scale_fill_gradientn(colours = c("#f0f0f0","#636363"), name = "Ploidy", na.value = "grey95",
                         guide = guide_colourbar(title.position = "top",
                                                 barwidth = unit(1.5,"cm"), barheight = unit(0.3,"cm"))),
    meta_theme)
  p_meta_ncomp <- make_meta_strip(
    sample_meta, "ncomponents", "N comp.",
    scale_fill_gradientn(colours = c("#ffffd4","#993404"), name = "N comp.", na.value = "grey95",
                         guide = guide_colourbar(title.position = "top",
                                                 barwidth = unit(1.5,"cm"), barheight = unit(0.3,"cm"))),
    meta_theme)
  p_meta_bestk <- make_meta_strip(
    sample_meta, "best_K", "Best K",
    scale_fill_gradientn(colours = c("#f7fbff","#2171b5"), name = "Best K", na.value = "grey95",
                         guide = guide_colourbar(title.position = "top",
                                                 barwidth = unit(1.5,"cm"), barheight = unit(0.3,"cm"))),
    meta_theme)
  
  meta_stack <- wrap_plots(
    p_meta_class, p_meta_ncna, p_meta_purity,
    p_meta_ploidy, p_meta_ncomp, p_meta_bestk,
    ncol = 1
  ) + plot_layout(heights = rep(1, 6))
  
  layout_design_rank <- "
  ABCD
  ##E#
"
  final_plot_rank <- (
    p_clock_rank_bar + p_max_rank_bar + p_mut_heat + p_cancer_bar + meta_stack
  ) +
    plot_layout(design = layout_design_rank,
                widths = c(0.13, 0.10, 0.30, 0.10),
                heights = c(4, 1), guides = "collect") +
    plot_annotation(
      title    = paste0("Driver CNA — clock rank summary — ", ct),
      subtitle = paste0(
        n_genes, " genes (\u2265", MIN_GENE_PROP * 100, "% of ", n_samples, " samples)  |  ",
        "Heatmap: black border = timed | X = mult 1 | * = mult 2  |  ",
        "Clock rank derived from Segments.rds"
      ),
      theme = theme(plot.title    = element_text(size = 13, face = "bold"),
                    plot.subtitle = element_text(size = 8,  colour = "grey40"))
    ) & two_row_legend_theme
  
  save_plot(final_plot_rank,
            file.path(out_dir_rank, paste0("drivers_clockrank_", ct, ".png")),
            n_samples, n_genes)
}

message("\n✔  All done.")