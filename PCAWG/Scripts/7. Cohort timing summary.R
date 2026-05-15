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

# ══════════════════════════════════════════════════════════════════════════════
# LOAD DATA
# ══════════════════════════════════════════════════════════════════════════════

Drivers <- readRDS("~/Documents/material-tickTack-2026/PCAWG/Data/Drivers.rds")

out_dir <- "~/Documents/material-tickTack-2026/PCAWG/Plot/"

# ══════════════════════════════════════════════════════════════════════════════
# PARAMETERS
# ══════════════════════════════════════════════════════════════════════════════

# Keep genes mutated in at least this proportion of samples within cancer type
MIN_GENE_PROP <- 0.30   # 5% — adjust as needed

# ══════════════════════════════════════════════════════════════════════════════
# LEGEND LAYOUT HELPER
# Splits all legends into two rows of adjacent guides using a custom grob
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

# Each individual guide wraps its keys into at most this many rows
LEGEND_NROW <- 2

# ══════════════════════════════════════════════════════════════════════════════
# PER-CANCER-TYPE LOOP
# ══════════════════════════════════════════════════════════════════════════════

cancer_type_list <- unique(na.omit(Drivers$IntoGen_cancer_type))

ct = "PRAD"

for (ct in cancer_type_list) {
  
  message("\n── Processing: ", ct, " ──")
  
  # ── Filter to cancer type ─────────────────────────────────────────────────
  drivers_ct <- Drivers %>%
    filter(IntoGen_cancer_type == ct)
  
  if (nrow(drivers_ct) == 0) {
    message("  No data, skipping.")
    next
  }
  
  # total samples for this cancer type
  n_total_samples <- n_distinct(drivers_ct$sample)
  message("  Total samples: ", n_total_samples)
  
  # ── Gene filter: keep genes present in >= MIN_GENE_PROP of samples ────────
  gene_props <- drivers_ct %>%
    group_by(gene) %>%
    summarise(prop = n_distinct(sample) / n_total_samples,
              .groups = "drop")
  
  genes_keep <- gene_props %>%
    filter(prop >= MIN_GENE_PROP) %>%
    pull(gene)
  
  if (length(genes_keep) == 0) {
    message("  No genes pass proportion filter, skipping.")
    next
  }
  
  drivers_ct <- drivers_ct %>% filter(gene %in% genes_keep)
  message("  Genes retained: ", length(genes_keep))
  
  # ── Gene order: by median clock_mean (early→late, bottom→top) ─────────────
  gene_order <- drivers_ct %>%
    filter(!is.na(clock_mean)) %>%
    group_by(gene) %>%
    summarise(median_clock = median(clock_mean, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(median_clock) %>%
    pull(gene)
  
  genes_notimed <- setdiff(unique(drivers_ct$gene), gene_order)
  gene_order    <- c(genes_notimed, gene_order)
  
  # ── Sample order: by class then wgd_status ────────────────────────────────
  sample_order <- drivers_ct %>%
    distinct(sample, class, wgd_status) %>%
    arrange(class, wgd_status) %>%
    pull(sample) %>%
    unique()
  
  n_genes   <- length(gene_order)
  n_samples <- length(sample_order)
  
  # ── heatmap_df: one row per gene × sample ─────────────────────────────────
  heatmap_df <- drivers_ct %>%
    mutate(
      gene   = factor(gene,   levels = gene_order),
      sample = factor(sample, levels = sample_order),
      karyotype_class = case_when(
        karyotype == "2:1" ~ "2:1",
        karyotype == "2:0" ~ "2:0",
        karyotype == "2:2" ~ "2:2",
        karyotype == "1:1" ~ "1:1",
        karyotype == "1:0" ~ "1:0",
        TRUE               ~ "other"
      ),
      karyotype_class = factor(karyotype_class, levels = names(karyotype_palette)),
      mutation_status = factor(mutation_status, levels = names(mutation_palette)),
      class           = factor(class,           levels = names(class_palette)),
      arm_level_class = factor(
        ifelse(is.na(arm_level_class), "other", arm_level_class),
        levels = names(arm_palette)
      )
    ) %>%
    group_by(gene, sample) %>%
    arrange(match(mutation_status, c("CNA_driver","CI_M","M","WT"))) %>%
    slice(1) %>%
    ungroup()
  
  # ── sample_meta ───────────────────────────────────────────────────────────
  sample_meta <- drivers_ct %>%
    distinct(sample, class, wgd_status, n_cna, purity,
             ncomponents, best_K, IntoGen_cancer_type, ploidy) %>%
    mutate(
      sample = factor(sample, levels = sample_order),
      class  = factor(class,  levels = names(class_palette))
    )
  
  # ── PANEL: karyotype stacked bar ──────────────────────────────────────────
  kar_bar_df <- heatmap_df %>%
    filter(!is.na(karyotype_class)) %>%
    count(gene, karyotype_class) %>%
    group_by(gene) %>%
    mutate(pct = n / sum(n)) %>%
    ungroup() %>%
    mutate(gene = factor(gene, levels = gene_order))
  
  p_kar_bar <- ggplot(kar_bar_df,
                      aes(x = pct, y = gene, fill = karyotype_class)) +
    geom_col(width = 0.7, position = "stack") +
    scale_fill_manual(values = karyotype_palette, name = "Karyotype",
                      drop = FALSE,
                      guide = guide_legend(nrow = LEGEND_NROW,
                                           title.position = "top")) +
    scale_x_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "% samples", y = NULL, title = "Karyotype\ndistribution") +
    theme_minimal(base_size = 9) +
    two_row_legend_theme +
    theme(
      axis.text.y        = element_text(size = 8, hjust = 1, face = "italic"),
      axis.text.x        = element_text(size = 7),
      axis.ticks.x       = element_line(linewidth = 0.3),
      panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(size = 8, face = "bold", hjust = 0.5),
      plot.margin        = margin(4, 4, 4, 4)
    )
  
  # ── PANEL: timing violin ──────────────────────────────────────────────────
  timing_df <- drivers_ct %>%
    filter(!is.na(clock_mean)) %>%
    mutate(gene = factor(gene, levels = gene_order))
  
  timing_summary <- timing_df %>%
    group_by(gene) %>%
    summarise(med = median(clock_mean),
              lo  = median(clock_low,  na.rm = TRUE),
              hi  = median(clock_high, na.rm = TRUE),
              .groups = "drop")
  
  p_timing <- ggplot(timing_df, aes(y = gene, x = clock_mean)) +
    geom_violin(fill = "#b2d8e8", colour = NA, alpha = 0.6,
                scale = "width", trim = TRUE) +
    geom_point(aes(colour = clock_mean), size = 1.2, alpha = 0.7,
               position = position_jitter(height = 0.15, seed = 42)) +
    geom_errorbarh(data = timing_summary,
                   aes(y = gene, xmin = lo, xmax = hi, x = NULL),
                   height = 0.25, colour = "grey30", linewidth = 0.5) +
    geom_point(data = timing_summary, aes(y = gene, x = med),
               shape = 21, fill = "white", colour = "grey20",
               size = 2.5, stroke = 0.8) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                       labels = c("0", "0.5", "1")) +
    scale_colour_gradientn(colours = c("#ffffcc","#41b6c4","#253494"),
                           limits = c(0, 1), guide = "none") +
    scale_y_discrete(drop = FALSE) +
    labs(x = "Clock mean", y = NULL, title = "Timing\ndistribution") +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      axis.text.x        = element_text(size = 7),
      axis.ticks.x       = element_line(linewidth = 0.3),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(size = 8, face = "bold", hjust = 0.5),
      plot.margin        = margin(4, 4, 4, 2)
    )
  
  # ── PANEL: cancer type stacked bar ────────────────────────────────────────
  # one colour per cancer type present across the full dataset
  all_cancer_types <- sort(unique(na.omit(Drivers$IntoGen_cancer_type)))
  cancer_palette   <- setNames(
    colorRampPalette(c(
      "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99",
      "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a",
      "#ffff99","#b15928","#8dd3c7","#ffffb3","#bebada"
    ))(length(all_cancer_types)),
    all_cancer_types
  )
  
  cancer_bar_df <- drivers_ct %>%
    filter(!is.na(IntoGen_cancer_type)) %>%
    mutate(gene = factor(gene, levels = gene_order),
           IntoGen_cancer_type = factor(IntoGen_cancer_type,
                                        levels = all_cancer_types)) %>%
    distinct(gene, sample, IntoGen_cancer_type) %>%
    count(gene, IntoGen_cancer_type) %>%
    group_by(gene) %>%
    mutate(pct = n / sum(n)) %>%
    ungroup()
  
  p_cancer_bar_all_ttype <- ggplot(cancer_bar_df,
                         aes(x = pct, y = gene, fill = IntoGen_cancer_type)) +
    geom_col(width = 0.7, position = "stack") +
    scale_fill_manual(values = cancer_palette, name = "Cancer type",
                      drop = TRUE,
                      guide = guide_legend(nrow = LEGEND_NROW,
                                           title.position = "top")) +
    scale_x_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "% samples", y = NULL, title = "Cancer type\ndistribution") +
    theme_minimal(base_size = 9) +
    two_row_legend_theme +
    theme(
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      axis.text.x        = element_text(size = 7),
      axis.ticks.x       = element_line(linewidth = 0.3),
      panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(size = 8, face = "bold", hjust = 0.5),
      plot.margin        = margin(4, 4, 4, 2)
    )
  
  
  # ── PANEL: sample count bar (absolute N, horizontal) ──────────────────────
  sample_count_df <- drivers_ct %>%
    distinct(gene, sample) %>%
    count(gene, name = "n_samples") %>%
    mutate(
      gene  = factor(gene, levels = gene_order),
      # proportion label for annotation
      # prop  = n_samples / n_total_samples,
      label = paste0(n_samples)
                     # ,"\n(", round(prop * 100, 1), "%)")
    )
  
  p_cancer_bar <- ggplot(sample_count_df,
                         aes(x = n_samples, y = gene)) +
    geom_col(width = 0.7, fill = cancer_palette[ct], colour = NA) +
    geom_text(aes(label = label),
              hjust  = -0.1,
              size   = 2.2,
              colour = "grey30",
              lineheight = 0.85) +
    scale_x_continuous(
      name   = "N samples",
      expand = expansion(mult = c(0, 0.35))   # room for text labels
    ) +
    scale_y_discrete(drop = FALSE) +
    labs(y = NULL,
         title = paste0("N samples\n(total = ", n_total_samples, ")")) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      axis.text.x        = element_text(size = 7),
      axis.ticks.x       = element_line(linewidth = 0.3),
      panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(size = 8, face = "bold", hjust = 0.5),
      plot.margin        = margin(4, 8, 4, 2)   # extra right margin for labels
    )
  
  # ── HEATMAP PANELS ────────────────────────────────────────────────────────
  heatmap_theme <- theme_minimal(base_size = 9) +
    two_row_legend_theme +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid   = element_blank(),
      plot.title   = element_text(size = 8, face = "bold", hjust = 0.5),
      plot.margin  = margin(4, 2, 4, 2)
    )
  
  p_kar_heat <- ggplot(heatmap_df,
                       aes(x = sample, y = gene, fill = karyotype_class)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = karyotype_palette, name = "Karyotype",
                      na.value = "grey95", drop = FALSE,
                      guide = guide_legend(nrow = LEGEND_NROW,
                                           title.position = "top")) +
    scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL, title = "Karyotype") + heatmap_theme
  
  p_mut_heat <- ggplot(heatmap_df,
                       aes(x = sample, y = gene, fill = mutation_status)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = mutation_palette, name = "Mutation",
                      na.value = "grey95", drop = FALSE,
                      guide = guide_legend(nrow = LEGEND_NROW,
                                           title.position = "top")) +
    scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL, title = "Mutation status") + heatmap_theme
  
  p_clock_heat <- ggplot(heatmap_df,
                         aes(x = sample, y = gene, fill = clock_mean)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_gradientn(colours  = c("#ffffcc","#41b6c4","#253494"),
                         name     = "Timing",
                         na.value = "grey95", limits = c(0, 1),
                         guide    = guide_colourbar(title.position = "top",
                                                    barwidth = unit(2, "cm"),
                                                    barheight = unit(0.3, "cm"))) +
    scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL, title = "Timing") + heatmap_theme
  
  p_arm_heat <- ggplot(heatmap_df,
                       aes(x = sample, y = gene, fill = arm_level_class)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = arm_palette, name = "Arm level",
                      na.value = "grey95", drop = FALSE,
                      guide = guide_legend(nrow = LEGEND_NROW,
                                           title.position = "top")) +
    scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL, title = "Arm level") + heatmap_theme
  
  # ── METADATA STRIPS ───────────────────────────────────────────────────────
  meta_theme <- theme_minimal(base_size = 9) +
    two_row_legend_theme +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y  = element_text(size = 7, hjust = 1),
      axis.ticks.y = element_blank(),
      panel.grid   = element_blank(),
      plot.margin  = margin(1, 2, 1, 2)
    )
  
  make_meta_strip <- function(data, fill_var, y_label, scale_layer) {
    ggplot(data, aes(x = sample, y = y_label, fill = .data[[fill_var]])) +
      geom_tile(colour = "white", linewidth = 0.2) +
      scale_layer +
      scale_x_discrete(drop = FALSE) +
      labs(x = NULL, y = NULL) +
      meta_theme
  }
  
  p_meta_class <- make_meta_strip(
    sample_meta, "class", "Class",
    scale_fill_manual(values = class_palette, name = "Class",
                      na.value = "grey95", drop = FALSE,
                      guide = guide_legend(nrow = LEGEND_NROW,
                                           title.position = "top"))
  )
  p_meta_ncna <- make_meta_strip(
    sample_meta, "n_cna", "N CNA",
    scale_fill_gradientn(colours = c("#fee8c8","#e34a33"), name = "N CNA",
                         na.value = "grey95",
                         guide = guide_colourbar(title.position = "top",
                                                 barwidth  = unit(1.5, "cm"),
                                                 barheight = unit(0.3, "cm")))
  )
  p_meta_purity <- make_meta_strip(
    sample_meta, "purity", "Purity",
    scale_fill_gradientn(colours = c("#edf8e9","#238b45"), name = "Purity",
                         na.value = "grey95", limits = c(0, 1),
                         guide = guide_colourbar(title.position = "top",
                                                 barwidth  = unit(1.5, "cm"),
                                                 barheight = unit(0.3, "cm")))
  )
  p_meta_ploidy <- make_meta_strip(
    sample_meta, "ploidy", "Ploidy",
    scale_fill_gradientn(colours = c("#f0f0f0","#636363"), name = "Ploidy",
                         na.value = "grey95",
                         guide = guide_colourbar(title.position = "top",
                                                 barwidth  = unit(1.5, "cm"),
                                                 barheight = unit(0.3, "cm")))
  )
  p_meta_ncomp <- make_meta_strip(
    sample_meta, "ncomponents", "N comp.",
    scale_fill_gradientn(colours = c("#ffffd4","#993404"), name = "N comp.",
                         na.value = "grey95",
                         guide = guide_colourbar(title.position = "top",
                                                 barwidth  = unit(1.5, "cm"),
                                                 barheight = unit(0.3, "cm")))
  )
  p_meta_bestk <- make_meta_strip(
    sample_meta, "best_K", "Best K",
    scale_fill_gradientn(colours = c("#f7fbff","#2171b5"), name = "Best K",
                         na.value = "grey95",
                         guide = guide_colourbar(title.position = "top",
                                                 barwidth  = unit(1.5, "cm"),
                                                 barheight = unit(0.3, "cm")))
  )
  
  p_meta_cancer <- make_meta_strip(
    sample_meta %>%
      mutate(IntoGen_cancer_type = factor(IntoGen_cancer_type,
                                          levels = all_cancer_types)),
    "IntoGen_cancer_type", "Cancer type",
    scale_fill_manual(values = cancer_palette, name = "Cancer type",
                      na.value = "grey95", drop = TRUE,
                      guide = guide_legend(nrow = LEGEND_NROW,
                                           title.position = "top"))
  )
  
  p_meta_arm_level_events <- make_meta_strip(
    sample_meta %>%
      mutate(arm_level_class = factor(arm_level_class,
                                          levels = arm_level_class)),
    "IntoGen_cancer_type", "Cancer type",
    scale_fill_manual(values = cancer_palette, name = "Cancer type",
                      na.value = "grey95", drop = TRUE,
                      guide = guide_legend(nrow = LEGEND_NROW,
                                           title.position = "top"))
  )
  
  
  # ── ASSEMBLE ──────────────────────────────────────────────────────────────
  meta_stack <- wrap_plots(
    p_meta_class, p_meta_ncna, p_meta_purity,
    p_meta_ploidy, p_meta_ncomp, p_meta_bestk, 
    # p_meta_cancer,
    ncol = 1
  ) + plot_layout(heights = rep(1, 7))
  
  heatmap_stack <- wrap_plots(
    # p_kar_heat, 
    p_mut_heat,
    # p_clock_heat, 
    # p_arm_heat,
    ncol = 1
  ) + plot_layout(heights = rep(1, 4))
  
  centre_block <- wrap_plots(
    meta_stack,
    heatmap_stack,
    ncol    = 1,
    heights = c(7, max(4, 4 * n_genes / 10))
  )
  
  # 1. Define the layout grid
  # Row 1: A, B, C, D (Your 4 main plots)
  # Row 2: Empty, Empty, E (Metadata), Empty
  layout_design <- "
  ABCD
  ##E#
"
  
  final_plot <- (p_kar_bar + p_timing + p_mut_heat + 
                   p_cancer_bar +
                   meta_stack) +
    plot_layout(
      design = layout_design,
      widths  = c(0.1, 0.1, 0.3, 0.1), # Relative widths of columns A, B, C, D
      # widths  = c(0.1, 0.1, 0.3), 
      heights = c(4, 1),               # Height ratio between main plots and meta_stack
      guides  = "collect"
    ) +
    plot_annotation(
      title    = paste0("Driver CNA annotation summary — ", ct),
      subtitle = paste0(
        n_genes, " genes (\u2265", MIN_GENE_PROP * 100, "% of ",
        n_samples, " samples)  |  ",
        "Left: karyotype distribution  |  Centre: timing + heatmaps  |  ",
        "Right: cancer type distribution"
      ),
      theme = theme(
        plot.title    = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 8,  colour = "grey40")
      )
    ) &
    two_row_legend_theme
  
  final_plot
  
  # ── SAVE ──────────────────────────────────────────────────────────────────
  fig_w <- max(10, n_samples * 0.28 + 12)
  fig_h <- max(7,  n_genes   * 0.40 + 8)
  
  base_font <- max(9, min(16, round(sqrt(fig_w * fig_h) * 0.55)))
  
  text_override <- theme(
    text              = element_text(size = base_font),
    axis.text.x       = element_blank(),   # keep sample IDs hidden after override
    axis.ticks.x      = element_blank(),
    axis.text.x.top   = element_text(size = base_font * 0.85),
    axis.text.y       = element_text(size = base_font * 0.85),
    axis.text.y.left  = element_text(size = base_font * 0.85, hjust = 1,
                                     face = "italic"),
    axis.title        = element_text(size = base_font * 0.90),
    axis.title.x      = element_text(size = base_font * 0.90),
    axis.title.y      = element_text(size = base_font * 0.90),
    plot.title        = element_text(size = base_font * 1.30, face = "bold",
                                     hjust = 0.5),
    plot.subtitle     = element_text(size = base_font * 0.80, colour = "grey40"),
    legend.text       = element_text(size = base_font * 0.75),
    legend.title      = element_text(size = base_font * 0.80, face = "bold"),
    strip.text        = element_text(size = base_font * 0.85)
  )
  
  out_file <- file.path(out_dir, paste0("drivers_summary_", ct, ".png"))
  ggsave(
    final_plot & text_override,
    filename  = out_file,
    width     = fig_w,
    height    = fig_h,
    dpi       = 180,
    limitsize = FALSE
  )
  message("  Saved \u2192 ", out_file)
}

