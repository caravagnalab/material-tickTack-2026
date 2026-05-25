mutation_palette <- c(
  "CI_M"       = "#d9473f",
  "M"          = "#ff9c57",
  "CNA_driver" = "#7b2d8b",
  "WT"         = "#eeeeee"
)

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

# Assumes these objects are already in the environment from enrichment_tests.R:
#   results_all, count_table, hm_signature_summary, hm_pancancer_sig
#   EVENT_LEVELS, CLASS_LEVELS, mutation_palette, class_palette, out_dir

out_dir_stats <- file.path(out_dir, "driver_enrichment")
dir.create(out_dir_stats, showWarnings = FALSE, recursive = TRUE)

base_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(face = "bold", size = 9),
    plot.title        = element_text(face = "bold", size = 11),
    plot.subtitle     = element_text(size = 8, colour = "grey40"),
    legend.title      = element_text(face = "bold", size = 8),
    legend.text       = element_text(size = 7)
  )

sig_label <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE       ~ ""
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 1 — HM enrichment heatmap
# log2(OR) for per_event_HM_vs_rest across cancer types × genes, faceted by event
# Only genes with at least one locally-significant result
# ══════════════════════════════════════════════════════════════════════════════

heat_df <- results_all %>%
  filter(
    IntoGen_cancer_type != "pan_cancer",
    test_type == "per_event_HM_vs_rest",
    !is.na(odds_ratio)
  ) %>%
  mutate(
    log2_OR   = log2(odds_ratio),
    log2_OR_c = pmax(pmin(log2_OR, 3), -3),   # cap at ±3 for colour scale
    sig       = sig_label(p_adj_local),
    event     = factor(event, levels = EVENT_LEVELS)
  )

sig_genes_heat <- heat_df %>%
  filter(p_adj_local < 0.05) %>%
  pull(gene) %>% unique()

if (length(sig_genes_heat) > 0) {

  # Order genes: most-significant-cancer-types first
  gene_ord_heat <- heat_df %>%
    filter(gene %in% sig_genes_heat, p_adj_local < 0.05) %>%
    count(gene, sort = TRUE) %>%
    pull(gene)

  heat_df_sig <- heat_df %>%
    filter(gene %in% sig_genes_heat) %>%
    mutate(gene = factor(gene, levels = rev(gene_ord_heat)))

  p1 <- ggplot(heat_df_sig,
               aes(x = IntoGen_cancer_type, y = gene, fill = log2_OR_c)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(aes(label = sig), size = 2.2, vjust = 0.75, colour = "grey15") +
    facet_wrap(~ event, nrow = 1) +
    scale_fill_gradientn(
      colours  = c("#2166ac", "#92c5de", "#f7f7f7", "#f4a582", "#d6604d"),
      limits   = c(-3, 3),
      na.value = "grey88",
      name     = "log₂(OR)\n[capped ±3]",
      guide    = guide_colourbar(
        barwidth  = unit(0.4, "cm"),
        barheight = unit(3.0, "cm"),
        title.position = "top"
      )
    ) +
    labs(
      title    = "HM enrichment heatmap — per event, per cancer type",
      subtitle = "log₂(OR) vs rest (Classic + WGD)  |  * pₐₐⱼ<0.05  ** <0.01  *** <0.001  (BH local)",
      x = NULL, y = NULL
    ) +
    base_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(face = "italic", size = 8),
      panel.spacing = unit(0.3, "cm")
    )

  n_ct   <- n_distinct(heat_df_sig$IntoGen_cancer_type)
  n_gene <- length(sig_genes_heat)
  ggsave(
    file.path(out_dir_stats, "01_hm_enrichment_heatmap.png"), p1,
    width  = max(10, n_ct * 0.45 * 4 + 5),   # 4 event panels
    height = max(4,  n_gene * 0.38 + 3),
    dpi = 180
  )
  message("Saved → 01_hm_enrichment_heatmap.png  (",
          n_gene, " genes, ", n_ct, " cancer types)")

} else {
  message("Plot 1 skipped: no locally-significant HM enrichments found.")
  p1 <- NULL
}

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 2 — Cross-cancer HM signature bubble plot
# From hm_signature_summary: median OR vs gene, bubble = n cancer types
# ══════════════════════════════════════════════════════════════════════════════

if (nrow(hm_signature_summary) > 0) {

  sig_sum_pl <- hm_signature_summary %>%
    mutate(
      event = factor(event, levels = EVENT_LEVELS),
      gene  = factor(gene,  levels = rev(unique(gene)))
    )

  p2 <- ggplot(sig_sum_pl,
               aes(x = median_OR, y = gene,
                   size = n_cancer_types_enriched,
                   colour = event)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               colour = "grey55", linewidth = 0.5) +
    geom_point(alpha = 0.82) +
    scale_x_log10(
      breaks = c(0.5, 1, 2, 5, 10, 20),
      labels = c("0.5", "1", "2", "5", "10", "20"),
      name   = "Median OR (log scale)"
    ) +
    scale_size_continuous(
      name   = "N cancer types\n(enriched in HM)",
      range  = c(2, 9),
      breaks = seq_len(max(sig_sum_pl$n_cancer_types_enriched))
    ) +
    scale_colour_manual(values = mutation_palette, name = "Event", drop = FALSE) +
    labs(
      title    = "Cross-cancer HM signature",
      subtitle = "Genes where event is enriched in HM (OR > 1, pₐₐⱼ < 0.05) in ≥1 cancer type",
      y = NULL
    ) +
    base_theme +
    theme(
      axis.text.y        = element_text(face = "italic"),
      axis.text.x        = element_text(angle = 0, hjust = 0.5),
      panel.grid.major.x = element_line(colour = "grey85"),
      panel.grid.major.y = element_line(colour = "grey92")
    )

  ggsave(
    file.path(out_dir_stats, "02_hm_signature_bubbleplot.png"), p2,
    width  = 7,
    height = max(4, nrow(sig_sum_pl) * 0.38 + 3),
    dpi = 180
  )
  message("Saved → 02_hm_signature_bubbleplot.png  (",
          nrow(sig_sum_pl), " gene×event combinations)")

} else {
  message("Plot 2 skipped: hm_signature_summary is empty.")
  p2 <- NULL
}

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 3 — Pan-cancer forest plot
# OR ± 95% CI for pan-cancer per_event_HM_vs_rest significant hits
# ══════════════════════════════════════════════════════════════════════════════

forest_df <- results_all %>%
  filter(
    IntoGen_cancer_type == "pan_cancer",
    test_type           == "per_event_HM_vs_rest",
    !is.na(odds_ratio),
    p_adj_local         < 0.05
  ) %>%
  mutate(
    event     = factor(event, levels = EVENT_LEVELS),
    direction = ifelse(odds_ratio > 1, "Enriched in HM", "Depleted in HM"),
    # clip CI to 0.05 so log scale doesn't break on near-zero lower bounds
    ci_low_c  = pmax(ci_low,  0.05),
    ci_high_c = pmin(ci_high, 200),
    gene      = reorder(gene, log2(odds_ratio))
  )

if (nrow(forest_df) > 0) {

  p3 <- ggplot(forest_df,
               aes(x = odds_ratio, y = gene,
                   colour = event, shape = direction)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               colour = "grey40", linewidth = 0.6) +
    geom_errorbarh(
      aes(xmin = ci_low_c, xmax = ci_high_c),
      height = 0.28, linewidth = 0.55, alpha = 0.65
    ) +
    geom_point(size = 2.4) +
    facet_wrap(~ event, scales = "free_y", nrow = 1) +
    scale_x_log10(
      breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10, 20),
      labels = c("0.1","0.25","0.5","1","2","5","10","20"),
      name   = "Odds ratio (log scale) — HM vs rest"
    ) +
    scale_colour_manual(values = mutation_palette, name = "Event", drop = FALSE) +
    scale_shape_manual(
      values = c("Enriched in HM" = 19, "Depleted in HM" = 1),
      name   = "Direction"
    ) +
    labs(
      title    = "Pan-cancer forest plot — HM vs rest",
      subtitle = "Significant hits (pₐₐⱼ < 0.05, BH local) | filled = enriched in HM, open = depleted",
      y = NULL
    ) +
    base_theme +
    theme(
      axis.text.y   = element_text(face = "italic", size = 8),
      axis.text.x   = element_text(angle = 0, hjust = 0.5),
      panel.spacing = unit(0.4, "cm")
    )

  n_ev_f <- n_distinct(forest_df$event)
  ggsave(
    file.path(out_dir_stats, "03_pancancer_forest.png"), p3,
    width  = max(8, n_ev_f * 3.2 + 3),
    height = max(4, n_distinct(forest_df$gene) * 0.42 + 3),
    dpi = 180
  )
  message("Saved → 03_pancancer_forest.png  (",
          n_distinct(forest_df$gene), " genes)")

} else {
  message("Plot 3 skipped: no pan-cancer significant hits.")
  p3 <- NULL
}

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 4 — Event proportion stacked bars for top HM-enriched genes
# Pan-cancer count_table for the top genes from hm_signature_summary
# (falls back to hm_pancancer_sig if summary is empty)
# ══════════════════════════════════════════════════════════════════════════════

top_genes_prop <- if (nrow(hm_signature_summary) > 0) {
  hm_signature_summary %>%
    arrange(desc(n_cancer_types_enriched), desc(median_OR)) %>%
    distinct(gene) %>%
    head(24) %>%
    pull(gene)
} else {
  hm_pancancer_sig %>%
    arrange(p_adj_local, desc(abs(log2(odds_ratio)))) %>%
    pull(gene) %>% unique() %>% head(24)
}

if (length(top_genes_prop) > 0) {

  prop_df <- count_table %>%
    filter(
      IntoGen_cancer_type == "pan_cancer",
      gene  %in% top_genes_prop
    ) %>%
    mutate(
      gene  = factor(gene,  levels = top_genes_prop),
      event = factor(event, levels = EVENT_LEVELS),
      class = factor(class, levels = CLASS_LEVELS)
    )

  p4 <- ggplot(prop_df, aes(x = class, y = prop, fill = event)) +
    geom_col(width = 0.78, position = "stack",
             colour = "white", linewidth = 0.2) +
    geom_text(
      data    = prop_df %>% filter(prop >= 0.07),
      aes(label = percent(prop, accuracy = 1)),
      position = position_stack(vjust = 0.5),
      size = 2.0, colour = "white", fontface = "bold"
    ) +
    facet_wrap(~ gene, ncol = 6) +
    scale_fill_manual(values = mutation_palette, name = "Event", drop = FALSE) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.03))
    ) +
    scale_x_discrete(drop = FALSE) +
    labs(
      title    = "Event proportions per class — top HM-enriched genes (pan-cancer)",
      subtitle = "Proportion of WGD / HM / Classic samples carrying each event",
      x = NULL, y = "Proportion of samples"
    ) +
    base_theme +
    theme(
      axis.text.x  = element_text(angle = 40, hjust = 1, size = 7),
      strip.text   = element_text(face = "bold.italic", size = 7.5),
      panel.spacing = unit(0.3, "cm")
    )

  n_rows_facet <- ceiling(length(top_genes_prop) / 6)
  ggsave(
    file.path(out_dir_stats, "04_event_proportions_top_genes.png"), p4,
    width  = 14,
    height = max(5, n_rows_facet * 3.0 + 2),
    dpi = 180
  )
  message("Saved → 04_event_proportions_top_genes.png  (",
          length(top_genes_prop), " genes)")

} else {
  message("Plot 4 skipped: no top genes to display.")
  p4 <- NULL
}

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 5 — Pairwise comparison summary heatmap
# -log10(p_adj_local) for HM_vs_Classic and HM_vs_WGD, per event × gene
# Helps see whether HM differs more from Classic or from WGD
# ══════════════════════════════════════════════════════════════════════════════

pairwise_df <- results_all %>%
  filter(
    IntoGen_cancer_type == "pan_cancer",
    test_type           == "per_event_pairwise",
    comparison          %in% c("HM_vs_Classic", "HM_vs_WGD"),
    !is.na(p_adj_local)
  ) %>%
  mutate(
    nlp       = pmin(-log10(p_adj_local), 4),   # cap at 10^-4
    log2_OR_c = pmax(pmin(log2(odds_ratio), 3), -3),
    sig       = sig_label(p_adj_local),
    event     = factor(event,      levels = EVENT_LEVELS),
    comparison = factor(comparison, levels = c("HM_vs_Classic", "HM_vs_WGD"))
  )

# Keep only genes significant in at least one pairwise comparison
sig_genes_pw <- pairwise_df %>%
  filter(p_adj_local < 0.05) %>%
  pull(gene) %>% unique()

if (length(sig_genes_pw) > 0) {

  gene_ord_pw <- pairwise_df %>%
    filter(gene %in% sig_genes_pw) %>%
    group_by(gene) %>%
    summarise(mean_nlp = mean(nlp, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_nlp) %>%
    pull(gene)

  pairwise_sig <- pairwise_df %>%
    filter(gene %in% sig_genes_pw) %>%
    mutate(gene = factor(gene, levels = gene_ord_pw))

  p5 <- ggplot(pairwise_sig,
               aes(x = comparison, y = gene, fill = log2_OR_c)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sig), size = 2.2, vjust = 0.75, colour = "grey15") +
    facet_wrap(~ event, nrow = 1) +
    scale_fill_gradientn(
      colours  = c("#2166ac", "#92c5de", "#f7f7f7", "#f4a582", "#d6604d"),
      limits   = c(-3, 3),
      na.value = "grey88",
      name     = "log₂(OR)\n[capped ±3]",
      guide    = guide_colourbar(
        barwidth  = unit(0.4, "cm"),
        barheight = unit(3.0, "cm"),
        title.position = "top"
      )
    ) +
    labs(
      title    = "Pan-cancer pairwise: HM vs Classic and HM vs WGD",
      subtitle = "Per-event log₂(OR) | OR > 0 = enriched in HM  |  * pₐₐⱼ<0.05  ** <0.01  *** <0.001",
      x = NULL, y = NULL
    ) +
    base_theme +
    theme(
      axis.text.x   = element_text(angle = 30, hjust = 1, size = 8),
      axis.text.y   = element_text(face = "italic", size = 8),
      panel.spacing = unit(0.3, "cm")
    )

  ggsave(
    file.path(out_dir_stats, "05_pancancer_pairwise_heatmap.png"), p5,
    width  = max(10, 4 * 2.5 + 4),
    height = max(4, length(sig_genes_pw) * 0.38 + 3),
    dpi = 180
  )
  message("Saved → 05_pancancer_pairwise_heatmap.png  (",
          length(sig_genes_pw), " genes)")

} else {
  message("Plot 5 skipped: no significant pairwise pan-cancer hits.")
  p5 <- NULL
}

message("\n✔  All enrichment plots saved to: ", out_dir_stats)
