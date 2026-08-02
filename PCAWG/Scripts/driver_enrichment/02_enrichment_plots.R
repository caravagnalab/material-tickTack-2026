mutation_palette <- c(
  "CI_M"       = "#d9473f",
  "M"          = "#ff9c57",
  "CNA_driver" = "#7b2d8b",
  "WT"         = "#eeeeee"
)

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

out_dir       <- "../../Plot/"
out_dir_stats <- paste0(out_dir, "driver_enrichment")

EVENT_LEVELS      <- c("CNA_driver", "CI_M", "M", "WT")
CLASS_LEVELS      <- c("WGD", "HM", "Classic")
COMPARISONS_PAIR  <- c("HM_vs_Classic", "HM_vs_WGD", "WGD_vs_Classic")
COMPARISONS_NO_WT <- c("HM_vs_WGD", "WGD_vs_Classic")

# Significance threshold for gene selection and annotation
SIG_THRESHOLD <- 0.1
# Maximum genes to show per plot
MAX_GENES <- 24

# ══════════════════════════════════════════════════════════════════════════════
# LOAD OBJECTS
# ══════════════════════════════════════════════════════════════════════════════

message("Loading objects ...")

results_all      <- readRDS(file.path(out_dir_stats, "enrichment_results.rds"))
count_table      <- readRDS(file.path(out_dir_stats, "event_counts.rds"))
hm_sig_summary   <- readRDS(file.path(out_dir_stats, "hm_signature_summary.rds"))
hm_pancancer_sig <- readRDS(file.path(out_dir_stats, "hm_pancancer_significant.rds"))
Drivers          <- readRDS("../../Data/Drivers.rds")

# Add zero_count column if missing (run enrichment script first to get it)
if (!"zero_count" %in% names(results_all)) {
  results_all <- results_all %>% mutate(zero_count = FALSE)
  message("Warning: zero_count missing — re-run enrichment script.")
}

message("Objects loaded.")

# ══════════════════════════════════════════════════════════════════════════════
# SHARED THEME
# ══════════════════════════════════════════════════════════════════════════════

dir.create(out_dir_stats, showWarnings = FALSE, recursive = TRUE)

base_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold", size = 9),
    plot.title       = element_text(face = "bold", size = 11),
    plot.subtitle    = element_text(size = 8, colour = "grey40"),
    legend.title     = element_text(face = "bold", size = 8),
    legend.text      = element_text(size = 7)
  )

# ══════════════════════════════════════════════════════════════════════════════
# ANNOTATION LOGIC — applied consistently in every plot
#
# Each (gene, event, comparison) cell can be one of three cases:
#
#  CASE 1 — normal:
#    both classes have samples, neither has zero event carriers
#    -> valid OR and CI, sig stars as annotation
#
#  CASE 2 — zero_count with valid OR  (zero_count == TRUE, !is.na(odds_ratio)):
#    one class has zero event carriers but still has samples
#    -> OR is 0 or Inf (computable but unreliable), CI is [0, Inf]
#    -> orange CI bar, "0" prefix on sig label
#
#  CASE 3 — NA OR  (is.na(odds_ratio)):
#    one class has zero samples for this gene entirely
#    -> OR cannot be computed
#    -> filled square at x=1, "NA" label
#
# ══════════════════════════════════════════════════════════════════════════════

# Significance stars (+ for weak, * ** *** for stronger)
sig_stars <- function(p) {
  case_when(
    is.na(p)          ~ "",
    p < 0.001         ~ "***",
    p < 0.01          ~ "**",
    p < 0.05          ~ "*",
    p < SIG_THRESHOLD ~ "+",
    TRUE              ~ ""
  )
}

# Text label shown inside/beside each data point or tile:
#   NA OR     -> "NA"
#   zero OR   -> "0 ***" (prefix + stars)
#   normal    -> "***" (stars only)
cell_label <- function(odds_ratio, zero_count, p_adj_local) {
  case_when(
    is.na(odds_ratio) ~ "NA",
    zero_count        ~ paste0("0 ", sig_stars(p_adj_local)),
    TRUE              ~ sig_stars(p_adj_local)
  )
}

# Drop WT for WGD-involving comparisons (trivially expected, not informative)
drop_wt_wgd <- function(df) {
  df %>% filter(!(comparison %in% COMPARISONS_NO_WT & event == "WT"))
}

# ══════════════════════════════════════════════════════════════════════════════
# GENE SELECTION
#
# For a given slice of results_all (one cancer type or pan-cancer):
# 1. Find genes significant (p_adj_local < SIG_THRESHOLD) in any
#    (comparison x event) cell with a valid OR
# 2. Find genes with at least one zero_count cell
# 3. Rank: significant first, then zero-count, then by median OR
# 4. Return top MAX_GENES names
# ══════════════════════════════════════════════════════════════════════════════

select_top_genes <- function(results_slice, max_genes = MAX_GENES) {
  
  # restrict to the tests and comparisons we care about
  df <- results_slice %>%
    filter(
      test_type  == "per_event_pairwise",
      comparison %in% COMPARISONS_PAIR
    ) %>%
    drop_wt_wgd()
  
  if (nrow(df) == 0) return(character(0))
  
  df %>%
    group_by(gene) %>%
    summarise(
      # TRUE if any cell is significant with a valid OR
      any_sig   = any(!is.na(odds_ratio) &
                        !is.na(p_adj_local) &
                        p_adj_local < SIG_THRESHOLD),
      # TRUE if any cell has a zero-count issue
      any_zero  = any(zero_count),
      # median OR across all cells for this gene (ignoring NAs)
      median_OR = median(odds_ratio, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    arrange(desc(any_sig), desc(any_zero), desc(median_OR)) %>%
    head(max_genes) %>%
    pull(gene)
}

# ══════════════════════════════════════════════════════════════════════════════
# PREPARE BASE DATA — pan-cancer gene selection done once here
# ══════════════════════════════════════════════════════════════════════════════

top_genes_pancancer <- select_top_genes(
  results_all %>% filter(IntoGen_cancer_type == "pan_cancer")
)

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 1 — Overview heatmap: all cancer types x all genes
#
# Each tile = one (gene, cancer_type, comparison, event) combination
# Colour    = log2(OR), capped at +/-3
# Border    = orange if zero_count (either case 2 or 3)
# Label     = "NA" / "0 ***" / "***" depending on case
#
# Only genes with at least one significant or zero-count result are shown
# to keep the plot readable.
# ══════════════════════════════════════════════════════════════════════════════

heat_raw <- results_all %>%
  filter(
    IntoGen_cancer_type != "pan_cancer",
    test_type  == "per_event_pairwise",
    comparison %in% COMPARISONS_PAIR
  ) %>%
  drop_wt_wgd() %>%
  mutate(
    log2_OR   = ifelse(is.na(odds_ratio), NA_real_,
                       pmax(pmin(log2(odds_ratio), 3), -3)),
    label     = cell_label(odds_ratio, zero_count, p_adj_local),
    event      = factor(event,      levels = EVENT_LEVELS),
    comparison = factor(comparison, levels = COMPARISONS_PAIR)
  )

# keep only genes with >=1 significant or zero-count result anywhere
genes_heat <- heat_raw %>%
  group_by(gene) %>%
  summarise(
    n_sig  = sum(!is.na(p_adj_local) & p_adj_local < SIG_THRESHOLD,
                 na.rm = TRUE),
    n_zero = sum(zero_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_sig > 0 | n_zero > 0) %>%
  arrange(desc(n_sig), desc(n_zero)) %>%
  pull(gene)

# fallback if nothing passes the filter
if (length(genes_heat) == 0) {
  genes_heat <- heat_raw %>%
    group_by(gene) %>%
    summarise(median_OR = median(odds_ratio, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(median_OR)) %>%
    head(MAX_GENES) %>%
    pull(gene)
}

heat_df <- heat_raw %>%
  filter(gene %in% genes_heat) %>%
  mutate(gene = factor(gene, levels = rev(genes_heat)))

if (nrow(heat_df) > 0) {
  
  n_ct        <- n_distinct(heat_df$IntoGen_cancer_type)
  n_gene      <- length(genes_heat)
  plot_width  <- min(120, max(10, n_ct * 0.45 * 4 + 5))
  plot_height <- min(120, max(8,  n_gene * 0.38 * 3 + 4))
  
  p1 <- ggplot(heat_df,
               aes(x = IntoGen_cancer_type, y = gene, fill = log2_OR)) +
    # normal tiles — white border
    geom_tile(data = filter(heat_df, !zero_count),
              colour = "white", linewidth = 0.25) +
    # zero-count tiles — orange border
    geom_tile(data = filter(heat_df, zero_count),
              colour = "#e08214", linewidth = 0.6) +
    geom_text(aes(label = label), size = 2.0, vjust = 0.75, colour = "grey15") +
    facet_grid(comparison ~ event) +
    scale_fill_gradientn(
      colours  = c("#2166ac", "#92c5de", "#f7f7f7", "#f4a582", "#d6604d"),
      limits   = c(-3, 3),
      na.value = "grey88",
      name     = "log2(OR)\n[+/-3]"
    ) +
    labs(
      title    = "Enrichment heatmap — all cancer types",
      subtitle = paste0(
        "Rows: HM vs Classic / HM vs WGD / WGD vs Classic  |  ",
        "Colour = log2(OR)  |  orange border = zero count in one class  |  ",
        "label: NA = class absent, 0 = event absent in one class, ",
        "+ p<", SIG_THRESHOLD, "  * p<0.05  ** p<0.01  *** p<0.001"
      ),
      x = NULL, y = NULL
    ) +
    base_theme +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y  = element_text(face = "italic", size = 8),
      panel.spacing = unit(0.3, "cm"),
      strip.text.y = element_text(angle = 0, hjust = 0, size = 8)
    )
  
  ggsave(file.path(out_dir_stats, "01_hm_enrichment_heatmap.png"), p1,
         width = plot_width, height = plot_height, dpi = 180, limitsize = FALSE)
  message("Saved -> 01_hm_enrichment_heatmap.png (",
          n_gene, " genes, ", n_ct, " cancer types)")
  
} else {
  message("Plot 1 skipped: no data.")
}

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 2 — Pan-cancer bubble plot
#
# Each point = one (gene, event, comparison) combination aggregated across
# all cancer types (using only rows with valid OR)
# x axis   = median OR across cancer types
# size     = number of cancer types where the gene is significantly enriched
# colour   = event type
# Orange triangle overlay = gene has >=1 cancer type with a zero-count cell
# ══════════════════════════════════════════════════════════════════════════════

bubble_df <- results_all %>%
  filter(
    IntoGen_cancer_type != "pan_cancer",
    test_type  == "per_event_pairwise",
    comparison %in% COMPARISONS_PAIR,
    !is.na(odds_ratio),                  # only valid OR for aggregation
    gene %in% top_genes_pancancer
  ) %>%
  drop_wt_wgd() %>%
  group_by(gene, event, comparison) %>%
  summarise(
    median_OR           = median(odds_ratio, na.rm = TRUE),
    n_sig_enriched      = sum(odds_ratio > 1 &
                                p_adj_local < SIG_THRESHOLD, na.rm = TRUE),
    any_zero_ct         = any(zero_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    event      = factor(event,      levels = EVENT_LEVELS),
    comparison = factor(comparison, levels = COMPARISONS_PAIR),
    gene       = factor(gene, levels = rev(top_genes_pancancer))
  )

if (nrow(bubble_df) > 0) {
  
  p2 <- ggplot(bubble_df,
               aes(x = median_OR, y = gene,
                   size = n_sig_enriched, colour = event)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               colour = "grey55", linewidth = 0.5) +
    geom_point(alpha = 0.82) +
    # orange open triangle for genes with zero-count in >=1 cancer type
    geom_point(data = filter(bubble_df, any_zero_ct),
               shape = 2, colour = "#e08214", size = 3, stroke = 0.8) +
    facet_grid(comparison ~ event, scales = "free_x") +
    scale_x_log10(
      breaks = c(0.25, 0.5, 1, 2, 5, 10, 20),
      labels = c("0.25","0.5","1","2","5","10","20"),
      name   = "Median OR (log scale)"
    ) +
    scale_size_continuous(
      name   = paste0("N cancer types\nenriched (p<", SIG_THRESHOLD, ")"),
      range  = c(1, 9),
      breaks = 0:max(bubble_df$n_sig_enriched)
    ) +
    scale_colour_manual(values = mutation_palette, name = "Event", drop = FALSE) +
    labs(
      title    = "Pan-cancer bubble plot — top genes",
      subtitle = paste0(
        "Top ", MAX_GENES, " genes selected for pan-cancer  |  ",
        "WT excluded from WGD comparisons  |  ",
        "orange triangle = zero-count in >=1 cancer type"
      ),
      y = NULL
    ) +
    base_theme +
    theme(
      axis.text.y  = element_text(face = "italic"),
      strip.text.y = element_text(angle = 0, hjust = 0, size = 8)
    )
  
  ggsave(file.path(out_dir_stats, "02_bubble_pancancer.png"), p2,
         width  = 14,
         height = max(8, n_distinct(bubble_df$gene) * 0.38 * 3 + 4),
         dpi = 180)
  message("Saved -> 02_bubble_pancancer.png")
  
} else {
  message("Plot 2 skipped: no data.")
}

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 3 — Pan-cancer forest plot
#
# Each row = one (gene, event, comparison) combination at pan-cancer level
# x axis   = OR on log scale
# CI bars  = 95% confidence interval from Fisher test
#            grey = normal, orange = zero_count with valid OR (case 2)
# Points   = filled circle (enriched) or open circle (depleted)
#            filled square = OR undefined (case 3)
# Labels   = significance stars; "NA" for undefined OR
# ══════════════════════════════════════════════════════════════════════════════

forest_df <- results_all %>%
  filter(
    IntoGen_cancer_type == "pan_cancer",
    test_type           == "per_event_pairwise",
    comparison          %in% COMPARISONS_PAIR,
    gene %in% top_genes_pancancer
  ) %>%
  drop_wt_wgd() %>%
  mutate(
    event      = factor(event,      levels = EVENT_LEVELS),
    comparison = factor(comparison, levels = COMPARISONS_PAIR),
    gene       = factor(gene, levels = rev(top_genes_pancancer)),
    label      = cell_label(odds_ratio, zero_count, p_adj_local),
    # shape encoding: filled circle, open circle, or square
    point_type = case_when(
      is.na(odds_ratio) ~ "absent",       # case 3: class has no samples
      odds_ratio >= 1   ~ "enriched",     # OR >= 1
      TRUE              ~ "depleted"      # OR < 1
    )
  )

# split into subsets for layered rendering
forest_normal  <- filter(forest_df, !is.na(odds_ratio), !zero_count)
forest_zero_ci <- filter(forest_df, !is.na(odds_ratio),  zero_count)
forest_na      <- filter(forest_df,  is.na(odds_ratio))

if (nrow(forest_df) > 0) {
  
  p3 <- ggplot(forest_df,
               aes(x = odds_ratio, y = gene, colour = event)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               colour = "grey40", linewidth = 0.6) +
    # CI bars: grey for normal rows
    geom_errorbarh(data = forest_normal,
                   aes(xmin = ci_low, xmax = ci_high),
                   height = 0.28, linewidth = 0.55, alpha = 0.65) +
    # CI bars: orange for zero-count rows with valid OR
    geom_errorbarh(data = forest_zero_ci,
                   aes(xmin = ci_low, xmax = ci_high),
                   colour = "#e08214",
                   height = 0.28, linewidth = 0.7, alpha = 0.85) +
    # points for rows with valid OR
    geom_point(data = filter(forest_df, !is.na(odds_ratio)),
               aes(shape = point_type), size = 2.4) +
    # square marker at x=1 for rows where OR is undefined
    geom_point(data = forest_na,
               aes(x = 1), shape = 15, size = 3.0) +
    # significance label at right end of CI for valid-OR rows
    geom_text(data = filter(forest_df, !is.na(odds_ratio)),
              aes(x = ci_high, label = label),
              hjust = -0.2, size = 2.3, colour = "grey20") +
    # "NA" label for undefined-OR rows
    geom_text(data = forest_na,
              aes(x = 1, label = "NA"),
              hjust = -0.3, size = 2.2, colour = "grey40") +
    facet_grid(comparison ~ event, scales = "free") +
    scale_x_log10(
      breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10, 20),
      labels = c("0.1","0.25","0.5","1","2","5","10","20"),
      limits = c(0.05, 200),
      oob    = scales::oob_squish,
      name   = "Odds ratio (log scale)"
    ) +
    scale_colour_manual(values = mutation_palette, name = "Event", drop = FALSE) +
    scale_shape_manual(
      values = c("enriched" = 19, "depleted" = 1, "absent" = 15),
      name   = "OR direction\n(square = class absent)"
    ) +
    labs(
      title    = "Pan-cancer forest plot — all pairwise comparisons",
      subtitle = paste0(
        "Top ", MAX_GENES, " genes  |  axis [0.05, 200]  |  ",
        "WT excluded from WGD comparisons  |  ",
        "orange CI = zero event carriers in one class  |  ",
        "square + NA = class absent for this gene"
      ),
      y = NULL
    ) +
    base_theme +
    theme(
      axis.text.y   = element_text(face = "italic", size = 8),
      panel.spacing = unit(0.4, "cm"),
      strip.text.y  = element_text(angle = 0, hjust = 0, size = 8)
    )
  
  n_genes_f <- n_distinct(forest_df$gene)
  ggsave(file.path(out_dir_stats, "03_forest_pancancer.png"), p3,
         width  = max(10, 4 * 3.0 + 3),
         height = max(8, n_genes_f * 0.42 * 3 + 4),
         dpi = 180, limitsize = FALSE)
  message("Saved -> 03_forest_pancancer.png")
  
} else {
  message("Plot 3 skipped: no data.")
}

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 4 — Event proportion bars (pan-cancer)
#
# Each facet = one gene (top pan-cancer genes)
# x axis    = sample class (WGD, HM, Classic)
# y axis    = proportion of samples in that class carrying each event
# n=        = number of samples above each bar
# seg=      = number of distinct CNA segments below each bar
# ══════════════════════════════════════════════════════════════════════════════

prop_df <- count_table %>%
  filter(
    IntoGen_cancer_type == "pan_cancer",
    gene %in% top_genes_pancancer
  ) %>%
  mutate(
    gene  = factor(gene,  levels = top_genes_pancancer),
    event = factor(event, levels = EVENT_LEVELS),
    class = factor(class, levels = CLASS_LEVELS)
  )

# sample counts per (gene, class) — one number per bar
n_samples_df <- prop_df %>%
  group_by(gene, class) %>%
  summarise(n_samples = sum(n), .groups = "drop") %>%
  mutate(
    gene  = factor(gene,  levels = top_genes_pancancer),
    class = factor(class, levels = CLASS_LEVELS),
    label = paste0("n=", n_samples)
  )

# segment counts per (gene, class) from the raw Drivers object
n_segments_df <- Drivers %>%
  filter(!is.na(class), class %in% CLASS_LEVELS,
         gene %in% top_genes_pancancer) %>%
  group_by(gene, class) %>%
  summarise(n_segments = n_distinct(segment_id), .groups = "drop") %>%
  mutate(
    gene  = factor(gene,  levels = top_genes_pancancer),
    class = factor(class, levels = CLASS_LEVELS),
    label = paste0("seg=", n_segments)
  )

if (nrow(prop_df) > 0) {
  
  p4 <- ggplot(prop_df, aes(x = class, y = prop, fill = event)) +
    geom_col(width = 0.78, position = "stack",
             colour = "white", linewidth = 0.2) +
    # proportion labels inside bars (only if >=7%)
    geom_text(data = filter(prop_df, prop >= 0.07),
              aes(label = percent(prop, accuracy = 1)),
              position = position_stack(vjust = 0.5),
              size = 1.8, colour = "white", fontface = "bold") +
    # n_samples above each bar
    geom_text(data = n_samples_df,
              aes(x = class, y = 1.04, label = label),
              inherit.aes = FALSE,
              size = 2.0, colour = "grey25", vjust = 0) +
    # n_segments below each bar
    geom_text(data = n_segments_df,
              aes(x = class, y = -0.07, label = label),
              inherit.aes = FALSE,
              size = 2.0, colour = "grey40", vjust = 1) +
    facet_wrap(~ gene, ncol = 6) +
    scale_fill_manual(values = mutation_palette, name = "Event", drop = FALSE) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(-0.12, 1.12),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_x_discrete(drop = FALSE) +
    labs(
      title    = paste0("Event proportions — top ", MAX_GENES,
                        " pan-cancer genes"),
      subtitle = "n = samples per class (above bar)  |  seg = CNA segments (below bar)",
      x = NULL, y = "Proportion of samples"
    ) +
    base_theme +
    theme(
      axis.text.x  = element_text(angle = 40, hjust = 1, size = 7),
      strip.text   = element_text(face = "bold.italic", size = 7.5),
      panel.spacing = unit(0.4, "cm")
    )
  
  n_rows <- ceiling(length(top_genes_pancancer) / 6)
  ggsave(file.path(out_dir_stats, "04_proportions_pancancer.png"), p4,
         width = 14, height = max(5, n_rows * 3.5 + 2), dpi = 180)
  message("Saved -> 04_proportions_pancancer.png")
  
} else {
  message("Plot 4 skipped: no data.")
}

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 5 — Pan-cancer pairwise heatmap
#
# Each tile = one (gene, comparison) combination, faceted by event
# Colour    = log2(OR) capped at +/-3
# Border    = orange if zero_count
# Label     = same as plot 1
# Gene ordering: significant genes on top, then zero-count, then by effect size
# ══════════════════════════════════════════════════════════════════════════════

pairwise_df <- results_all %>%
  filter(
    IntoGen_cancer_type == "pan_cancer",
    test_type           == "per_event_pairwise",
    comparison          %in% COMPARISONS_PAIR
  ) %>%
  drop_wt_wgd() %>%
  mutate(
    log2_OR    = ifelse(is.na(odds_ratio), NA_real_,
                        pmax(pmin(log2(odds_ratio), 3), -3)),
    label      = cell_label(odds_ratio, zero_count, p_adj_local),
    event      = factor(event,      levels = EVENT_LEVELS),
    comparison = factor(comparison, levels = COMPARISONS_PAIR)
  )

if (nrow(pairwise_df) > 0) {
  
  # order genes: sig on top, then zero-count, then by mean |log2OR|
  gene_order <- pairwise_df %>%
    group_by(gene) %>%
    summarise(
      any_sig      = any(!is.na(p_adj_local) & p_adj_local < SIG_THRESHOLD,
                         na.rm = TRUE),
      any_zero     = any(zero_count, na.rm = TRUE),
      mean_abs_lor = mean(abs(log2_OR), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(any_sig, any_zero, mean_abs_lor) %>%   # ascending so top of plot = most sig
    pull(gene)
  
  pairwise_df <- mutate(pairwise_df,
                        gene = factor(gene, levels = gene_order))
  
  p5 <- ggplot(pairwise_df,
               aes(x = comparison, y = gene, fill = log2_OR)) +
    geom_tile(data = filter(pairwise_df, !zero_count),
              colour = "white", linewidth = 0.3) +
    geom_tile(data = filter(pairwise_df, zero_count),
              colour = "#e08214", linewidth = 0.8) +
    geom_text(aes(label = label), size = 2.2, vjust = 0.75, colour = "grey15") +
    facet_wrap(~ event, nrow = 1) +
    scale_fill_gradientn(
      colours  = c("#2166ac", "#92c5de", "#f7f7f7", "#f4a582", "#d6604d"),
      limits   = c(-3, 3),
      na.value = "grey88",
      name     = "log2(OR)\n[+/-3]"
    ) +
    labs(
      title    = "Pan-cancer pairwise heatmap",
      subtitle = paste0(
        "OR > 0 = enriched in class 1 of comparison  |  ",
        "WT excluded from WGD comparisons  |  ",
        "orange border = zero count  |  ",
        "NA = class absent  |  0 = event absent in one class"
      ),
      x = NULL, y = NULL
    ) +
    base_theme +
    theme(
      axis.text.x  = element_text(angle = 30, hjust = 1, size = 8),
      axis.text.y  = element_text(face = "italic", size = 8),
      panel.spacing = unit(0.3, "cm")
    )
  
  ggsave(file.path(out_dir_stats, "05_heatmap_pancancer.png"), p5,
         width  = max(10, 3 * 2.5 + 4),
         height = max(4, n_distinct(pairwise_df$gene) * 0.38 + 3),
         dpi = 180, limitsize = FALSE)
  message("Saved -> 05_heatmap_pancancer.png")
  
} else {
  message("Plot 5 skipped: no data.")
}

message("\nPan-cancer plots done.")

# ══════════════════════════════════════════════════════════════════════════════
# PER-CANCER-TYPE PLOTS
# Plots 2-5 repeated for each tumor type using per-CT gene selection.
# ══════════════════════════════════════════════════════════════════════════════

cancer_types <- results_all %>%
  filter(IntoGen_cancer_type != "pan_cancer") %>%
  pull(IntoGen_cancer_type) %>%
  unique() %>% sort()

message("Running per-cancer-type plots for ", length(cancer_types), " tumor types ...")

for (ct in cancer_types) {
  
  ct_dir <- file.path(out_dir_stats, gsub("[^A-Za-z0-9_]", "_", ct))
  dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)
  
  # subset all objects to this cancer type
  res_ct  <- results_all %>% filter(IntoGen_cancer_type == ct)
  cnt_ct  <- count_table  %>% filter(IntoGen_cancer_type == ct)
  drv_ct  <- Drivers      %>% filter(IntoGen_cancer_type == ct)
  
  # gene selection specific to this cancer type
  top_ct <- select_top_genes(res_ct)
  
  if (length(top_ct) == 0) {
    message("  ", ct, " — no genes selected, skipping.")
    next
  }
  
  # ── PLOT 2: bubble plot ────────────────────────────────────────────────────
  bub_ct <- res_ct %>%
    filter(
      test_type  == "per_event_pairwise",
      comparison %in% COMPARISONS_PAIR,
      !is.na(odds_ratio),
      gene %in% top_ct
    ) %>%
    drop_wt_wgd() %>%
    mutate(
      event      = factor(event,      levels = EVENT_LEVELS),
      comparison = factor(comparison, levels = COMPARISONS_PAIR),
      gene       = factor(gene, levels = rev(top_ct)),
      label      = cell_label(odds_ratio, zero_count, p_adj_local)
    )
  
  if (nrow(bub_ct) > 0) {
    p2_ct <- ggplot(bub_ct, aes(x = odds_ratio, y = gene, colour = event)) +
      geom_vline(xintercept = 1, linetype = "dashed",
                 colour = "grey55", linewidth = 0.5) +
      geom_point(aes(size = odds_ratio), alpha = 0.85) +
      # orange open ring for zero-count points
      geom_point(data = filter(bub_ct, zero_count),
                 aes(size = odds_ratio),
                 shape = 1, colour = "#e08214", stroke = 1.0) +
      geom_text(aes(label = label),
                hjust = -0.4, size = 2.3, colour = "grey20",
                show.legend = FALSE) +
      facet_grid(comparison ~ event, scales = "free_x") +
      scale_x_log10(
        breaks = c(0.25, 0.5, 1, 2, 5, 10, 20),
        labels = c("0.25","0.5","1","2","5","10","20"),
        name   = "OR (log scale)"
      ) +
      scale_size_continuous(name = "OR", range = c(2, 8)) +
      scale_colour_manual(values = mutation_palette, name = "Event",
                          drop = FALSE) +
      labs(
        title    = paste0("Enrichment — ", ct),
        subtitle = paste0(
          "Top ", MAX_GENES, " genes  |  ",
          "WT excluded from WGD comparisons  |  ",
          "orange ring = zero event carriers in one class"
        ),
        y = NULL
      ) +
      base_theme +
      theme(
        axis.text.y  = element_text(face = "italic"),
        strip.text.y = element_text(angle = 0, hjust = 0, size = 8)
      )
    
    ggsave(file.path(ct_dir, "02_bubble.png"), p2_ct,
           width  = 14,
           height = max(5, n_distinct(bub_ct$gene) * 0.4 * 3 + 3),
           dpi = 180)
  }
  
  # ── PLOT 3: forest plot ────────────────────────────────────────────────────
  fst_ct <- res_ct %>%
    filter(
      test_type  == "per_event_pairwise",
      comparison %in% COMPARISONS_PAIR,
      gene %in% top_ct
    ) %>%
    drop_wt_wgd() %>%
    mutate(
      event      = factor(event,      levels = EVENT_LEVELS),
      comparison = factor(comparison, levels = COMPARISONS_PAIR),
      gene       = factor(gene, levels = rev(top_ct)),
      label      = cell_label(odds_ratio, zero_count, p_adj_local),
      point_type = case_when(
        is.na(odds_ratio) ~ "absent",
        odds_ratio >= 1   ~ "enriched",
        TRUE              ~ "depleted"
      )
    )
  
  if (nrow(fst_ct) > 0) {
    fst_ok <- filter(fst_ct, !is.na(odds_ratio), !zero_count)
    fst_zc <- filter(fst_ct, !is.na(odds_ratio),  zero_count)
    fst_na <- filter(fst_ct,  is.na(odds_ratio))
    
    p3_ct <- ggplot(fst_ct,
                    aes(x = odds_ratio, y = gene, colour = event)) +
      geom_vline(xintercept = 1, linetype = "dashed",
                 colour = "grey40", linewidth = 0.6) +
      geom_errorbarh(data = fst_ok,
                     aes(xmin = ci_low, xmax = ci_high),
                     height = 0.28, linewidth = 0.55, alpha = 0.65) +
      geom_errorbarh(data = fst_zc,
                     aes(xmin = ci_low, xmax = ci_high),
                     colour = "#e08214",
                     height = 0.28, linewidth = 0.7, alpha = 0.85) +
      geom_point(data = filter(fst_ct, !is.na(odds_ratio)),
                 aes(shape = point_type), size = 2.4) +
      geom_point(data = fst_na,
                 aes(x = 1), shape = 15, size = 3.0) +
      geom_text(data = filter(fst_ct, !is.na(odds_ratio)),
                aes(x = ci_high, label = label),
                hjust = -0.2, size = 2.3, colour = "grey20") +
      geom_text(data = fst_na,
                aes(x = 1, label = "NA"),
                hjust = -0.3, size = 2.2, colour = "grey40") +
      facet_grid(comparison ~ event, scales = "free") +
      scale_x_log10(
        breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10, 20),
        labels = c("0.1","0.25","0.5","1","2","5","10","20"),
        limits = c(0.05, 200),
        oob    = scales::oob_squish,
        name   = "Odds ratio (log scale)"
      ) +
      scale_colour_manual(values = mutation_palette, name = "Event",
                          drop = FALSE) +
      scale_shape_manual(
        values = c("enriched" = 19, "depleted" = 1, "absent" = 15),
        name   = "OR direction\n(square = class absent)"
      ) +
      labs(
        title    = paste0("Forest plot — ", ct),
        subtitle = paste0(
          "Top ", MAX_GENES, " genes  |  axis [0.05, 200]  |  ",
          "WT excluded from WGD comparisons  |  ",
          "orange CI = zero event carriers in one class  |  ",
          "square + NA = class absent"
        ),
        y = NULL
      ) +
      base_theme +
      theme(
        axis.text.y   = element_text(face = "italic", size = 8),
        panel.spacing = unit(0.4, "cm"),
        strip.text.y  = element_text(angle = 0, hjust = 0, size = 8)
      )
    
    ggsave(file.path(ct_dir, "03_forest.png"), p3_ct,
           width  = max(10, 4 * 3.0 + 3),
           height = max(6, n_distinct(fst_ct$gene) * 0.42 * 3 + 4),
           dpi = 180)
  }
  
  # ── PLOT 4: event proportions ──────────────────────────────────────────────
  prp_ct <- cnt_ct %>%
    filter(gene %in% top_ct) %>%
    mutate(
      gene  = factor(gene,  levels = top_ct),
      event = factor(event, levels = EVENT_LEVELS),
      class = factor(class, levels = CLASS_LEVELS)
    )
  
  if (nrow(prp_ct) > 0) {
    
    nsamp_ct <- prp_ct %>%
      group_by(gene, class) %>%
      summarise(n = sum(n), .groups = "drop") %>%
      mutate(
        gene  = factor(gene,  levels = top_ct),
        class = factor(class, levels = CLASS_LEVELS),
        label = paste0("n=", n)
      )
    
    nseg_ct <- drv_ct %>%
      filter(!is.na(class), class %in% CLASS_LEVELS, gene %in% top_ct) %>%
      group_by(gene, class) %>%
      summarise(n_seg = n_distinct(segment_id), .groups = "drop") %>%
      mutate(
        gene  = factor(gene,  levels = top_ct),
        class = factor(class, levels = CLASS_LEVELS),
        label = paste0("seg=", n_seg)
      )
    
    n_cols <- min(6, length(top_ct))
    
    p4_ct <- ggplot(prp_ct, aes(x = class, y = prop, fill = event)) +
      geom_col(width = 0.78, position = "stack",
               colour = "white", linewidth = 0.2) +
      geom_text(data = filter(prp_ct, prop >= 0.07),
                aes(label = percent(prop, accuracy = 1)),
                position = position_stack(vjust = 0.5),
                size = 1.8, colour = "white", fontface = "bold") +
      geom_text(data = nsamp_ct,
                aes(x = class, y = 1.04, label = label),
                inherit.aes = FALSE,
                size = 2.0, colour = "grey25", vjust = 0) +
      geom_text(data = nseg_ct,
                aes(x = class, y = -0.07, label = label),
                inherit.aes = FALSE,
                size = 2.0, colour = "grey40", vjust = 1) +
      facet_wrap(~ gene, ncol = n_cols) +
      scale_fill_manual(values = mutation_palette, name = "Event",
                        drop = FALSE) +
      scale_y_continuous(
        labels = percent_format(accuracy = 1),
        limits = c(-0.12, 1.12),
        expand = expansion(mult = c(0, 0))
      ) +
      scale_x_discrete(drop = FALSE) +
      labs(
        title    = paste0("Event proportions — top ", MAX_GENES,
                          " genes — ", ct),
        subtitle = "n = samples (above bar)  |  seg = CNA segments (below bar)",
        x = NULL, y = "Proportion of samples"
      ) +
      base_theme +
      theme(
        axis.text.x  = element_text(angle = 40, hjust = 1, size = 7),
        strip.text   = element_text(face = "bold.italic", size = 7.5),
        panel.spacing = unit(0.4, "cm")
      )
    
    n_rows <- ceiling(length(top_ct) / n_cols)
    ggsave(file.path(ct_dir, "04_proportions.png"), p4_ct,
           width = 14, height = max(4, n_rows * 3.5 + 2), dpi = 180)
  }
  
  # ── PLOT 5: pairwise heatmap ───────────────────────────────────────────────
  pw_ct <- res_ct %>%
    filter(
      test_type  == "per_event_pairwise",
      comparison %in% COMPARISONS_PAIR
    ) %>%
    drop_wt_wgd() %>%
    mutate(
      log2_OR    = ifelse(is.na(odds_ratio), NA_real_,
                          pmax(pmin(log2(odds_ratio), 3), -3)),
      label      = cell_label(odds_ratio, zero_count, p_adj_local),
      event      = factor(event,      levels = EVENT_LEVELS),
      comparison = factor(comparison, levels = COMPARISONS_PAIR)
    )
  
  if (nrow(pw_ct) > 0) {
    
    gene_ord_ct <- pw_ct %>%
      group_by(gene) %>%
      summarise(
        any_sig      = any(!is.na(p_adj_local) & p_adj_local < SIG_THRESHOLD,
                           na.rm = TRUE),
        any_zero     = any(zero_count, na.rm = TRUE),
        mean_abs_lor = mean(abs(log2_OR), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(any_sig, any_zero, mean_abs_lor) %>%
      pull(gene)
    
    pw_ct <- mutate(pw_ct, gene = factor(gene, levels = gene_ord_ct))
    
    p5_ct <- ggplot(pw_ct,
                    aes(x = comparison, y = gene, fill = log2_OR)) +
      geom_tile(data = filter(pw_ct, !zero_count),
                colour = "white", linewidth = 0.3) +
      geom_tile(data = filter(pw_ct, zero_count),
                colour = "#e08214", linewidth = 0.8) +
      geom_text(aes(label = label), size = 2.2, vjust = 0.75,
                colour = "grey15") +
      facet_wrap(~ event, nrow = 1) +
      scale_fill_gradientn(
        colours  = c("#2166ac", "#92c5de", "#f7f7f7", "#f4a582", "#d6604d"),
        limits   = c(-3, 3),
        na.value = "grey88",
        name     = "log2(OR)\n[+/-3]"
      ) +
      labs(
        title    = paste0("Pairwise heatmap — ", ct),
        subtitle = paste0(
          "OR > 0 = enriched in class 1  |  ",
          "WT excluded from WGD comparisons  |  ",
          "orange border = zero count  |  ",
          "NA = class absent  |  0 = event absent in one class"
        ),
        x = NULL, y = NULL
      ) +
      base_theme +
      theme(
        axis.text.x  = element_text(angle = 30, hjust = 1, size = 8),
        axis.text.y  = element_text(face = "italic", size = 8),
        panel.spacing = unit(0.3, "cm")
      )
    
    ggsave(file.path(ct_dir, "05_heatmap.png"), p5_ct,
           width  = max(10, 3 * 2.5 + 4),
           height = max(4, n_distinct(pw_ct$gene) * 0.38 + 3),
           dpi = 180)
  }
  
  message("  ", ct, " — done")
}

message("\nAll plots saved to: ", out_dir_stats)
