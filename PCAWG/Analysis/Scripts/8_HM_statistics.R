# ══════════════════════════════════════════════════════════════════════════════
#  Segment characteristics — WGD class analysis
#  Plots are produced:
#    1. Once globally  (all tumor types combined)
#    2. Once per IntoGen_cancer_type
#  Output: ~/Documents/material-tickTack-2026/PCAWG/Plot/HM segments distribution/
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(ggridges)
library(patchwork)
library(rstatix)
library(ggpubr)
library(ggsignif)

# ── Load data ─────────────────────────────────────────────────
Segments <- readRDS("~/Documents/material-tickTack-2026/PCAWG/Data/Segments.rds")

# ── Output directory ──────────────────────────────────────────
out_dir <- "~/Documents/material-tickTack-2026/PCAWG/Plot/HM segments distribution"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Colour palettes ───────────────────────────────────────────
pal_class <- c(Classic = "#4393C3", WGD = "#74C476", HM = "#D6604D")

# ── Chromosome level order (numeric, not alphabetical) ────────
# paste0("chr", 1:22) preserves 1,2,...,10,11,...,22 ordering.
# This single vector is reused everywhere chr is factored so that
# chr2 always sorts before chr10 regardless of R's default string sort.
chr_levels <- paste0("chr", 1:22)

# ── Theme ─────────────────────────────────────────────────────
theme_pub <- function() {
  theme_classic(base_size = 11) +
    theme(
      strip.background   = element_rect(fill = "#F0F0F0", colour = NA),
      strip.text         = element_text(face = "bold", size = 10),
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(size = 9, colour = "grey40"),
      axis.title         = element_text(size = 10),
      legend.position    = "bottom",
      panel.grid.major.y = element_line(colour = "grey92")
    )
}

# ── hg19/GRCh37 chromosome lengths (Mb) ──────────────────────
chr_lengths <- tibble(
  chr = chr_levels,
  chr_len_mb = c(
    249.250621, 243.199373, 198.022430, 191.154276, 180.915260,
    171.115067, 159.138663, 146.364022, 141.213431, 135.534747,
    135.006516, 133.851895, 115.169878, 107.349540, 102.531392,
    90.354753,  81.195210,  78.077248,  59.128983,  63.025520,
    48.129895,  51.304566
  )
)

# ── Pre-processing (global, done once) ───────────────────────
Segments <- Segments %>%
  mutate(
    seg_start  = as.numeric(str_extract(segment_name, "(?<=_)\\d+(?=_)")),
    seg_end    = as.numeric(str_extract(segment_name, "(?<=_\\d{1,12}_)\\d+")),
    seg_length = (seg_end - seg_start) / 1e6,
    chr        = factor(chr, levels = c(chr_levels, "chrX", "chrY")),
    seg_mid    = (seg_start + seg_end) / 2 / 1e6,
    wgd_status = class,
    class      = arm_level_class
  ) %>%
  filter(wgd_status %in% c("Classic", "WGD", "HM"))

# ── Clock rank (computed globally once, then filtered per subset) ─
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

Segments <- Segments %>% select(-c(clock_rank)) %>% 
  left_join(segment_event_ranks, by = c("sample", "clock_mean"))

# ── Axis label helpers ────────────────────────────────────────
log_len_breaks <- c(log10(0.1), log10(1), log10(10), log10(100))
log_len_labels <- c("0.1 Mb", "1 Mb", "10 Mb", "100 Mb")
log_snv_breaks <- c(log10(1), log10(10), log10(100), log10(1000), log10(10000), log10(1e5))
log_snv_labels <- c("1", "10", "100", "1K", "10K", "100K")


# ══════════════════════════════════════════════════════════════════════════════
#  CORE PLOTTING FUNCTION
#  Takes a (pre-filtered) Segments subset and a label string,
#  produces all three parts and saves them.
# ══════════════════════════════════════════════════════════════════════════════

run_segment_plots <- function(seg, label, out_dir) {
  
  # Minimum class diversity guard: skip if fewer than 2 classes present
  classes_present <- unique(na.omit(seg$wgd_status))
  if (length(classes_present) < 2) {
    message("  Skipping ", label, " — fewer than 2 WGD classes present.")
    return(invisible(NULL))
  }
  
  # Minimum segment guard
  if (nrow(seg) < 30) {
    message("  Skipping ", label, " — too few segments (n=", nrow(seg), ").")
    return(invisible(NULL))
  }
  
  n_samples <- n_distinct(seg$sample)
  message("  Running: ", label, "  (n segments=", nrow(seg),
          ", n samples=", n_samples, ")")
  
  seg_log <- seg %>%
    mutate(log_len = log10(seg_length),
           log_snv = log10(n_snvs + 1))
  
  # ── STATISTICS — Part 1: across classes ────────────────────
  
  kw_len  <- tryCatch(kruskal_test(seg_log, log_len ~ wgd_status), error = function(e) NULL)
  kw_snv  <- tryCatch(kruskal_test(seg_log, log_snv ~ wgd_status), error = function(e) NULL)
  
  dunn_len <- tryCatch(
    dunn_test(seg_log, log_len ~ wgd_status, p.adjust.method = "BH") %>%
      add_significance("p.adj"),
    error = function(e) NULL
  )
  dunn_snv <- tryCatch(
    dunn_test(seg_log, log_snv ~ wgd_status, p.adjust.method = "BH") %>%
      add_significance("p.adj"),
    error = function(e) NULL
  )
  
  max_log_len <- max(seg_log$log_len, na.rm = TRUE)
  max_log_snv <- max(seg_log$log_snv, na.rm = TRUE)
  step <- 0.12
  
  dunn_len_sig <- if (!is.null(dunn_len)) {
    dunn_len %>%
      filter(p.adj.signif != "ns") %>%
      arrange(p.adj) %>%
      mutate(y.position = max_log_len + step * row_number(),
             xmin = as.numeric(factor(group1, levels = names(pal_class))),
             xmax = as.numeric(factor(group2, levels = names(pal_class))))
  } else tibble()
  
  dunn_snv_sig <- if (!is.null(dunn_snv)) {
    dunn_snv %>%
      filter(p.adj.signif != "ns") %>%
      arrange(p.adj) %>%
      mutate(y.position = max_log_snv + step * row_number(),
             xmin = as.numeric(factor(group1, levels = names(pal_class))),
             xmax = as.numeric(factor(group2, levels = names(pal_class))))
  } else tibble()
  
  # ── STATISTICS — Part 2: within class × clock rank ─────────
  
  seg_log_ranked <- seg_log %>% filter(!is.na(clock_rank))
  
  # Only run if at least 2 rank levels are present within a class
  classes_with_multi_rank <- seg_log_ranked %>%
    group_by(wgd_status) %>%
    summarise(n_ranks = n_distinct(clock_rank), .groups = "drop") %>%
    filter(n_ranks >= 2) %>%
    pull(wgd_status)
  
  dunn_len_clock <- if (length(classes_with_multi_rank) > 0) {
    tryCatch(
      seg_log_ranked %>%
        filter(wgd_status %in% classes_with_multi_rank) %>%
        group_by(wgd_status) %>%
        dunn_test(log_len ~ clock_rank, p.adjust.method = "BH") %>%
        add_significance("p.adj") %>%
        ungroup(),
      error = function(e) NULL
    )
  } else NULL
  
  dunn_snv_clock <- if (length(classes_with_multi_rank) > 0) {
    tryCatch(
      seg_log_ranked %>%
        filter(wgd_status %in% classes_with_multi_rank) %>%
        group_by(wgd_status) %>%
        dunn_test(log_snv ~ clock_rank, p.adjust.method = "BH") %>%
        add_significance("p.adj") %>%
        ungroup(),
      error = function(e) NULL
    )
  } else NULL
  
  kw_len_clock <- tryCatch(
    seg_log_ranked %>%
      filter(wgd_status %in% classes_with_multi_rank) %>%
      group_by(wgd_status) %>%
      kruskal_test(log_len ~ clock_rank) %>%
      mutate(lab = paste0("KW p=", formatC(p, digits = 2, format = "e"))),
    error = function(e) NULL
  )
  
  kw_snv_clock <- tryCatch(
    seg_log_ranked %>%
      filter(wgd_status %in% classes_with_multi_rank) %>%
      group_by(wgd_status) %>%
      kruskal_test(log_snv ~ clock_rank) %>%
      mutate(lab = paste0("KW p=", formatC(p, digits = 2, format = "e"))),
    error = function(e) NULL
  )
  
  max_log_len_clk <- max(seg_log$log_len, na.rm = TRUE)
  max_log_snv_clk <- max(seg_log$log_snv, na.rm = TRUE)
  
  dunn_len_clk_sig <- if (!is.null(dunn_len_clock)) {
    dunn_len_clock %>%
      filter(p.adj.signif != "ns") %>%
      group_by(wgd_status) %>%
      arrange(p.adj) %>%
      mutate(y.position = max_log_len_clk + 0.10 * row_number()) %>%
      ungroup()
  } else tibble()
  
  dunn_snv_clk_sig <- if (!is.null(dunn_snv_clock)) {
    dunn_snv_clock %>%
      filter(p.adj.signif != "ns") %>%
      group_by(wgd_status) %>%
      arrange(p.adj) %>%
      mutate(y.position = max_log_snv_clk + 0.10 * row_number()) %>%
      ungroup()
  } else tibble()
  
  # ── PART 1 ─────────────────────────────────────────────────
  
  kw_len_p  <- if (!is.null(kw_len))  formatC(kw_len$p,  digits = 2, format = "e") else "NA"
  kw_snv_p  <- if (!is.null(kw_snv))  formatC(kw_snv$p,  digits = 2, format = "e") else "NA"
  
  p1a <- ggplot(seg_log, aes(x = wgd_status, y = log_len, fill = wgd_status)) +
    geom_violin(alpha = 0.55, linewidth = 0.4, trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
                 linewidth = 0.5, coef = 0) +
    {if (nrow(dunn_len_sig) > 0)
      geom_signif(
        data = dunn_len_sig,
        aes(xmin = xmin, xmax = xmax,
            y_position = y.position, annotations = p.adj.signif),
        manual = TRUE, tip_length = 0.01, textsize = 3.5, inherit.aes = FALSE
      )} +
    scale_y_continuous(breaks = log_len_breaks, labels = log_len_labels) +
    scale_fill_manual(values = pal_class, drop = FALSE) +
    labs(title    = "Segment length by class",
         subtitle = paste0("KW p=", kw_len_p, " | Dunn BH-corrected brackets"),
         x = NULL, y = "Segment length", fill = NULL) +
    theme_pub() + theme(legend.position = "none")
  
  p1b <- ggplot(seg_log, aes(x = wgd_status, y = log_snv, fill = wgd_status)) +
    geom_violin(alpha = 0.55, linewidth = 0.4, trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
                 linewidth = 0.5, coef = 0) +
    {if (nrow(dunn_snv_sig) > 0)
      geom_signif(
        data = dunn_snv_sig,
        aes(xmin = xmin, xmax = xmax,
            y_position = y.position, annotations = p.adj.signif),
        manual = TRUE, tip_length = 0.01, textsize = 3.5, inherit.aes = FALSE
      )} +
    scale_y_continuous(breaks = log_snv_breaks, labels = log_snv_labels) +
    scale_fill_manual(values = pal_class, drop = FALSE) +
    labs(title    = "n SNVs per segment by class",
         subtitle = paste0("KW p=", kw_snv_p, " | Dunn BH-corrected brackets"),
         x = NULL, y = "n SNVs", fill = NULL) +
    theme_pub() + theme(legend.position = "none")
  
  seg_summary <- seg %>%
    group_by(wgd_status) %>%
    summarise(mean_len = mean(seg_length, na.rm = TRUE),
              sd_len   = sd(seg_length,   na.rm = TRUE),
              mean_snv = mean(n_snvs,     na.rm = TRUE),
              sd_snv   = sd(n_snvs,       na.rm = TRUE),
              n = n(), .groups = "drop")
  
  p1c <- ggplot(seg_summary, aes(x = wgd_status, y = mean_len, fill = wgd_status)) +
    geom_col(width = 0.5, alpha = 0.85) +
    geom_errorbar(aes(ymin = mean_len - sd_len, ymax = mean_len + sd_len),
                  width = 0.18, linewidth = 0.8) +
    geom_text(aes(label = paste0("n=", scales::comma(n))),
              vjust = -0.5, size = 3.2) +
    scale_fill_manual(values = pal_class, drop = FALSE) +
    labs(title = "Mean segment length ± SD", x = NULL, y = "Length (Mb)", fill = NULL) +
    theme_pub() + theme(legend.position = "none")
  
  p1d <- ggplot(seg_summary, aes(x = wgd_status, y = mean_snv, fill = wgd_status)) +
    geom_col(width = 0.5, alpha = 0.85) +
    geom_errorbar(aes(ymin = pmax(0, mean_snv - sd_snv), ymax = mean_snv + sd_snv),
                  width = 0.18, linewidth = 0.8) +
    scale_fill_manual(values = pal_class, drop = FALSE) +
    labs(title = "Mean n SNVs ± SD", x = NULL, y = "n SNVs", fill = NULL) +
    theme_pub() + theme(legend.position = "none")
  
  p1e <- ggplot(seg, aes(x = seg_mid, y = chr, colour = wgd_status)) +
    geom_jitter(height = 0.25, alpha = 0.08, size = 0.5) +
    stat_summary(aes(group = wgd_status), fun = median, geom = "point",
                 shape = 21, fill = "white", size = 2.5, stroke = 1) +
    scale_colour_manual(values = pal_class, drop = FALSE) +
    scale_y_discrete(limits = rev) +
    facet_wrap(~wgd_status) +
    labs(title = "Genomic position of segments by class",
         x = "Chromosomal position (Mb)", y = NULL, colour = NULL) +
    theme_pub() + theme(legend.position = "none")
  
  part1 <- (p1a | p1b) / (p1c | p1d) / p1e +
    plot_annotation(
      title    = paste0("Part 1 — Segment characteristics by WGD class  [", label, "]"),
      subtitle = "Significance: * p<0.05  ** p<0.01  *** p<0.001  (Dunn, BH-corrected)",
      theme    = theme(plot.title = element_text(face = "bold", size = 14))
    )
  
  # ── PART 2 ─────────────────────────────────────────────────
  
  p2a <- ggplot(seg %>% filter(!is.na(clock_rank)),
                aes(x = seg_length, y = clock_rank,
                    fill = wgd_status, colour = wgd_status)) +
    geom_density_ridges(alpha = 0.45, scale = 0.9,
                        quantile_lines = TRUE, quantiles = 2,
                        rel_min_height = 0.01) +
    scale_x_log10(labels = scales::label_number(suffix = " Mb")) +
    scale_fill_manual(values = pal_class, drop = FALSE) +
    scale_colour_manual(values = pal_class, drop = FALSE) +
    facet_wrap(~wgd_status, ncol = 3) +
    labs(title    = "Segment length by within-sample clock rank",
         subtitle = "1 = earliest event in sample | 3+ = third or later",
         x = "Length (Mb, log scale)", y = "Within-sample rank",
         fill = NULL, colour = NULL) +
    theme_pub()
  
  p2b <- ggplot(seg %>% filter(!is.na(clock_rank)),
                aes(x = n_snvs + 1, y = clock_rank,
                    fill = wgd_status, colour = wgd_status)) +
    geom_density_ridges(alpha = 0.45, scale = 0.9,
                        quantile_lines = TRUE, quantiles = 2,
                        rel_min_height = 0.01) +
    scale_x_log10(labels = scales::label_number()) +
    scale_fill_manual(values = pal_class, drop = FALSE) +
    scale_colour_manual(values = pal_class, drop = FALSE) +
    facet_wrap(~wgd_status, ncol = 3) +
    labs(title    = "n SNVs by within-sample clock rank",
         subtitle = "1 = earliest event in sample | 3+ = third or later",
         x = "n SNVs + 1 (log scale)", y = "Within-sample rank",
         fill = NULL, colour = NULL) +
    theme_pub()
  
  p2c_box <- ggplot(seg_log %>% filter(!is.na(clock_rank)),
                    aes(x = clock_rank, y = log_len, fill = wgd_status)) +
    geom_violin(alpha = 0.45, linewidth = 0.3, trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white",
                 linewidth = 0.4, coef = 0) +
    {if (nrow(dunn_len_clk_sig) > 0)
      stat_pvalue_manual(dunn_len_clk_sig,
                         label = "p.adj.signif", xmin = "group1", xmax = "group2",
                         y.position = "y.position", tip.length = 0.01,
                         size = 2.8, hide.ns = TRUE)} +
    {if (!is.null(kw_len_clock))
      geom_text(data = kw_len_clock,
                aes(x = 1, y = Inf, label = lab, colour = wgd_status),
                vjust = 1.5, hjust = 0, size = 2.8, inherit.aes = FALSE)} +
    scale_y_continuous(breaks = log_len_breaks, labels = log_len_labels) +
    scale_fill_manual(values = pal_class, drop = FALSE) +
    scale_colour_manual(values = pal_class, drop = FALSE) +
    facet_wrap(~wgd_status, ncol = 3, scales = "free_y") +
    labs(title    = "Segment length across within-sample clock ranks",
         subtitle = "Label = KW p per class | brackets = Dunn BH-corrected",
         x = "Within-sample rank", y = "Segment length", fill = NULL) +
    theme_pub() + theme(legend.position = "none")
  
  p2d_box <- ggplot(seg_log %>% filter(!is.na(clock_rank)),
                    aes(x = clock_rank, y = log_snv, fill = wgd_status)) +
    geom_violin(alpha = 0.45, linewidth = 0.3, trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white",
                 linewidth = 0.4, coef = 0) +
    {if (nrow(dunn_snv_clk_sig) > 0)
      stat_pvalue_manual(dunn_snv_clk_sig,
                         label = "p.adj.signif", xmin = "group1", xmax = "group2",
                         y.position = "y.position", tip.length = 0.01,
                         size = 2.8, hide.ns = TRUE)} +
    {if (!is.null(kw_snv_clock))
      geom_text(data = kw_snv_clock,
                aes(x = 1, y = Inf, label = lab, colour = wgd_status),
                vjust = 1.5, hjust = 0, size = 2.8, inherit.aes = FALSE)} +
    scale_y_continuous(breaks = log_snv_breaks, labels = log_snv_labels) +
    scale_fill_manual(values = pal_class, drop = FALSE) +
    scale_colour_manual(values = pal_class, drop = FALSE) +
    facet_wrap(~wgd_status, ncol = 3, scales = "free_y") +
    labs(title    = "n SNVs across within-sample clock ranks",
         subtitle = "Label = KW p per class | brackets = Dunn BH-corrected",
         x = "Within-sample rank", y = "n SNVs", fill = NULL) +
    theme_pub() + theme(legend.position = "none")
  
  part2 <- (p2a / p2b) / (p2c_box / p2d_box) +
    plot_annotation(
      title    = paste0("Part 2 — Segment characteristics by class × clock rank  [", label, "]"),
      subtitle = paste(
        "1 = earliest CNA in sample | 2 = second | 3+ = third or later",
        "KW = Kruskal-Wallis | Dunn BH-corrected (* <0.05, ** <0.01, *** <0.001)",
        sep = "\n"
      ),
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
  
  # ── PART 3 ─────────────────────────────────────────────────
  
  n_samples_per_cell <- seg %>%
    filter(!is.na(clock_rank)) %>%
    group_by(wgd_status, clock_rank) %>%
    summarise(n_samples = n_distinct(sample), .groups = "drop")
  
  seg_cov <- seg %>%
    filter(!is.na(clock_rank),
           as.character(chr) %in% chr_lengths$chr) %>%
    group_by(wgd_status, clock_rank, chr) %>%
    summarise(total_seg_mb = sum(seg_length, na.rm = TRUE),
              n_seg = n(), .groups = "drop") %>%
    left_join(chr_lengths,        by = "chr") %>%
    left_join(n_samples_per_cell, by = c("wgd_status", "clock_rank")) %>%
    mutate(
      coverage_density = (total_seg_mb / chr_len_mb) / n_samples,
      chr = factor(chr, levels = chr_levels)
    )
  
  # Only keep rank levels actually present in this subset
  ranks_present <- levels(droplevels(seg$clock_rank[!is.na(seg$clock_rank)]))
  
  full_grid <- expand_grid(
    wgd_status = classes_present,
    clock_rank = factor(ranks_present, levels = rank_level_order),
    chr        = factor(chr_levels, levels = chr_levels)
  ) %>%
    left_join(chr_lengths, by = "chr") %>%
    left_join(seg_cov %>% select(wgd_status, clock_rank, chr,
                                 total_seg_mb, n_seg, coverage_density),
              by = c("wgd_status", "clock_rank", "chr")) %>%
    mutate(across(c(total_seg_mb, n_seg, coverage_density), ~replace_na(.x, 0)))
  
  wilcox_results <- full_grid %>%
    group_by(wgd_status, clock_rank) %>%
    mutate(row_median = median(coverage_density, na.rm = TRUE)) %>%
    group_by(wgd_status, clock_rank, chr) %>%
    summarise(
      density    = unique(coverage_density),
      row_median = unique(row_median),
      p_value    = tryCatch(
        wilcox.test(coverage_density, mu = unique(row_median),
                    exact = FALSE, conf.int = FALSE)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    group_by(wgd_status) %>%
    mutate(p_adj  = p.adjust(p_value, method = "BH"),
           signif = case_when(p_adj < 0.001 ~ "***",
                              p_adj < 0.01  ~ "**",
                              p_adj < 0.05  ~ "*",
                              TRUE          ~ "")) %>%
    ungroup()
  
  plot_data <- full_grid %>%
    left_join(wilcox_results %>% select(wgd_status, clock_rank, chr, p_adj, signif),
              by = c("wgd_status", "clock_rank", "chr")) %>%
    group_by(wgd_status, clock_rank) %>%
    mutate(row_median = median(coverage_density, na.rm = TRUE),
           log2_fc    = log2((coverage_density + 1e-4) / (row_median + 1e-4))) %>%
    ungroup()
  
  p3a <- ggplot(plot_data,
                aes(x = clock_rank,
                    y = fct_rev(factor(chr, levels = chr_levels)),
                    fill = coverage_density)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = signif), colour = "black", size = 2.8, vjust = 0.75) +
    scale_fill_gradientn(
      colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
      name    = "Mean coverage\ndensity / sample",
      labels  = scales::label_percent(accuracy = 0.1),
      limits  = c(0, NA)
    ) +
    facet_wrap(~wgd_status, ncol = 3) +
    labs(title    = "Chromosome coverage density by class × clock rank",
         subtitle = "Colour = sum(seg length) / chr length / n_samples | Stars = Wilcoxon vs row median (BH)",
         x = "Within-sample rank", y = NULL) +
    theme_pub() +
    theme(axis.text.y = element_text(size = 7),
          axis.text.x = element_text(angle = 30, hjust = 1),
          legend.key.height = unit(0.5, "cm"),
          legend.key.width  = unit(1.2, "cm"))
  
  p3b <- ggplot(plot_data,
                aes(x = clock_rank,
                    y = fct_rev(factor(chr, levels = chr_levels)),
                    fill = log2_fc)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = signif), colour = "grey20", size = 2.8, vjust = 0.75) +
    scale_fill_gradient2(
      low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0,
      name = "log2(density /\nrow median)",
      limits = c(-3, 3), oob = scales::squish
    ) +
    facet_wrap(~wgd_status, ncol = 3) +
    labs(title    = "Relative coverage density (log2 ratio vs row median)",
         subtitle = "Red = depleted | Blue = enriched relative to class × rank median",
         x = "Within-sample rank", y = NULL) +
    theme_pub() +
    theme(axis.text.y = element_text(size = 7),
          axis.text.x = element_text(angle = 30, hjust = 1),
          legend.key.height = unit(0.5, "cm"),
          legend.key.width  = unit(1.2, "cm"))
  
  part3 <- p3a / p3b +
    plot_annotation(
      title    = paste0("Part 3 — Chromosome coverage density  [", label, "]"),
      subtitle = "Normalised by chromosome length and n samples per cell",
      theme    = theme(plot.title = element_text(face = "bold", size = 14))
    )
  
  # ── Save stats CSV ─────────────────────────────────────────
  safe_label <- gsub("[^A-Za-z0-9_]", "_", label)
  
  bind_rows(
    if (!is.null(dunn_len))       dunn_len       %>% mutate(metric = "seg_length", contrast = "across_classes")      else NULL,
    if (!is.null(dunn_snv))       dunn_snv       %>% mutate(metric = "n_snvs",     contrast = "across_classes")      else NULL,
    if (!is.null(dunn_len_clock)) dunn_len_clock %>% mutate(metric = "seg_length", contrast = "within_class_clock")  else NULL,
    if (!is.null(dunn_snv_clock)) dunn_snv_clock %>% mutate(metric = "n_snvs",     contrast = "within_class_clock")  else NULL
  ) %>%
    write_csv(file.path(out_dir, paste0("stats_", safe_label, ".csv")))
  
  wilcox_results %>%
    arrange(wgd_status, clock_rank, chr) %>%
    write_csv(file.path(out_dir, paste0("coverage_wilcox_", safe_label, ".csv")))
  
  # ── Save PDFs ──────────────────────────────────────────────
  ggsave(file.path(out_dir, paste0("part1_class_", safe_label, ".pdf")),
         part1, width = 12, height = 14)
  ggsave(file.path(out_dir, paste0("part2_clockrank_", safe_label, ".pdf")),
         part2, width = 14, height = 26)
  ggsave(file.path(out_dir, paste0("part3_coverage_", safe_label, ".pdf")),
         part3, width = 14, height = 18)
  
  message("  ✔  Saved 3 PDFs + 2 CSVs for: ", label)
  invisible(NULL)
}


# ══════════════════════════════════════════════════════════════════════════════
#  RUN 1 — GLOBAL (all tumor types combined)
# ══════════════════════════════════════════════════════════════════════════════

message("\n════ GLOBAL (all tumor types) ════")
run_segment_plots(seg = Segments, label = "ALL", out_dir = out_dir)


# ══════════════════════════════════════════════════════════════════════════════
#  RUN 2 — PER IntoGen_cancer_type
# ══════════════════════════════════════════════════════════════════════════════

cancer_types <- sort(unique(na.omit(Segments$IntoGen_cancer_type)))
message("\n════ Per-tumor-type loop: ", length(cancer_types), " types ════")

for (ct in cancer_types) {
  seg_ct <- Segments %>% filter(IntoGen_cancer_type == ct)
  run_segment_plots(seg = seg_ct, label = ct, out_dir = out_dir)
}

message("\n✔  All done. Files saved to:\n  ", out_dir)