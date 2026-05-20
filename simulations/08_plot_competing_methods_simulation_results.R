#!/usr/bin/env Rscript

# "tickTack H produces better clustering than applying Mclust post-hoc to the outputs of AmpTimeR, 
# MutTimeR, and the no-clustering baseline, across a range of simulation conditions."

rm(list = ls())

# ── setup ─────────────────────────────────────────────────────────────────────
base_dir <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/simulations/"
.libPaths("~/R/orfeo_R_4.4/")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(ggbeeswarm)
  library(patchwork)
})

# ── load data ─────────────────────────────────────────────────────────────────
results_all <- readRDS(
  paste0(base_dir, "results_summary/07_results_competing_methods_BIC.rds")
)

# ── Nature Genetics theme ─────────────────────────────────────────────────────
# Helvetica-like, minimal axes, no background, small but readable
theme_ng <- function(base_size = 12) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      # axes
      axis.line        = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = unit(2, "pt"),
      axis.title       = element_text(size = base_size, face = "plain"),
      axis.text        = element_text(size = base_size - 1, colour = "black"),
      axis.text.x      = element_text(angle = 0, hjust = 1),
      # strips
      strip.background = element_blank(),
      strip.text       = element_text(size = base_size, face = "bold"),
      strip.placement  = "outside",
      # legend
      legend.position  = "right",
      legend.key.size  = unit(6, "pt"),
      legend.text      = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size - 1, face = "bold"),
      legend.background = element_blank(),
      # plot
      plot.title       = element_text(size = base_size + 1, face = "bold"),
      plot.subtitle    = element_text(size = base_size - 1, colour = "grey40"),
      panel.spacing    = unit(3, "pt"),
      plot.margin      = margin(4, 4, 4, 4, "pt")
    )
}

# Nature Genetics palette — distinct, print-safe, colourblind-friendly
method_colors <- c(
  "ticktack_Single" = "#377EB8",   # blue
  "AmpTimeR"        = "#E41A1C",   # red
  "MutTimeR"        = "#FF7F00"    # orange
)

method_labels <- c(
  "ticktack_Single" = "tickTack (no clust.)",
  "AmpTimeR"        = "AmpTimeR",
  "MutTimeR"        = "MutTimeR"
)

ref_method      <- "tickTack H"
methods_compare <- c("ticktack_Single", "AmpTimeR", "MutTimeR")

# ── data prep ─────────────────────────────────────────────────────────────────
df <- results_all %>%
  # filter(status == "OK", mutation_density != 2e-07) %>%
  mutate(
    method           = recode(method, "Var-BIC" = "tickTack H"),
    n_events         = as.numeric(n_events),
    mutation_density = round(as.numeric(mutation_density), 8)
  )

md_levels <- sort(unique(df$mutation_density))
df <- df %>%
  mutate(
    md_label = factor(
      paste0("mut. density: ", mutation_density),
      levels = paste0("mut. density: ", md_levels)
    )
  )

# ── helper: delta + wilcoxon for any metric ───────────────────────────────────
make_delta <- function(df, metric) {
  
  df_wide <- df %>%
    select(sim_name, iter, method, n_events, md_label, value = all_of(metric)) %>%
    pivot_wider(names_from = method, values_from = value)
  
  df_delta <- bind_rows(lapply(methods_compare, function(m) {
    tibble(
      sim_name = df_wide$sim_name,
      iter     = df_wide$iter,
      n_events = df_wide$n_events,
      md_label = df_wide$md_label,
      method   = m,
      delta    = df_wide[[ref_method]] - df_wide[[m]]
    )
  })) %>%
    filter(!is.na(delta)) %>%
    mutate(method = factor(method, levels = methods_compare))
  
  df_summary <- df_delta %>%
    group_by(method, n_events, md_label) %>%
    summarise(
      mean_delta = mean(delta),
      se         = sd(delta) / sqrt(n()),
      .groups    = "drop"
    ) %>%
    mutate(method = factor(method, levels = methods_compare))
  
  wilcox_res <- df_delta %>%
    group_by(method, n_events, md_label) %>%
    summarise(
      p = tryCatch(
        wilcox.test(delta, mu = 0, exact = FALSE)$p.value,
        error   = function(e) NA,
        warning = function(w) wilcox.test(delta, mu = 0, exact = FALSE)$p.value
      ),
      .groups = "drop"
    ) %>%
    mutate(
      p_adj  = p.adjust(p, method = "BH"),
      signif = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE          ~ ""
      )
    )
  
  y_star <- df_summary %>%
    group_by(n_events, md_label) %>%
    summarise(
      y = max(mean_delta + 2 * se, na.rm = TRUE) + diff(range(df_delta$delta, na.rm = TRUE)) * 0.06,
      .groups = "drop"
    )
  
  wilcox_res <- wilcox_res %>%
    left_join(y_star, by = c("n_events", "md_label")) %>%
    mutate(method = factor(method, levels = methods_compare))
  
  list(delta = df_delta, summary = df_summary, wilcox = wilcox_res)
}

# ── helper: build panel ───────────────────────────────────────────────────────
make_panel <- function(d, title, ylab) {
  
  ggplot(d$delta, aes(x = method, y = delta, colour = method)) +
    
    geom_hline(yintercept = 0, linetype = "dashed",
               linewidth = 0.3, colour = "grey60") +
    
    geom_beeswarm(size = 0.55, alpha = 0.2, cex = 0.9, show.legend = FALSE) +
    
    geom_errorbar(
      data    = d$summary,
      aes(y = mean_delta, ymin = mean_delta - 2 * se, ymax = mean_delta + 2 * se),
      width   = 0.15, linewidth = 0.6, colour = "black"
    ) +
    geom_point(
      data   = d$summary,
      aes(y  = mean_delta),
      size   = 2.8, shape = 21, colour = "black", fill = "white", stroke = 0.8
    ) +
    
    geom_text(
      data         = d$wilcox %>% filter(signif != ""),
      aes(x = method, y = y, label = signif),
      colour       = "black", size = 4.2, inherit.aes = FALSE
    ) +
    
    scale_colour_manual(values = method_colors, labels = method_labels,
                        name = NULL) +
    scale_x_discrete(labels = method_labels) +
    
    facet_grid(md_label ~ n_events, labeller = labeller(n_events = label_both),
               switch = "y") +
    
    labs(title = title, x = NULL, y = ylab) +
    
    theme_ng() +
    theme(legend.position = "bottom")
}

# ── panel A: delta adj. Rand Index ────────────────────────────────────────────
d_ari <- make_delta(df, "adj_RI")

pA <- make_panel(
  d_ari,
  title = "a  Clustering accuracy",
  ylab  = expression(Delta~"adj. Rand Index  (tickTack H \u2212 method)")
)

# ── shared legend ─────────────────────────────────────────────────────────────
legend_plot <- ggplot(
  tibble(method = factor(methods_compare, levels = methods_compare),
         x = 1, y = 1),
  aes(x = x, y = y, colour = method)
) +
  geom_point(size = 2) +
  scale_colour_manual(values = method_colors, labels = method_labels, name = NULL) +
  theme_ng() +
  theme(legend.position = "right",
        legend.key.size = unit(8, "pt"),
        legend.text     = element_text(size = 6))

leg <- cowplot::get_legend(legend_plot)

# ── assemble ──────────────────────────────────────────────────────────────────
# patchwork: stack A over B, attach legend on the right
library(cowplot)

combined <- plot_grid(
  plot_grid(pA, ncol = 1, align = "v", labels = NULL),
  # leg,
  ncol        = 1,
  rel_widths  = c(1)
)

# ── save ──────────────────────────────────────────────────────────────────────
n_col <- length(unique(df$n_events))
n_row <- length(unique(df$md_label))

ggsave(
  paste0(base_dir, "plot/08_ticktackH_adjRI_correlation_delta.pdf"),
  combined,
  width  = 15.2 * n_col,   # ~183 mm max for NG single column = 89 mm, double = 183 mm
  height = 8.0 * n_row ,     # two stacked panels
  units  = "cm",
  device = cairo_pdf             # embeds Helvetica correctly
)

message("Saved: plot/08_ticktackH_adjRI_correlation_delta.pdf")





########################## option 2 to visualize comparison #######################
# ── setup ─────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(patchwork)
})

# ── load data ─────────────────────────────────────────────────────────────────
results_all <- readRDS(
  paste0(base_dir, "results_summary/07_results_competing_methods_BIC.rds")
)

# ── shared theme ──────────────────────────────────────────────────────────────
theme_main <- function() {
  theme_bw(base_size = 13) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey88"),
      panel.grid.minor   = element_blank(),
      strip.background   = element_rect(fill = "grey94", colour = NA),
      strip.text         = element_text(face = "bold", size = 10),
      axis.title         = element_text(face = "bold", size = 12),
      axis.text          = element_text(colour = "black", size = 9),
      axis.text.x        = element_text(angle = 0, hjust = 0.5, size = 8),
      legend.position    = "right",
      legend.title       = element_text(face = "bold", size = 10),
      plot.title         = element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle      = element_text(size = 9, hjust = 0, colour = "grey40"),
      panel.spacing      = unit(0.6, "lines")
    )
}

# ── palettes ──────────────────────────────────────────────────────────────────
method_colors <- c(
  "ticktack_Single" = "#777EB8",
  "AmpTimeR"        = "#E41A1C",
  "MutTimeR"        = "#AE7F00"
)

method_labels <- c(
  "ticktack_Single" = "tickTack (no clust.)",
  "AmpTimeR"        = "AmpTimeR",
  "MutTimeR"        = "MutTimeR"
)

ref_method      <- "tickTack H"
methods_compare <- c("ticktack_Single", "AmpTimeR", "MutTimeR")

# ── data prep ─────────────────────────────────────────────────────────────────
df <- results_all %>%
  mutate(
    method           = recode(method, "Var-BIC" = "tickTack H"),
    n_events         = as.numeric(n_events),
    mutation_density = round(as.numeric(mutation_density), 8)
  )

md_levels <- sort(unique(df$mutation_density), decreasing = TRUE)
ne_levels  <- sort(unique(df$n_events))

df <- df %>%
  mutate(
    md_factor = factor(
      mutation_density,
      levels = md_levels,
      labels = formatC(md_levels, format = "e", digits = 0)
    ),
    ne_factor = factor(n_events, levels = ne_levels),
    is_hard   = mutation_density <= quantile(md_levels, 0.4)
  )

# ── hard-zone shading ─────────────────────────────────────────────────────────
n_hard    <- sum(md_levels <= quantile(md_levels, 0.4))
hard_rect <- data.frame(
  xmin = length(md_levels) - n_hard + 1.5,
  xmax = length(md_levels) + 0.5,
  ymin = -Inf,
  ymax = Inf
)

hard_bg <- function() {
  geom_rect(
    data        = hard_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill        = "#fff0e6",
    alpha       = 0.7
  )
}

palette_ne <- setNames(
  colorRampPalette(c("#f4a97f", "#084594"))(length(ne_levels)),
  as.character(ne_levels)
)

# ── delta + wilcoxon ──────────────────────────────────────────────────────────
make_delta <- function(df, metric) {
  
  df_wide <- df %>%
    select(sim_name, iter, method, n_events, md_factor, ne_factor,
           is_hard, value = all_of(metric)) %>%
    pivot_wider(names_from = method, values_from = value)
  
  df_delta <- bind_rows(lapply(methods_compare, function(m) {
    tibble(
      sim_name  = df_wide$sim_name,
      iter      = df_wide$iter,
      n_events  = df_wide$n_events,
      md_factor = df_wide$md_factor,
      ne_factor = df_wide$ne_factor,
      is_hard   = df_wide$is_hard,
      method    = m,
      delta     = df_wide[[ref_method]] - df_wide[[m]]
    )
  })) %>%
    filter(!is.na(delta)) %>%
    mutate(method = factor(method, levels = methods_compare))
  
  # Wilcoxon test per method x md_factor
  wilcox_res <- df_delta %>%
    group_by(method, md_factor) %>%
    summarise(
      p = tryCatch(
        wilcox.test(delta, mu = 0, exact = FALSE)$p.value,
        error   = function(e) NA_real_,
        warning = function(w) wilcox.test(delta, mu = 0, exact = FALSE)$p.value
      ),
      y_max = max(delta, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      p_adj  = p.adjust(p, method = "BH"),
      signif = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE          ~ ""
      ),
      # place stars just above the top of each violin
      y_star = y_max + diff(range(df_delta$delta, na.rm = TRUE)) * 0.04,
      method = factor(method, levels = methods_compare)
    )
  
  list(delta = df_delta, wilcox = wilcox_res)
}

# ── violin panel builder ───────────────────────────────────────────────────────
make_violin_delta <- function(data, y_lab, ref_line = 0, title = NULL) {
  
  ggplot(data$delta, aes(x = md_factor, y = delta)) +
    hard_bg() +
    geom_violin(
      aes(fill = method),
      alpha     = 0.15,
      linewidth = 0.3,
      trim      = TRUE
    ) +
    geom_jitter(
      aes(colour = ne_factor),
      width = 0.08, size = 0.8, alpha = 0.2
    ) +
    geom_hline(
      yintercept = ref_line,
      linetype   = "dashed",
      colour     = "grey40",
      linewidth  = 0.5
    ) +
    # significance stars
    geom_text(
      data         = data$wilcox %>% filter(signif != ""),
      aes(x = md_factor, y = y_star, label = signif),
      colour       = "black",
      size         = 4.5,
      inherit.aes  = FALSE
    ) +
    scale_fill_manual(
      values = method_colors,
      labels = method_labels,
      name   = "Method"
    ) +
    scale_colour_manual(
      values = palette_ne,
      name   = "N events",
      guide  = guide_legend(
        order        = 2,
        override.aes = list(size = 2, alpha = 1)
      )
    ) +
    facet_wrap(
      ~ method,
      nrow     = 1,
      labeller = labeller(method = method_labels)
    ) +
    labs(x = "Mutation density", y = y_lab, title = title) +
    theme_main()
}

# ── build panel ───────────────────────────────────────────────────────────────
d_ari <- make_delta(df, "adj_RI")

pA <- make_violin_delta(
  data  = d_ari,
  y_lab = expression(Delta ~ "adj. Rand Index  (tickTack H \u2212 method)"),
  title = "a  Clustering accuracy"
)

# ── assemble + save ───────────────────────────────────────────────────────────
combined <- pA +
  plot_layout(guides = "collect") +
  plot_annotation(theme = theme(legend.position = "right"))

ggsave(
  paste0(base_dir, "plot/08_ticktackH_adjRI_delta_violin.pdf"),
  combined,
  width  = 10 * length(methods_compare),
  height = 14,
  units  = "cm",
  device = cairo_pdf
)

message("Saved: plot/08_ticktackH_adjRI_delta_violin.pdf")



