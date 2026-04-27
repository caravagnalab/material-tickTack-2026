rm(list=ls())
base_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/simulations/"
#!/usr/bin/env Rscript
.libPaths("~/R/orfeo_R_4.4/")
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)



results_all <- readRDS(paste0(base_dir, "results_summary/03_A_results_simulations.rds"))

ALPHA = 0.8
palette_crit <- c(
  "BIC" = "#1b9e77",
  "AIC" = "#d95f02",
  "ICL" = "#7570b3",
  "LOO" = "#e7298a"
)

theme_plot <- function() {
  theme_bw(base_size = 13) +
    theme(
      panel.grid.major = element_line(linewidth = 0.2, colour = "grey85"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text = element_text(face = "bold", size = 10),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, colour = "grey40")
    )
}

results_plot <- results_all %>%
  filter(status %in% c("OK", "ONLY_VARIATIONAL")) %>%
  filter(mutation_density != 2e-07) %>%
  filter(criterion %in% c("BIC", "AIC")) %>%
  mutate(
    criterion         = factor(criterion, levels = c("BIC", "AIC", "ICL", "LOO")),
    n_clusters        = as.numeric(n_clusters),
    n_events          = as.numeric(n_events),
    mutation_density  = round(as.numeric(mutation_density), 8),
    K                 = as.numeric(K),
    correct_K         = as.integer(K == n_clusters)
  )

# ── friendly labels ────────────────────────────────────────────────────────────
md_labels <- sort(unique(results_plot$mutation_density))
ne_labels  <- sort(unique(results_plot$n_events))

results_plot <- results_plot %>%
  mutate(
    md_label = factor(
      paste0("mut. density: ", mutation_density),
      levels = paste0("mut. density: ", md_labels)
    ),
    ne_label = factor(
      paste0("n events: ", n_events),
      levels = paste0("n events: ", ne_labels)
    )
  )


ggplot(results_plot, aes(x = n_clusters, y = K)) +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_grid(purity + coverage ~ mutation_density + n_events) +
  labs(
    title = "Inferred vs true number of clusters",
    x = "True number of clusters",
    y = "Inferred number of clusters"
  ) +
  theme_bw()

p1 <- ggplot(results_plot, aes(x = n_clusters, y = K, color = criterion)) +
  geom_jitter(width = 0.15, height = 0.15, alpha = 0.4, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  scale_color_manual(values = palette_crit) +
  facet_grid(n_events ~ mutation_density, labeller = label_both) +
  labs(
    title = "Recovery of the number of clusters",
    # subtitle = "Dashed line indicates perfect recovery (K inferred = K true)",
    x = "True number of clusters",
    y = "Inferred number of clusters"
  ) +
  theme_plot()

p1


# ══════════════════════════════════════════════════════════════════════════════
# PLOT A option 2 –  recovery heatmap
# ══════════════════════════════════════════════════════════════════════════════
# ── compute signed error ──────────────────────────────────────────────────────
recovery_heat <- results_plot %>%
  group_by(criterion, n_clusters, md_label, ne_label) %>%
  summarise(
    recovery_rate = mean(K == n_clusters, na.rm = TRUE),
    mean_error    = mean(K - n_clusters, na.rm = TRUE),   # signed: negative = underestimate
    n_sim         = n(),
    .groups       = "drop"
  )

# ── symmetric color scale centered at 0 ──────────────────────────────────────
max_err <- max(abs(recovery_heat$mean_error), na.rm = TRUE)

p1 <- ggplot(recovery_heat,
             aes(x = factor(n_clusters), y = criterion, fill = mean_error)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(
    aes(label = sprintf("%+.1f\n(%s)", mean_error,
                        scales::percent(recovery_rate, accuracy = 1))),
    size = 2.5, lineheight = 1.1, colour = "grey15"
  ) +
  scale_fill_gradient2(
    low      = "#d95f02",   # orange-red  → underestimation
    mid      = "#f7f7f7",   # white       → perfect
    high     = "#1b9e77",   # teal-green  → overestimation
    midpoint = 0,
    limits   = c(-max_err, max_err),
    name     = "Mean\nK error"
  ) +
  facet_grid(ne_label ~ md_label) +
  labs(
    title    = "Bias in cluster number estimation",
    # subtitle = "Cell: mean signed error (K inferred \u2212 K true) · parentheses: exact recovery rate",
    x        = "True number of clusters",
    y        = NULL
  ) +
  theme_plot() +
  theme(
    legend.position  = "right",
    axis.text.y      = element_text(face = "bold"),
    panel.grid       = element_blank()
  )

p1

###### Parameter accuracy: correlation plot #####
results_corr <- results_plot %>%
  filter(K > 1, !is.na(correlation), n_clusters > 1)

p2 <- ggplot(results_corr,
             aes(x = factor(n_clusters), y = correlation, fill = criterion)) +
  # geom_boxplot(alpha = 0.85, outlier.alpha = 0.25, outlier.size = 0.8,
  #              linewidth = 0.4) +
  scale_fill_manual(values = palette_crit) +
  facet_grid(ne_label ~ md_label) +
  labs(
    title    = "Correlation between true and inferred timing",
    # subtitle = "Only simulations with inferred K > 1",
    x        = "True number of clusters",
    y        = "Pearson correlation"
  ) +
  theme_plot()

p2

# ══════════════════════════════════════════════════════════════════════════════
# PLOT C  –  ARI boxplot
# ══════════════════════════════════════════════════════════════════════════════
p3 <- ggplot(results_plot,
             aes(x = factor(n_clusters), y = adj_RI, fill = criterion)) +
  # geom_boxplot(alpha = 0.85, outlier.alpha = 0.25, outlier.size = 0.8,
  #              linewidth = 0.4) +
  scale_fill_manual(values = palette_crit) +
  facet_grid(ne_label ~ md_label) +
  labs(
    title = "Clustering performance (Adjusted Rand Index)",
    x     = "True number of clusters",
    y     = "Adjusted Rand Index"
  ) +
  theme_plot()

p3

# ══════════════════════════════════════════════════════════════════════════════
# COMBINED FIGURE  (a / b / c stacked)
# ══════════════════════════════════════════════════════════════════════════════
p_final <- (p1 / p2 / p3) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(
  paste0(base_dir, "/plot/04_A_results_sim_all_scores_v2.pdf"),
  plot   = p_final,
  width  = 14,
  height = 24
)

###### percentage of correct true k #####
results_plot %>%
  group_by(inference, criterion, n_clusters, mutation_density, n_events) %>%
  summarise(acc = mean(correct_K, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = factor(n_clusters), y = acc)) +
  geom_col(position = "dodge") +
  facet_grid(inference + criterion ~ mutation_density + n_events) +
  labs(
    title = "Accuracy of K Recovery",
    x = "True number of clusters",
    y = "Proportion correct"
  ) +
  theme_bw()

p4 <- results_plot %>%
  group_by(inference, criterion, n_clusters, mutation_density, n_events) %>%
  summarise(acc = mean(correct_K, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = factor(n_clusters), y = acc, fill = criterion)) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = palette_crit) +
  facet_grid(n_events ~ mutation_density, labeller = label_both) +
  labs(
    title = "Accuracy of Cluster Number Recovery",
    x = "True number of clusters",
    y = "Proportion correctly inferred"
  ) +
  theme_plot()

p4










# ── prepare data ──────────────────────────────────────────────────────────────
results_bic <- results_all %>%
  filter(status %in% c("OK", "ONLY_VARIATIONAL")) %>%
  filter(criterion == "BIC",
         mutation_density != 2e-07) %>%
  mutate(
    n_clusters       = as.numeric(n_clusters),
    n_events         = as.numeric(n_events),
    mutation_density = round(as.numeric(mutation_density), 8),
    K                = as.numeric(K),
    error            = K - n_clusters           # signed error column
  )

# ordered factors
md_levels <- sort(unique(results_bic$mutation_density),decreasing = TRUE)
ne_levels  <- sort(unique(results_bic$n_events))

results_bic <- results_bic %>%
  mutate(
    md_factor = factor(mutation_density,
                       levels = md_levels,
                       labels = formatC(md_levels, format = "e", digits = 0)),
    ne_factor = factor(n_events, levels = ne_levels),
    is_hard   = mutation_density <= quantile(md_levels, 0.4)
  )

# ── palettes ──────────────────────────────────────────────────────────────────
palette_ne   <- setNames(
  colorRampPalette(c("#f4a97f", "#084594"))(length(ne_levels)),
  as.character(ne_levels)
)
violin_fill  <- c("TRUE" = "#f4a97f", "FALSE" = "#084594")

# hard-zone shading rectangle (x positions of hard conditions)
n_hard    <- sum(md_levels <= quantile(md_levels, 0.4))
hard_rect <- data.frame(xmin = n_hard + 1.5, xmax = n_hard + 2.5, ymin = -Inf, ymax = Inf)

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

# ── violin helper ─────────────────────────────────────────────────────────────
hard_bg <- function() {
  geom_rect(data = hard_rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "#fff0e6", alpha = 0.7)
}

make_violin <- function(data, y_var, y_lab, ref_line = NULL) {
  p <- ggplot(data, aes(x = md_factor, y = .data[[y_var]])) +
    hard_bg() +
    geom_violin(aes(fill = is_hard),
                alpha = 0.45, linewidth = 0.3, trim = TRUE) +
    # geom_boxplot(aes(fill = is_hard),
    #              width = 0.15, outlier.shape = NA,
    #              linewidth = 0.35, alpha = 0.85) +
    geom_jitter(aes(colour = ne_factor),
                width = 0.08, size = 0.8, alpha = 0.2) +
    scale_fill_manual(
      values = violin_fill,
      name   = "Condition",
      labels = c("TRUE" = "Hard", "FALSE" = "Easy"),
      guide  = guide_legend(order = 1)
    ) +
    scale_colour_manual(
      values = palette_ne,
      name   = "N events",
      guide  = guide_legend(order = 2,
                            override.aes = list(size = 2, alpha = 1))
    ) +
    facet_wrap(~ n_clusters, nrow = 1,
               labeller = labeller(n_clusters = function(x) paste0("K = ", x))) +
    labs(x = "Mutation density", y = y_lab) +
    theme_main()
  
  if (!is.null(ref_line))
    p <- p + geom_hline(yintercept = ref_line, linetype = "dashed",
                        colour = "grey40", linewidth = 0.5)
  p
}

# ── panel a: signed error ─────────────────────────────────────────────────────
pa <- make_violin(
  data     = results_bic,
  y_var    = "error",
  y_lab    = "K inferred \u2212 K true",
  ref_line = 0
) +
  labs(
    title    = "a  Cluster number estimation",
    # subtitle = "Signed error (K inferred \u2212 K true) \u00b7 dashed line: unbiased recovery"
  )

# ── panel b: tau correlation ──────────────────────────────────────────────────
pb <- make_violin(
  data     = results_bic %>% filter(n_clusters > 1, K > 1, !is.na(correlation)),
  y_var    = "correlation",
  y_lab    = "Correlation",
  ref_line = 0.5
) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    title    = "b  Timing accuracy (correlation between tau real and tau inferred)",
    # subtitle = "Pearson correlation between true and inferred \u03c4 \u00b7 dashed line: 0.5"
  )

# ── panel c: ARI ──────────────────────────────────────────────────────────────
pc <- make_violin(
  data     = results_bic %>% filter(!is.na(adj_RI)),
  y_var    = "adj_RI",
  y_lab    = "Adjusted Rand Index",
  ref_line = 0.5
) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    title    = "c  Clustering quality (Adjusted Rand Index)",
    # subtitle = "ARI between true and inferred assignments \u00b7 dashed line: 0.5"
  )

# ── assemble ──────────────────────────────────────────────────────────────────
p_main <- (pa / pb / pc) +
  plot_layout(guides = "collect") +
  plot_annotation(
    # caption = "Orange shading: hard conditions (lowest mutation densities). Fill: hard (orange) vs easy (green) conditions. Points coloured by number of mutational events per cluster.",
    theme   = theme(
      plot.caption    = element_text(size = 9, colour = "grey40", hjust = 0),
      legend.position = "right"
    )
  )

p_main


ggsave(
  paste0(base_dir, "/plot/04_main_figure_BIC_violin_v3.pdf"),
  plot   = p_main,
  width  = 12,
  height = 10
)
