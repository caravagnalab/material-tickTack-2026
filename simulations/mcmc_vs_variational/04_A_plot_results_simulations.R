base_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack_2026/simulations/"
#!/usr/bin/env Rscript
.libPaths("~/R/orfeo_R_4.4/")
library(dplyr)
library(ggplot2)
library(patchwork)


results_all <- readRDS(paste0(base_dir,"results_summary/results_all_simulations.rds"))
results_all_var <- readRDS(paste0(base_dir,"results_summary/results_all_simulations_variational.rds"))

theme_plot <- function() {
  theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_line(size = 0.2, colour = "grey85"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90", colour = NA),
      strip.text = element_text(face = "bold", size = 11),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
}

ALPHA = 0.8
palette_crit <- c(
  "BIC" = alpha("#1b9e77", alpha = ALPHA),
  "AIC" = alpha("#d95f02", alpha = ALPHA),
  "ICL" = alpha("#7570b3", alpha = ALPHA),
  "LOO" = alpha("#e7298a", alpha = ALPHA)
)



######## Model selection: retrieval of true number of clusters ############
results_plot <- results_all %>%
  filter(status %in% c("OK", "ONLY_VARIATIONAL")) %>%   # keep usable results
  mutate(
    inference = factor(inference, levels = c("var", "mcmc")),
    criterion = factor(criterion, levels = c("BIC", "AIC", "ICL", "LOO")),
    n_clusters = as.numeric(n_clusters),
    n_events = as.numeric(n_events)
  )

ggplot(results_plot, aes(x = n_clusters, y = K)) +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_grid(inference + criterion ~ mutation_density + n_events) +
  labs(
    title = "Inferred vs True Number of Clusters",
    x = "True number of clusters",
    y = "Inferred number of clusters"
  ) +
  theme_bw()

p1 <- ggplot(results_plot, aes(x = n_clusters, y = K, color = criterion)) +
  geom_jitter(width = 0.15, height = 0.15, alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(values = palette_crit) +
  facet_grid(inference ~ mutation_density + n_events, labeller = label_both) +
  labs(
    title = "Recovery of the Number of Clusters",
    subtitle = "Dashed line indicates perfect recovery (K inferred = K true)",
    x = "True number of clusters",
    y = "Inferred number of clusters"
  ) +
  theme_plot()

p1

###### Parameter accuracy: correlation plot #####
results_corr <- results_plot %>%
  filter(K > 1, !is.na(correlation), n_clusters > 1)

ggplot(results_corr, aes(x = factor(n_clusters), y = correlation)) +
  geom_boxplot(outlier.alpha = 0.3) +
  facet_grid(inference + criterion ~ mutation_density + n_events, scales = "free") +
  labs(
    title = "Correlation between True and Inferred τ",
    x = "True number of clusters",
    y = "Correlation"
  ) +
  theme_bw()

p2 <- ggplot(results_corr, aes(x = factor(n_clusters), y = correlation, fill = criterion)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.2) +
  scale_fill_manual(values = palette_crit) +
  facet_grid(inference ~ mutation_density + n_events, labeller = label_both) +
  labs(
    title = "Correlation between True and Inferred Timing (τ)",
    subtitle = "Only simulations with inferred K > 1",
    x = "True number of clusters",
    y = "Correlation"
  ) +
  theme_plot()

p2

###### Clustering quality: Adj RI plot ######
ggplot(results_plot, aes(x = factor(n_clusters), y = adj_RI)) +
  geom_boxplot(outlier.alpha = 0.3) +
  facet_grid(inference + criterion ~ mutation_density + n_events) +
  labs(
    title = "Adjusted Rand Index",
    x = "True number of clusters",
    y = "Adjusted Rand Index"
  ) +
  theme_bw()

p3 <- ggplot(results_plot, aes(x = factor(n_clusters), y = adj_RI, fill = criterion)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.2) +
  scale_fill_manual(values = palette_crit) +
  facet_grid(inference ~ mutation_density + n_events, labeller = label_both) +
  labs(
    title = "Clustering Performance (Adjusted Rand Index)",
    x = "True number of clusters",
    y = "Adjusted Rand Index"
  ) +
  theme_plot()

p3

p_final <- (p1 / p2 / p3) + plot_annotation(tag_levels = "A")
ggplot2::ggsave(paste0(base_dir,"/plot/results_sim_all_scores_fullinference.png"),plot = p_final, width = 10, height = 20)


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
  facet_grid(inference ~ mutation_density + n_events, labeller = label_both) +
  labs(
    title = "Accuracy of Cluster Number Recovery",
    x = "True number of clusters",
    y = "Proportion correctly inferred"
  ) +
  theme_plot()

p4
