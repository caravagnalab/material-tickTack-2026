base_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack_2026/simulations/"
#!/usr/bin/env Rscript
.libPaths("~/R/orfeo_R_4.4/")
library(ggplot2)
library(dplyr)
library(viridis)
library(patchwork)

results_all_K <- readRDS(paste0(base_dir, "results_summary/03_B_results_comptime_simulations_per_K.rds"))

# mean time
mean_time_K <- results_all_K %>%
  filter(!is.na(K)) %>% 
  filter(inference != "mcmc") %>%
  group_by(inference, K, n_events, mutation_density) %>%
  summarise(mean_time = mean(time_seconds, na.rm = TRUE), .groups = "drop")

p1 <- ggplot(mean_time_K, aes(x = factor(K), y = mean_time)) +
  geom_col(position = "dodge") +
  facet_grid(mutation_density ~ n_events, scales = "free_y") +
  scale_fill_viridis(discrete = TRUE, option = "E") +
  labs(
    x = "Inferred K",
    y = "Mean time (s)",
    fill = "Method",
    title = "Mean Time per K for Variational vs MCMC"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold")
  )
p1

# failed K inference 
failed_K_summary <- results_all_K %>%
  filter(failed == TRUE, !is.na(K)) %>%
  group_by(inference, K, n_clusters, n_events) %>%
  summarise(failed_count = n(), .groups = "drop") %>%
  group_by(inference, n_clusters, n_events) %>%
  mutate(failed_percent = 100 * failed_count / sum(failed_count)) %>%
  ungroup()

p2 <- ggplot(failed_K_summary, aes(x = factor(K), y = failed_percent, fill = inference)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = failed_count), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3) +
  facet_grid(n_clusters ~ n_events, scales = "free_y") +
  scale_fill_viridis(discrete = TRUE, option = "E") +
  labs(
    x = "K that failed",
    y = "Percentage of failures (%)",
    fill = "Method",
    title = "Percentage of failed K values per method"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold")
  )

p2

#  Distribution of total times per inference method
total_time <- results_all_K %>%
  filter(inference != "mcmc") %>%
  group_by(sim_name, iter, inference, n_events, mutation_density) %>%
  summarise(total_time = sum(time_seconds, na.rm = TRUE), .groups = "drop")

p3 <- ggplot(total_time, aes(x = total_time, fill = inference)) +
  geom_histogram(alpha = 0.8, bins = 20, position = "dodge") +
  facet_grid(mutation_density ~ n_events, scales = "free") +
  scale_fill_viridis(discrete = TRUE, option = "E") +
  labs(
    x = "Total time per iteration (s)",
    y = "Count",
    fill = "Method",
    title = "Distribution of Total Inference Time"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold")
  )

p3

p_final <- (p1 / p2 / p3) + plot_annotation(tag_levels = "a")
ggplot2::ggsave(paste0(base_dir,"/plot/04_B_results_comptime.png"),plot = p_final, width = 10, height = 20)

