rm(list=ls())
.libPaths("~/R/orfeo_R_4.4/")
library(tickTack)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Distributions were truncated at the 10th–90th percentiles to reduce the influence of extreme outliers
base_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/simulations/"
data_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/simulations/data/"

info_parameters_PCAWG <- readRDS(paste0(data_dir,"/00_A_info_parameters_singleseg_PCAWG_BeforeSmoothing.rds"))
info_parameters_PCAWG <- info_parameters_PCAWG %>% filter(n_mut > 3 ) %>% 
  group_by(sample_id) %>% 
  mutate(n_seg = n()) %>% ungroup()
#                                                           len > 1e+06)
info_parameters_PCAWG %>% dplyr::select(sample_id) %>% distinct()


# info_parameters_PCAWG <- info_parameters_PCAWG %>% select(sample_id, median_DP, purity, median_seg_length, median_n_mut, n_seg) %>% 
#   dplyr::filter(purity > 0.4)

df_long <- info_parameters_PCAWG %>% 
  # dplyr::filter(purity > 0.4, median_DP >= 40) %>%
  dplyr::select(-c(sample_id,segment_id)) %>% pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

percentile_bounds <- df_long %>%
  group_by(variable) %>%
  summarise(
    p10 = quantile(value, 0.10, na.rm = TRUE),
    p90 = quantile(value, 0.90, na.rm = TRUE),
    .groups = "drop"
  )


df_trimmed <- df_long %>%
  left_join(percentile_bounds, by = "variable") %>%
  filter(value >= p10, value <= p90)

# df_trimmed <- df_long

stats_df <- df_trimmed %>%
  group_by(variable) %>%
  summarise(
    min = min(value, na.rm = TRUE),
    q1 = quantile(value, 0.25, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q3 = quantile(value, 0.75, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  )


p_main <- ggplot(df_trimmed, aes(x = variable, y = value)) +
  
  geom_boxplot(fill = "#4C72B0", alpha = 0.3, width = 0.5, outlier.shape = NA) +
  
  # Median
  geom_crossbar(
    data = stats_df,
    aes(y = median, ymin = median, ymax = median),
    color = "darkred", width = 0.5, linetype = "dashed"
  ) +
  
  # Q1 / Q3
  geom_crossbar(data = stats_df,
                aes(y = q1, ymin = q1, ymax = q1),
                linetype = "dashed", width = 0.3) +
  
  geom_crossbar(data = stats_df,
                aes(y = q3, ymin = q3, ymax = q3),
                linetype = "dashed", width = 0.3) +
  
  # Min / Max points
  # geom_point(data = stats_df, aes(y = min), shape = 16, size = 2) +
  # geom_point(data = stats_df, aes(y = max), shape = 16, size = 2) +
  
  # Labels
  geom_text(data = stats_df,
            aes(y = median, label = paste0("M=", signif(median, 3))),
            color = "darkred", vjust = -1, size = 3) +
  
  # geom_text(data = stats_df,
  #           aes(y = q1, label = paste0("Q1=", signif(q1, 3))),
  #           vjust = 1.5, size = 2.8) +
  
  # geom_text(data = stats_df,
  #           aes(y = q3, label = paste0("Q3=", signif(q3, 3))),
  #           vjust = -1, size = 2.8) +
  # 
  # geom_text(data = stats_df,
  #           aes(y = min, label = paste0("Min=", signif(min, 3))),
  #           vjust = 1.5, size = 2.5) +
  # 
  # geom_text(data = stats_df,
  #           aes(y = max, label = paste0("Max=", signif(max, 3))),
  #           vjust = -1, size = 2.5) +
  
  facet_wrap(~variable, scales = "free", ncol = 3) +
  
  labs(
    title = "Distribution of PCAWG samples parameters (min number of mutations = 3, trimmed at 10%-90%)",
    subtitle = "Min, Q1, Median, Q3, Max shown",
    x = "",
    y = "Value"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


p_main <- p_main +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(10, 30, 10, 10)
  )



p_density <- ggplot(df_trimmed, aes(x = value)) +
  
  geom_density(fill = "#4C72B0", alpha = 0.3) +
  
  # Quartiles
  # geom_vline(data = stats_df, aes(xintercept = q1), linetype = "dashed") +
  geom_vline(data = stats_df, aes(xintercept = median),
             color = "darkred", linewidth = 1, linetype = "dashed") +
  # geom_vline(data = stats_df, aes(xintercept = q3), linetype = "dashed") +
  
  # Min / Max
  # geom_vline(data = stats_df, aes(xintercept = min), linetype = "dotted") +
  # geom_vline(data = stats_df, aes(xintercept = max), linetype = "dotted") +
  
  # Labels
  # geom_text(data = stats_df,
  #           aes(x = min, y = Inf, label = paste0("Min=", signif(min, 3))),
  #           vjust = 2, size = 2.8) +
  # 
  # geom_text(data = stats_df,
  #           aes(x = q1, y = Inf, label = paste0("Q1=", signif(q1, 3))),
  #           vjust = 3.5, size = 3) +
  # 
  geom_text(data = stats_df,
            aes(x = median, y = Inf, label = paste0("M=", signif(median, 3))),
            color = "darkred", vjust = 5, size = 3) +
  
  # geom_text(data = stats_df,
  #           aes(x = q3, y = Inf, label = paste0("Q3=", signif(q3, 3))),
  #           vjust = 3.5, size = 3) +
  # 
  # geom_text(data = stats_df,
  #           aes(x = max, y = Inf, label = paste0("Max=", signif(max, 3))),
  #           vjust = 2, size = 2.8) +
  
  facet_wrap(~variable, scales = "free", ncol = 3) +
  
  labs(
    title = "Distribution of PCAWG samples parameters (min number of mutations = 3, trimmed at 10%-90%)",
    subtitle = "Dotted = Min/Max, Dashed = Quartiles, Red = Median",
    x = "Value",
    y = "Density"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p_density

p_density <- p_density +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(10, 30, 10, 10)
  )

############# check for correlation between number of mutations and number of segments #################
plot(info_parameters_PCAWG$n_seg, log(info_parameters_PCAWG$median_n_mut))
fit <- lm(log(median_n_mut) ~ n_seg, data = info_parameters_PCAWG)
abline(fit, col = "red", lwd = 2)

plot(info_parameters_PCAWG$n_seg, log(info_parameters_PCAWG$median_density_mut))
info_parameters_PCAWG$purity_group <- ifelse(info_parameters_PCAWG$median_seg_length > 0.6e+07, "High", "Low")
info_parameters_PCAWG$purity_group <- as.factor(info_parameters_PCAWG$purity_group)
ggplot(info_parameters_PCAWG, aes(x = n_seg, 
                                  y = log(median_density_mut), 
                                  color = purity_group)) +
  geom_point(alpha = 0.6) +                  # Scatter plot points
  geom_smooth(method = "lm", se = FALSE) +   # Fits and draws the two separate lines
  labs(title = "Bimodal Regression Lines by Purity",
       x = "n_seg",
       y = "log(median_density_mut)",
       color = "Purity Cluster") +
  theme_minimal()

plot(info_parameters_PCAWG$n_seg, info_parameters_PCAWG$purity)

plot(info_parameters_PCAWG$n_seg, info_parameters_PCAWG$median_DP)
plot(info_parameters_PCAWG$n_seg, info_parameters_PCAWG$median_seg_length)
fit <- lm((median_seg_length) ~ n_seg, data = info_parameters_PCAWG)
abline(fit, col = "red", lwd = 2)

plot(info_parameters_PCAWG$median_seg_length, log(info_parameters_PCAWG$median_density_mut))
fit <- lm(log(median_density_mut) ~ median_seg_length, data = info_parameters_PCAWG)
abline(fit, col = "red", lwd = 2)



# p_final <- (p_main / p_density ) + plot_annotation(tag_levels = "a")
p_final <- ( p_density ) + plot_annotation(tag_levels = "a")
ggplot2::ggsave(paste0(base_dir,"/plot/00_B_PCAWG_params_BeforeSmoothing_trimmed.pdf"),plot = p_final, width = 12, height = 8)
