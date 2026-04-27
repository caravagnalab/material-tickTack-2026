# Proportion of ttypers

# Average number of events and events' length by ttype



# TMB plot
rm(list=ls())
library(tidyverse)
install.packages("ggplot2")
library(ggplot2)
library(patchwork)
library(dplyr)

# TMB is typically calculated as mutations per megabase
Samples <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Samples.rds")
GENOME_SIZE_MB <- 3000

tmb_data <- Samples %>%
  filter(!is.na(IntoGen_cancer_type)) %>%
  mutate(
    n_snvs_filtered = map_int(cnaqc_path, function(path) {
      x <- readRDS(path)
      x$snvs %>%
        filter(mutation_type == "single base substitution") %>%  # adjust if needed
        nrow()
    }),
    TMB = n_snvs_filtered / GENOME_SIZE_MB
  )

saveRDS(tmb_data, "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/TMB_data.rds")


tmb_data <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/TMB_data.rds")

tmb_data <- tmb_data %>%
  group_by(IntoGen_cancer_type) %>%
  mutate(median_snvs = median(n_snvs)) %>%
  ungroup() %>%
  arrange((median_snvs)) %>%
  mutate(IntoGen_cancer_type = factor(IntoGen_cancer_type, levels = unique(IntoGen_cancer_type)))


Segments <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Segments.rds")
tmb_data <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/TMB_data.rds")

tmb_data <- tmb_data %>% 
  left_join(Segments %>% select(sample, best_K) %>% distinct(), by = join_by(sample))

library(tidyverse)

# Step 1: compute rank on a clean copy
tmb_plot_data <- tmb_data %>%
  group_by(IntoGen_cancer_type, best_K) %>%
  mutate(sample_rank = dplyr::row_number(TMB)) %>%
  ungroup()

# Step 2: fixed width + ordered cancer levels by median n_snvs
cancer_levels <- tmb_plot_data %>%
  group_by(IntoGen_cancer_type) %>%
  summarise(median_snvs = median(n_snvs), .groups = "drop") %>%
  arrange(median_snvs) %>%
  pull(IntoGen_cancer_type)

FIXED_WIDTH <- tmb_plot_data %>%
  group_by(IntoGen_cancer_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  pull(n) %>%
  max()

# Step 3: offsets per cancer type in correct order
offsets <- tibble(IntoGen_cancer_type = cancer_levels) %>%
  mutate(
    max_rank = FIXED_WIDTH,
    offset   = cumsum(lag(max_rank + 2, default = 0))
  )

K_SPACING <- 0.1  # increase if you want more separation

# Step 4: join offsets and compute global x position with K separation
tmb_plot_data <- tmb_plot_data %>%
  left_join(offsets, by = "IntoGen_cancer_type") %>%
  group_by(IntoGen_cancer_type) %>%
  mutate(
    k_index  = as.numeric(factor(best_K)),          # numeric index per K
    k_center = mean(unique(k_index)),               # center groups
    k_shift  = (k_index - k_center) * K_SPACING     # horizontal shift
  ) %>%
  ungroup() %>%
  mutate(
    x_pos = sample_rank + offset + k_shift
  )

# Step 5: label positions in correct order
label_pos <- tmb_plot_data %>%
  group_by(IntoGen_cancer_type) %>%
  summarise(x_center = mean(x_pos), .groups = "drop") %>%
  arrange(x_center)

# Step 6: vline separators
vline_positions <- offsets %>%
  mutate(vline = offset + max_rank + 1) %>%
  pull(vline)
vline_positions <- vline_positions[-length(vline_positions)]

# Step 7: median segments
medians <- tmb_plot_data %>%
  group_by(IntoGen_cancer_type, best_K) %>%
  summarise(
    median_TMB = median(TMB),
    xmin       = min(x_pos) - 0.4,
    xmax       = max(x_pos) + 0.4,
    .groups    = "drop"
  )

# Step 8: plot
p <- ggplot(tmb_plot_data, aes(x = x_pos, y = TMB, color = factor(best_K), shape = class)) +
  geom_vline(xintercept = vline_positions, color = "grey60", linewidth = 0.4, linetype = "dashed") +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_segment(
    data = medians,
    aes(x = xmin, xend = xmax, y = median_TMB, yend = median_TMB, color = factor(best_K)),
    linewidth = 0.8, inherit.aes = FALSE
  ) +
  scale_x_continuous(
    breaks = label_pos$x_center,
    labels = label_pos$IntoGen_cancer_type
  ) +
  scale_color_manual(
    values = c("1" = "#4DAF4A", "2" = "#377EB8", "3" = "#FF7F00", "4" = "#984EA3", "5" = "#E41A1C"),
    name = "Best k"
  ) +
  scale_shape_manual(
    values = c("HM" = 17, "WGD" = 16, "Classic" = 15),
    name = "Class"
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = c("0.01", "0.1", "1", "10", "100")
  ) +
  annotation_logticks(sides = "l") +
  labs(
    title = "Tumor Mutational Burden across Cancer Types",
    subtitle = paste0("n = ", nrow(tmb_plot_data), " samples | WGS (~", GENOME_SIZE_MB, " Mb)"),
    x = "Cancer Type",
    y = "TMB (mutations / Mb)",
    caption = "Ordered by median n_snvs | Lines = median per k group"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
    legend.position    = "bottom",
    plot.title         = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90")
  )

p

ggsave("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Plot/TMB_plot_bestk_ranked.pdf", p, width = 22, height = 7)





# =========================
# CNA + FGA (FIXED)
# =========================

# 1. correct tau_rank
Segments <- Segments %>%
  group_by(sample) %>%
  mutate(
    tau_rank = match(clock_mean, sort(unique(clock_mean)))
  ) %>%
  ungroup()

# 2. FIXED aggregation (no double counting)
sample_summary <- Segments %>%
  group_by(sample, IntoGen_cancer_type, tau_rank) %>%
  summarise(
    frac_genome_tau = first(frac_genome_affected_per_timing_group),  # 🔑 key fix
    n_cna_tau       = n(),
    .groups = "drop"
  ) %>%
  group_by(sample, IntoGen_cancer_type) %>%
  summarise(
    frac_genome_affected = sum(frac_genome_tau),
    n_cna_per_sample     = sum(n_cna_tau),
    .groups = "drop"
  )

# 3. same ordering + x positions
sample_summary <- sample_summary %>%
  mutate(IntoGen_cancer_type = factor(IntoGen_cancer_type, levels = cancer_levels)) %>%
  group_by(IntoGen_cancer_type) %>%
  mutate(sample_rank = row_number(n_cna_per_sample)) %>%
  ungroup() %>%
  left_join(offsets, by = "IntoGen_cancer_type") %>%
  mutate(
    x_pos = sample_rank + offset
  )

# =========================
# PLOTS
# =========================

p_cna <- ggplot(sample_summary, aes(x = x_pos, y = n_cna_per_sample)) +
  geom_vline(xintercept = vline_positions, color = "grey60", linewidth = 0.4, linetype = "dashed") +
  geom_point(alpha = 0.5, size = 1.2) +
  # scale_y_log10() +
  scale_x_continuous(breaks = label_pos$x_center, labels = label_pos$IntoGen_cancer_type) +
  labs(y = "Number of CNAs", x = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

p_frac <- ggplot(sample_summary, aes(x = x_pos, y = frac_genome_affected)) +
  geom_vline(xintercept = vline_positions, color = "grey60", linewidth = 0.4, linetype = "dashed") +
  geom_point(alpha = 0.5, size = 1.2, color = "steelblue") +
  scale_x_continuous(breaks = label_pos$x_center, labels = label_pos$IntoGen_cancer_type) +
  labs(y = "Fraction genome affected", x = "Cancer type") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

# =========================
# FINAL STACKED FIGURE
# =========================

p_final <- p / p_cna / p_frac + plot_layout(heights = c(1.2, 1, 1))

p_final

ggsave(
  "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Plot/TMB_CNA_FGA_combined.pdf",
  p_final,
  width = 22,
  height = 12
)


# verify fractionGA plot
Segments %>%
  group_by(sample, clock_mean) %>%
  summarise(n_unique = n_distinct(frac_genome_affected_per_timing_group)) %>%
  filter(n_unique > 1)
