library(dplyr)
library(tidyverse)
library(ggpubr)
# Segments <- readRDS("~/Docs/GitHub/material-tickTack-2026/PCAWG/Data/Segments.rds")
# Samples <- readRDS("~/Docs/GitHub/material-tickTack-2026/PCAWG/Data/Samples.rds")
Segments <- readRDS("../Data/Segments.rds")
Samples <- readRDS("../Data/Samples.rds")


counts <- Samples %>% filter(!is.na(IntoGen_cancer_type)) %>%
  count(IntoGen_cancer_type, class) %>%
  group_by(IntoGen_cancer_type) %>%
  mutate(total = sum(n)) %>%
  ungroup()

class_order <- counts %>%
  distinct(IntoGen_cancer_type, total) %>%
  arrange(desc(total)) %>%
  pull(IntoGen_cancer_type)

counts$IntoGen_cancer_type <- factor(counts$IntoGen_cancer_type,
                       levels = class_order)

ggplot(counts,
       aes(IntoGen_cancer_type, n)) +
  geom_col(aes(fill = class, color = class), alpha = 0.5) +
  geom_text(aes(label = n, color = class),
            position = position_stack(vjust = 0.5), ) +
  theme_bw() + xlab("IntoGen cancer type") + ylab("")+
  scale_fill_manual(
    values = alpha(c("Classic" = "darkgrey", "HM" = "steelblue", "WGD" = "firebrick"), 0.2),
    name = "Class"
  ) +
  scale_color_manual(
    values = alpha(c("Classic" = "darkgrey", "HM" = "steelblue", "WGD" = "firebrick"), 1),
    name = "Class"
  ) + theme(axis.text.x = element_text(angle=90, hjust=1))


# make donuts
# --- Data (replace with real counts) ---
data_all <- tibble(
  class = factor(c("Classic", "WGD", "HM"), levels = c("Classic", "WGD", "HM")),
  n     = c(Samples %>% filter(class == "Classic") %>% nrow(),
            Samples %>% filter(class == "WGD") %>% nrow(),
            Samples %>% filter(class == "HM") %>% nrow())
)

# --- Colors ---
class_colors <- c("Classic" = "#888780", "WGD" = "#B03A2E", "HM" = "#2E86C1")

# --- Donut function ---
make_donut <- function(df, title, subtitle) {
  df <- df %>%
    mutate(
      pct     = n / sum(n),
      pct_lab = paste0(round(pct * 100, 1), "%"),
      ymax    = cumsum(pct),
      ymin    = lag(ymax, default = 0),
      ymid    = (ymin + ymax) / 2,
      x_lab   = 1.15 * cos(2 * pi * (1 - ymid)),
      y_lab   = 1.15 * sin(2 * pi * (1 - ymid))
    )
  
  ggplot(df) +
    geom_rect(
      aes(xmin = 3, xmax = 4.2, ymin = ymin, ymax = ymax, fill = class),
      color = "white", linewidth = 0.6
    ) +
    geom_text(
      aes(x = 3.6, y = ymid, label = pct_lab),
      size = 3.2, fontface = "bold", color = "white"
    ) +
    annotate("text", x = 0, y = 0, label = subtitle,
             size = 3, color = "grey50", hjust = 0.5, vjust = 0.5) +
    scale_fill_manual(values = class_colors, name = "Class") +
    coord_polar(theta = "y", start = 0, direction = -1) +
    xlim(0, 4.5) +
    labs(title = title) +
    theme_void(base_size = 12) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 13,
                                   margin = margin(b = 4)),
      legend.position = "none",
      plot.margin   = margin(10, 20, 10, 20)
    )
}

# --- Build plots ---
n_tot = data_all$n %>% sum()
p_all <- make_donut(data_all, "All tumor types", paste0("n = ", n_tot))
# p_brca <- make_donut(data_brca, "BRCA",          "n = 207")

# --- Shared legend ---
p_legend <- ggplot(data_all, aes(x = 1, y = n, fill = class)) +
  geom_col() +
  scale_fill_manual(
    values = class_colors,
    name   = "Class",
    guide  = guide_legend(
      keywidth  = unit(0.35, "cm"),
      keyheight = unit(0.35, "cm"),
      label.theme = element_text(size = 10)
    )
  ) +
  theme_void() +
  theme(legend.position = "right")

legend_only <- cowplot::get_legend(p_legend)

# --- Combine ---
# final <- (p_all | p_brca | patchwork::wrap_elements(legend_only)) +
final <- (p_all | patchwork::wrap_elements(legend_only)) +
  plot_layout(widths = c(1, 1, 0.35))

final



counts %>% filter(total > 10) %>%
  mutate(f = n / total) %>% filter(class == "HM") %>%
  arrange(desc(f))

counts %>% filter(total > 10) %>%
  mutate(f = n / total) %>% filter(class == "WGD") %>%
  arrange(desc(f))

df1 = Segments %>% 
  filter(!is.na(IntoGen_cancer_type)) %>% 
  group_by(sample, clock_mean, class, IntoGen_cancer_type, HM_cluster) %>%
  summarise(n_segments = n()) %>% ungroup() %>%
  group_by(sample) %>% mutate(cluster_per_sample = n()) %>%
  ungroup() %>% group_by(sample) %>%
  mutate(
    to_keep = ifelse(
      class %in% c("WGD", "HM") & cluster_per_sample > 1 & n_segments != max(n_segments), F, T
    )
  ) %>% filter(to_keep)

## Compare timing significance - WGS vs Classic

library(tidyverse)
library(ggpubr)      # for stat_compare_means and significance bars
library(rstatix)     # for tidy statistical tests and filtering by group size

# ── 1. Filter to the two classes of interest ──────────────────────────────────
# Keep only "Classic" and "WGD" (and "HM" where present)
df_filt <- df1 %>%
  filter(to_keep == TRUE) %>%                          # optional: respect your existing flag
  filter(class %in% c("Classic", "WGD", "HM")) # "HM"

# ── 2. Drop tumour-type × class combinations with too few samples ─────────────
# Threshold: require at least N samples per group. Adjust MIN_N as needed.
MIN_N <- 5

group_counts <- df_filt %>%
  group_by(IntoGen_cancer_type, class) %>%
  summarise(n = n(), .groups = "drop")

# Keep only cancer types where EVERY class present has >= MIN_N samples
# (so we don't run a test with one arm having 2 observations)
valid_cancer_types <- group_counts %>%
  filter(n >= MIN_N) %>%
  group_by(IntoGen_cancer_type) %>%
  filter(n_distinct(class) >= 2) %>%   # need at least 2 classes to compare
  pull(IntoGen_cancer_type) %>%
  unique()

df_plot <- df_filt %>%
  filter(IntoGen_cancer_type %in% valid_cancer_types)

# ── 3. Build the comparison pairs for each cancer type ────────────────────────
# For tumour types with HM, we do all three pairwise comparisons;
# for the rest, just Classic vs WGD.
make_comparisons <- function(cancer_type, data) {
  classes <- data %>%
    filter(IntoGen_cancer_type == cancer_type) %>%
    pull(class) %>%
    unique() %>%
    sort()
  combn(classes, 2, simplify = FALSE)   # all pairwise pairs
}

comparisons_by_type <- setNames(
  lapply(valid_cancer_types, make_comparisons, data = df_plot),
  valid_cancer_types
)

# ── 4. Plot ───────────────────────────────────────────────────────────────────
# We facet by cancer type. Because each facet may have different comparison
# pairs, we use ggpubr::stat_compare_means() which handles this gracefully.

# Colour palette — tweak to taste
class_colours <- c(
  "Classic" = "darkgrey",
  "WGD"     = "firebrick",
  "HM"      = "steelblue"
)

# p-value label style: shows stars + exact p-value. Change to "p.signif"
# for stars only, or "p.format" for the number only.
LABEL_TYPE <- "p.signif"   # "****", "***", "**", "*", "ns"

p <- ggplot(df_plot, aes(x = class, y = clock_mean, fill = class)) +
  
  # ── Boxplot layer ──────────────────────────────────────────────────────────
  geom_boxplot(
    outlier.shape  = 21,
    outlier.size   = 1.5,
    outlier.alpha  = 0.5,
    width          = 0.55,
    colour         = "grey30"
  ) +
  
  # ── Individual points (optional but useful for small groups) ───────────────
  geom_jitter(
    width = 0.15, size = 0.8, alpha = 0.35, colour = "grey20"
  ) +
  
  # ── Significance bars ──────────────────────────────────────────────────────
  # stat_compare_means uses a Wilcoxon rank-sum test by default
  # (non-parametric equivalent of t-test; safer given unknown distributions).
  # Switch to method = "t.test" if you prefer a parametric test.
  stat_compare_means(
    comparisons = comparisons_by_type[[1]],  # placeholder — overridden by facet
    method      = "wilcox.test",
    label       = LABEL_TYPE,
    hide.ns     = TRUE,   # set TRUE to suppress "ns" bars if you prefer cleaner plots
    step.increase = 0.08   # vertical spacing between stacked bars
  ) +
  
  # ── Facet: one panel per cancer type ──────────────────────────────────────
  facet_wrap(
    ~ IntoGen_cancer_type,
    scales = "free_x",    # each facet shows only its own classes
    nrow   = NULL         # let ggplot choose the grid layout
  ) +
  
  # ── Aesthetics ─────────────────────────────────────────────────────────────
  scale_fill_manual(values = class_colours) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +  # headroom for bars
  labs(
    x     = NULL,
    y     = "Clock mean",
    fill  = "Class",
    title = "Clock mean by class across tumour types",
    subtitle = paste0(
      "Wilcoxon rank-sum test |",
      " samples per group | stars: * p<0.05, ** p<0.01, *** p<0.001, **** p<0.0001"
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92", colour = "grey70"),
    strip.text       = element_text(face = "bold", size = 9),
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank()
  )

# ── 5. Render / save ──────────────────────────────────────────────────────────
print(p)

# Save to file — adjust width/height to taste (more cancer types = wider plot)
ggsave(
  "clock_mean_by_class_and_tumour.pdf",
  plot   = p,
  width  = 14,
  height = 8,
  device = cairo_pdf   # handles unicode and special chars cleanly
)

#### Compare timing
df <- df1 %>%
  group_by(IntoGen_cancer_type) %>%
  mutate(median_clock = median(clock_mean)) %>%
  ungroup() %>%
  arrange((median_clock)) %>%
  mutate(IntoGen_cancer_type = factor(IntoGen_cancer_type, levels = unique(IntoGen_cancer_type)))

# Step 1: compute rank on a clean copy
df_plot_data <- df %>%
  group_by(IntoGen_cancer_type) %>%
  mutate(sample_rank = dplyr::row_number(clock_mean)) %>%
  ungroup()

# Step 2: fixed width + ordered cancer levels by median n_snvs
cancer_levels <- df_plot_data %>%
  group_by(IntoGen_cancer_type) %>%
  summarise(median_clock = median(clock_mean), .groups = "drop") %>%
  arrange(median_clock) %>%
  pull(IntoGen_cancer_type)

FIXED_WIDTH <- df_plot_data %>%
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
df_plot_data <- df_plot_data %>%
  left_join(offsets, by = "IntoGen_cancer_type") %>%
  group_by(IntoGen_cancer_type) %>%
  # mutate(
  #   k_index  = as.numeric(factor(best_K)),          # numeric index per K
  #   k_center = mean(unique(k_index)),               # center groups
  #   k_shift  = (k_index - k_center) * K_SPACING     # horizontal shift
  # ) %>%
  ungroup() %>%
  mutate(
    x_pos = sample_rank + offset #+ k_shift
  )

# Step 5: label positions in correct order
label_pos <- df_plot_data %>%
  group_by(IntoGen_cancer_type) %>%
  summarise(x_center = mean(x_pos), .groups = "drop") %>%
  arrange(x_center)

# Step 6: vline separators
vline_positions <- offsets %>%
  mutate(vline = offset + max_rank + 1) %>%
  pull(vline)
vline_positions <- vline_positions[-length(vline_positions)]

# Step 7: median segments
medians <- df_plot_data %>%
  group_by(IntoGen_cancer_type) %>%
  summarise(
    median_clock = median(clock_mean),
    xmin       = min(x_pos) - 0.4,
    xmax       = max(x_pos) + 0.4,
    .groups    = "drop"
  )

# Step 8: plot
p <- ggplot(df_plot_data, aes(x = x_pos, y = clock_mean, color = class)) +
  geom_vline(xintercept = vline_positions, color = "grey60", linewidth = 0.4, linetype = "dashed") +
  geom_point(size = 1.2, alpha = 0.5) +
  geom_segment(
    data = medians,
    aes(x = xmin, xend = xmax, y = median_clock, yend = median_clock),
    linewidth = 0.8, inherit.aes = FALSE
  ) +
  scale_x_continuous(
    breaks = label_pos$x_center,
    labels = label_pos$IntoGen_cancer_type
  ) +
  scale_color_manual(
    values = c("Classic" = "grey", "HM" = "steelblue", "WGD" = "firebrick"),
    name = "Class"
  ) +
  # scale_shape_manual(
  #   values = c("HM" = 17, "WGD" = 16, "Classic" = 15),
  #   name = "Class"
  # ) +
  # scale_y_log10(
  #   breaks = c(0.01, 0.1, 1, 10, 100),
  #   labels = c("0.01", "0.1", "1", "10", "100")
  # ) +
  # annotation_logticks(sides = "l") +
  labs(
    title = "Time of the copy number events across Cancer Types (fpr HM we are timing the HM cluster)",
    #subtitle = paste0("n = ", nrow(tmb_plot_data), " samples | WGS (~", GENOME_SIZE_MB, " Mb)"),
    x = "Cancer Type",
    y = "Median Clock",
    # caption = "Ordered by median n_snvs | Lines = median per k group"
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

# ggsave("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Plot/TMB_plot_bestk_ranked.pdf", p, width = 22, height = 7)



