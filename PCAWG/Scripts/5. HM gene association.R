
rm(list = ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(broom)
library(purrr)

dir.create("../Analysis_results/HM_gene_associtation", recursive = T)
res.dir = "../Analysis_results/HM_gene_associtation"

# --- Data loading -------------------------------------------------------------
df <- readRDS("../Data/Drivers.rds") %>%
  dplyr::mutate(mutation_status = ifelse(mutation_status == "WT", "WT", "M")) %>%
  dplyr::select(gene, mutation_status, class, ttype, sample_id)

CLASS_COLORS <- c(Classic = "#4393C3", HM = "#D6604D", WGD = "#74C476")

# One row per sample × gene
df <- df %>%
  distinct(sample_id, gene, .keep_all = TRUE) %>%
  mutate(
    mutated = mutation_status == "M",
    class   = factor(class, levels = c("Classic", "HM", "WGD")),
    ttype   = factor(ttype)
  )


# =============================================================================
# PART 1 — Raw mutation frequencies (for heatmap tiles)
# =============================================================================

# Complete sample x gene matrix: missing rows = WT (not missing data)
# Every sample must appear for every gene before computing frequencies
all_samples <- df %>% distinct(sample_id, class, ttype)
all_genes   <- df %>% distinct(gene)

df_complete <- all_samples %>%
  crossing(all_genes) %>%
  left_join(
    df %>% select(sample_id, gene, mutated),
    by = c("sample_id", "gene")
  ) %>%
  mutate(mutated = replace_na(mutated, FALSE))

# Recompute frequencies on the complete matrix
samples_per_class <- df_complete %>%
  group_by(class) %>%
  summarise(total_samples = n_distinct(sample_id), .groups = "drop")

gene_freq <- df_complete %>%
  filter(mutated) %>%
  group_by(class, gene) %>%
  summarise(n_samples = n_distinct(sample_id), .groups = "drop") %>%
  right_join(
    samples_per_class %>% crossing(all_genes),
    by = c("class", "gene")
  ) %>%
  mutate(
    n_samples = replace_na(n_samples, 0L),
    frequency = n_samples / total_samples
  )

# Keep genes where at least one class exceeds the frequency threshold
min_freq <- 0.05

genes_keep <- gene_freq %>%
  group_by(gene) %>%
  filter(any(frequency > min_freq)) %>%
  pull(gene) %>% unique()

gene_freq_filtered <- gene_freq %>%
  filter(gene %in% genes_keep)


# =============================================================================
# PART 2 — Tumor-type-adjusted significance (three pairwise logistic regressions)
#
# Rather than one model with a single reference class, we run three explicit
# pairwise comparisons, each on just the two relevant classes:
#   Classic vs HM  → star on Classic tile and HM tile
#   Classic vs WGD → star on Classic tile and WGD tile
#   HM vs WGD      → star on HM tile and WGD tile
#
# Each tile therefore gets a star from its OWN direct comparison:
#   Classic tile → Classic vs HM  and  Classic vs WGD
#   HM tile      → Classic vs HM  and  HM vs WGD
#   WGD tile     → Classic vs WGD and  HM vs WGD
# The most significant of the relevant comparisons is shown on each tile.
# =============================================================================

run_pairwise_logistic <- function(gene_name, data, class_a, class_b) {
  d <- data %>%
    filter(gene == gene_name, class %in% c(class_a, class_b)) %>%
    mutate(class = factor(class, levels = c(class_a, class_b)))
  
  if (n_distinct(d$mutated) < 2) return(NULL)
  if (any(table(d$class) == 0))  return(NULL)
  
  # Drop ttype levels not present in this subset
  d <- d %>% mutate(ttype = droplevels(ttype))
  
  fit  <- tryCatch(
    suppressWarnings(glm(mutated ~ class + ttype, data = d, family = binomial)),
    error = function(e) NULL
  )
  fit0 <- tryCatch(
    suppressWarnings(glm(mutated ~ ttype,          data = d, family = binomial)),
    error = function(e) NULL
  )
  if (is.null(fit) | is.null(fit0)) return(NULL)
  
  lrt_p <- tryCatch(
    anova(fit0, fit, test = "LRT")$`Pr(>Chi)`[2],
    error = function(e) NA_real_
  )
  
  coef_row <- tryCatch({
    tidy(fit, exponentiate = TRUE, conf.int = FALSE) %>%
      filter(startsWith(term, "class")) %>%
      slice(1) %>%
      transmute(
        gene       = gene_name,
        class_a    = class_a,
        class_b    = class_b,
        comparison = paste(class_a, "vs", class_b),
        OR         = estimate,          # odds of class_b vs class_a
        CI_low     = exp(log(estimate) - 1.96 * std.error),
        CI_high    = exp(log(estimate) + 1.96 * std.error),
        p          = lrt_p              # LRT p (tumor-type adjusted)
      )
  }, error = function(e) NULL)
  
  coef_row
}

cat("Running three pairwise tumor-type-adjusted logistic regressions...\n")

pairs <- list(
  c("Classic", "HM"),
  c("Classic", "WGD"),
  c("HM",      "WGD")
)

pairwise_res <- map_dfr(pairs, function(p) {
  map_dfr(genes_keep, ~run_pairwise_logistic(.x, df_complete, p[1], p[2]))
})

# FDR correction within each comparison independently
pairwise_res <- pairwise_res %>%
  group_by(comparison) %>%
  mutate(q = p.adjust(p, method = "BH")) %>%
  ungroup()

# Summary per comparison
for (comp in unique(pairwise_res$comparison)) {
  n_sig <- sum(pairwise_res$q[pairwise_res$comparison == comp] < 0.05,
               na.rm = TRUE)
  cat(sprintf("  %s: %d genes q < 0.05\n", comp, n_sig))
}


# =============================================================================
# PART 3 — Significance stars per tile
#
# Each tile shows the star from the most relevant direct comparison:
#   Classic tile → min q of (Classic vs HM, Classic vs WGD)
#   HM tile      → min q of (Classic vs HM, HM vs WGD)
#   WGD tile     → min q of (Classic vs WGD, HM vs WGD)
# =============================================================================

make_stars <- function(q) {
  case_when(
    q < 0.001 ~ "***",
    q < 0.01  ~ "**",
    q < 0.05  ~ "*",
    q < 0.10  ~ "·",
    TRUE      ~ ""
  )
}

# For each tile, collect the q-values from comparisons that involve that class
# and take the minimum (most significant direct comparison)
tile_stars_full <- bind_rows(
  # Classic tile: involved in Classic vs HM and Classic vs WGD
  pairwise_res %>%
    filter(comparison %in% c("Classic vs HM", "Classic vs WGD")) %>%
    group_by(gene) %>%
    slice_min(q, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(gene, class = "Classic", q, OR),
  
  # HM tile: involved in Classic vs HM and HM vs WGD
  pairwise_res %>%
    filter(comparison %in% c("Classic vs HM", "HM vs WGD")) %>%
    group_by(gene) %>%
    slice_min(q, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(gene, class = "HM", q,
              # flip OR so it always represents "this class vs the other"
              OR = ifelse(comparison == "Classic vs HM", 1 / OR, OR)),
  
  # WGD tile: involved in Classic vs WGD and HM vs WGD
  pairwise_res %>%
    filter(comparison %in% c("Classic vs WGD", "HM vs WGD")) %>%
    group_by(gene) %>%
    slice_min(q, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(gene, class = "WGD", q, OR)
) %>%
  mutate(
    star  = make_stars(q),
    class = factor(class, levels = c("Classic", "HM", "WGD"))
  )


# =============================================================================
# PART 4 — Gene ordering by hierarchical clustering
# =============================================================================

freq_wide <- gene_freq_filtered %>%
  select(gene, class, frequency) %>%
  pivot_wider(names_from = class, values_from = frequency, values_fill = 0) %>%
  tibble::column_to_rownames("gene")

gene_levels <- rownames(freq_wide)[hclust(dist(freq_wide))$order]

gene_freq_filtered <- gene_freq_filtered %>%
  mutate(
    gene  = factor(gene,  levels = gene_levels),
    class = factor(class, levels = c("Classic", "HM", "WGD"))
  )

# Mark dominant class (highest frequency per gene) for border overlay
dominant_class <- gene_freq_filtered %>%
  group_by(gene) %>%
  slice_max(frequency, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(gene, dominant_class = class)

# Final plotting data: merge frequencies + stars
plot_data <- gene_freq_filtered %>%
  left_join(dominant_class, by = "gene") %>%
  mutate(is_dominant = class == dominant_class) %>%
  left_join(tile_stars_full, by = c("gene", "class")) %>%
  mutate(
    tile_label = case_when(
      frequency > 0 & !is.na(star) & star != "" ~
        paste0(percent(frequency, accuracy = 1), "\n", star),
      frequency > 0 ~
        percent(frequency, accuracy = 1),
      TRUE ~ ""
    )
  )


# =============================================================================
# PART 5 — PLOTS
# =============================================================================

heatmap_theme <- function() {
  theme_minimal(base_size = 11) +
    theme(
      axis.text.x     = element_text(face = "bold",
                                     color = CLASS_COLORS[c("Classic","HM","WGD")]),
      axis.text.y     = element_text(size = 8),
      panel.grid      = element_blank(),
      legend.position = "right",
      plot.title      = element_text(face = "bold"),
      plot.subtitle   = element_text(color = "gray40", size = 9)
    )
}

star_caption <- "Stars = tumor-type-adjusted FDR  ·q<0.10  *q<0.05  **q<0.01  ***q<0.001"

# --- Plot 1: Absolute frequency -----------------------------------------------

p1 <- ggplot(plot_data, aes(x = class, y = gene, fill = frequency)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_tile(
    data  = filter(plot_data, is_dominant & frequency > 0),
    color = "gray20", linewidth = 0.9, fill = NA
  ) +
  geom_text(aes(label = tile_label), size = 2.5, color = "gray20",
            lineheight = 0.85) +
  scale_fill_gradient(
    low = "white", high = "steelblue",
    labels = percent, name = "Mutation\nfrequency"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    title    = "Driver mutation frequency per evolutionary class",
    subtitle = paste0("Genes mutated in >", min_freq * 100,
                      "% of samples in ≥1 class  |  Border = dominant class"),
    caption  = star_caption,
    x = NULL, y = NULL
  ) +
  heatmap_theme()
p1

# --- Plot 2: Relative enrichment (deviation from gene mean) -------------------

plot_data_rel <- plot_data %>%
  group_by(gene) %>%
  mutate(
    mean_freq     = mean(frequency),
    relative_freq = frequency - mean_freq
  ) %>%
  ungroup()

max_dev <- max(abs(plot_data_rel$relative_freq), na.rm = TRUE)

p2 <- ggplot(plot_data_rel, aes(x = class, y = gene, fill = relative_freq)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = tile_label), size = 2.5, color = "gray20",
            lineheight = 0.85) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC",
    midpoint = 0, limits = c(-max_dev, max_dev),
    labels = percent, name = "Frequency\nvs gene mean"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    title    = "Driver enrichment relative to gene mean",
    subtitle = "Blue = enriched vs average  |  Red = depleted  |  Labels = absolute frequency",
    caption  = star_caption,
    x = NULL, y = NULL
  ) +
  heatmap_theme()

# --- Save ---------------------------------------------------------------------
h <- max(4, length(gene_levels) * 0.28 + 2)
ggsave(file.path(res.dir, "driver_heatmap_absolute.pdf"), p1, width = 6.5, height = h)
ggsave(file.path(res.dir, "driver_heatmap_relative.pdf"), p2, width = 6.5, height = h)

print(p1)
print(p2)


# =============================================================================
# PART 6 — Volcano plots (one panel per pairwise comparison)
# =============================================================================
# x-axis : log2(OR) — effect size and direction
#           positive = more mutated in class_b, negative = more in class_a
# y-axis : -log10(q) — significance
# color  : which class the gene is enriched in (or non-significant)
# labels : genes passing significance + effect size threshold
# =============================================================================

# Comparison-specific color palettes:
# each comparison gets two colors (one per class) + gray for non-significant
volcano_colors <- list(
  "Classic vs HM"  = c(Classic = "#4393C3", HM      = "#D6604D", NS = "grey75"),
  "Classic vs WGD" = c(Classic = "#4393C3", WGD     = "#74C476", NS = "grey75"),
  "HM vs WGD"      = c(HM      = "#D6604D", WGD     = "#74C476", NS = "grey75")
)

q_thresh    <- 0.05   # significance threshold
or_thresh   <- log2(1.5)  # minimum effect size to label (log2 OR)
top_n_label <- 15     # max genes to label per direction

volcano_data <- pairwise_res %>%
  filter(!is.na(OR), !is.na(q)) %>%
  mutate(
    log2OR   = log2(OR),
    neg_logq = -log10(q),
    # Enrichment direction: which class has higher mutation rate?
    enriched_in = case_when(
      q < q_thresh & OR > 1 ~ class_b,
      q < q_thresh & OR < 1 ~ class_a,
      TRUE                  ~ "NS"
    )
  )

# Select genes to label: significant + large effect, capped at top_n per side
label_genes <- volcano_data %>%
  filter(enriched_in != "NS", abs(log2OR) >= or_thresh) %>%
  group_by(comparison, enriched_in) %>%
  slice_min(q, n = top_n_label, with_ties = FALSE) %>%
  ungroup() %>%
  select(gene, comparison)

volcano_data <- volcano_data %>%
  left_join(label_genes %>% mutate(do_label = TRUE),
            by = c("gene", "comparison")) %>%
  mutate(
    do_label = replace_na(do_label, FALSE),
    label    = ifelse(do_label, as.character(gene), NA_character_)
  )

# Build one plot per comparison, then combine with patchwork
library(ggrepel)
library(patchwork)

make_volcano <- function(comp, data, colors) {
  d      <- data %>% filter(comparison == comp)
  cols   <- colors[[comp]]
  cls    <- names(cols)     # e.g. c("Classic", "HM", "NS")
  class_a <- unique(d$class_a)
  class_b <- unique(d$class_b)
  
  # Count sig genes per side for subtitle
  n_a <- sum(d$enriched_in == class_a, na.rm = TRUE)
  n_b <- sum(d$enriched_in == class_b, na.rm = TRUE)
  
  ggplot(d, aes(x = log2OR, y = neg_logq, color = enriched_in)) +
    # Non-significant points first (behind)
    geom_point(data = filter(d, enriched_in == "NS"),
               size = 1.5, alpha = 0.4) +
    # Significant points on top
    geom_point(data = filter(d, enriched_in != "NS"),
               size = 2.2, alpha = 0.85) +
    # Gene labels (only significant + large effect)
    geom_text_repel(
      aes(label = label),
      size          = 2.8,
      max.overlaps  = 20,
      min.segment.length = 0.2,
      segment.color = "grey50",
      segment.size  = 0.3,
      show.legend   = FALSE
    ) +
    # Threshold lines
    geom_hline(yintercept = -log10(q_thresh),
               linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = c(-or_thresh, or_thresh),
               linetype = "dotted", color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey30", linewidth = 0.3) +
    scale_color_manual(
      values = cols,
      labels = c(class_a, class_b, "NS"),
      name   = "Enriched in"
    ) +
    scale_x_continuous(
      labels = function(x) ifelse(x > 0,
                                  paste0("+", round(x, 1)),
                                  round(x, 1))
    ) +
    annotate("text", x = -Inf, y = Inf,
             label = paste0(n_a, " genes"), hjust = -0.2, vjust = 1.5,
             color = cols[class_a], size = 3, fontface = "bold") +
    annotate("text", x =  Inf, y = Inf,
             label = paste0(n_b, " genes"), hjust = 1.2, vjust = 1.5,
             color = cols[class_b], size = 3, fontface = "bold") +
    labs(
      title = comp,
      x     = paste0("log2 OR  (← more in ", class_a,
                     "  |  more in ", class_b, " →)"),
      y     = "-log10(q)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, color = "grey92"),
      plot.title       = element_text(face = "bold", size = 11),
      legend.position  = "bottom",
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 9)
    )
}

p_volc_list <- map(
  names(volcano_colors),
  ~make_volcano(.x, volcano_data, volcano_colors)
)

# Combine three panels side by side
p_volcano <- wrap_plots(p_volc_list, nrow = 1) +
  plot_annotation(
    title    = "Tumor-type-adjusted driver enrichment across evolutionary classes",
    subtitle = paste0("Dashed line: q = ", q_thresh,
                      "  |  Dotted lines: |log2 OR| = ", round(or_thresh, 2),
                      "  |  Labeled: top ", top_n_label, " significant genes per side"),
    theme = theme(plot.title    = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(color = "grey40", size = 9))
  )
p_volcano

ggsave(file.path(res.dir, "driver_volcano.pdf"), p_volcano, width = 15, height = 6)
print(p_volcano)
cat("  driver_volcano.pdf\n")

# --- Summary table ----
summary_table <- gene_freq_filtered %>%
  select(gene, class, frequency) %>%
  pivot_wider(names_from = class, values_from = frequency,
              names_prefix = "freq_") %>%
  left_join(
    pairwise_res %>%
      select(gene, comparison, OR, q) %>%
      pivot_wider(names_from  = comparison,
                  values_from = c(OR, q),
                  names_sep   = "_"),
    by = "gene"
  ) %>%
  arrange(`q_Classic vs HM`)

write.table(summary_table, file.path(res.dir, "driver_signature_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nGenes retained:", length(gene_levels), "\n")
cat("Threshold: any class >", percent(min_freq), "\n")
cat("Outputs: driver_heatmap_absolute.pdf\n")
cat("         driver_heatmap_relative.pdf\n")
cat("         driver_signature_summary.tsv\n")
