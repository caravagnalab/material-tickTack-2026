library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

Drivers <- readRDS("~/Docs/GitHub/material-tickTack-2026/PCAWG/Data/Drivers.rds") %>% 
  filter(mutation_status != "WT") %>%
  select(segment_id, gene, 
         karyotype, sample_id,
         mutation_status, mult_estimate,
         best_K, IntoGen_cancer_type, class) 

# ── Data preparation ───────────────────────────────────────────────────────────

# Remove rows where best_K is NA (these samples can't be assigned to a group),
# then create a binary grouping variable: K1 for samples with best_K == 1,
# K_gt1 for all others (best_K 2–5). This pools the higher-K samples together
# to increase statistical power.
df_filt <- Drivers |>
  filter(!is.na(best_K)) |>
  mutate(K_group = ifelse(best_K >= 2, "K1", "K_gt1"))


# ── Fisher test function ───────────────────────────────────────────────────────

# This function receives a subset of the data (one cancer_type × class stratum)
# and a gene name, and returns a one-row tibble with the association statistics.
run_fisher <- function(data, gene_name) {
  
  # For each sample in this stratum, determine whether the gene of interest
  # is present (i.e. appears in at least one row for that sample).
  # Using distinct() on sample_id + K_group + gene_present ensures that
  # samples appearing multiple times (duplicate rows) are counted only once.
  samples <- data |>
    mutate(gene_present = gene == gene_name) |>
    distinct(sample_id, K_group, gene_present) |>
    group_by(sample_id, K_group) |>
    # Collapse to one row per sample: TRUE if the gene was present in any row
    summarise(gene_present = any(gene_present), .groups = "drop")
  
  # Build the 2×2 contingency table:
  #   rows    → K_group (K1 vs K_gt1)
  #   columns → gene_present (TRUE vs FALSE)
  tab <- table(samples$K_group, samples$gene_present)
  
  # If one of the groups is entirely absent (e.g. no K1 samples in this stratum,
  # or the gene never appears), the table won't be 2×2 and Fisher's test would
  # fail. Return NULL to skip this gene/stratum combination safely.
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)
  
  # Run Fisher's exact test. simulate.p.value = TRUE uses Monte Carlo simulation
  # to handle tables with zero or very small cell counts, which is common with
  # small sample sizes. This avoids p-value inflation from the normal approximation.
  ft <- fisher.test(tab, simulate.p.value = TRUE)
  
  # Return a tidy one-row summary for this gene in this stratum.
  # Counts of gene-positive samples in each group are included to help interpret
  # the OR and flag cases where significance is driven by very few events.
  tibble(
    gene        = gene_name,
    n_K1        = sum(samples$K_group == "K1"),          # total samples in K1 group
    n_Kgt1      = sum(samples$K_group == "K_gt1"),       # total samples in K_gt1 group
    n_gene_K1   = sum(samples$K_group == "K1"   & samples$gene_present), # gene-positive in K1
    n_gene_Kgt1 = sum(samples$K_group == "K_gt1" & samples$gene_present), # gene-positive in K_gt1
    OR          = ft$estimate,   # odds ratio (how much more likely the gene is in one group)
    p_value     = ft$p.value
  )
}


# ── Run tests across all strata ────────────────────────────────────────────────

# Split the filtered data by cancer_type × class, then within each stratum
# loop over all genes present and apply run_fisher().
# group_modify() passes each sub-dataframe as .x to the anonymous function.
# map_dfr() iterates over the gene names and row-binds the results into one tibble.
results <- df_filt |>
  group_by(IntoGen_cancer_type, class) |>
  group_modify(~ {
    genes_in_stratum <- unique(.x$gene)               # genes observed in this stratum
    map_dfr(genes_in_stratum, \(g) run_fisher(.x, g)) # test each gene, bind rows
  }) |>
  ungroup()


# ── Multiple testing correction ────────────────────────────────────────────────

# Apply Benjamini-Hochberg FDR correction separately within each stratum.
# Correcting within rather than globally is appropriate here because each
# cancer_type × class combination is an independent analysis with its own
# set of hypotheses (genes tested). A global correction would be too conservative.
# Genes with q_value < 0.05 can be considered statistically significant
# after accounting for multiple comparisons.
results <- results |>
  group_by(IntoGen_cancer_type, class) |>
  mutate(q_value = p.adjust(p_value, method = "BH")) |>
  ungroup() |>
  arrange(IntoGen_cancer_type, class, p_value) # sort by stratum then by evidence strength

results %>% filter(p_value <= 0.1)

