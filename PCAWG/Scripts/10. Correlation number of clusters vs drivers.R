library(dplyr)
library(purrr)
library(tidyr)

Drivers <- readRDS("/orfeo/scratch/cdslab/antonelloa/material-tickTack-2026/PCAWG/Data/Drivers.rds") %>% 
  filter(mutation_status != "WT") %>%
  select(segment_id, gene, 
         karyotype, sample_id,
         mutation_status, mult_estimate,
         best_K, IntoGen_cancer_type, class) 

# Step 1: prep
df_filt <- Drivers %>%
  filter(!is.na(best_K)) |>
  mutate(K_group = ifelse(best_K == 1, "K1", "K_gt1"))

# Step 2-3: Fisher test per cancer_type × class × gene
run_fisher <- function(data, gene_name) {
  # One row per sample: was this gene altered?
  samples <- data %>%
    mutate(gene_present = gene == gene_name) %>%
    distinct(sample_id, K_group, gene_present) %>%
    # A sample "has" the gene if any row for it is TRUE
    group_by(sample_id, K_group) %>%
    summarise(gene_present = any(gene_present), .groups = "drop")
  
  tab <- table(samples$K_group, samples$gene_present)
  
  # Need at least 2x2 with both groups represented
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)
  
  ft <- fisher.test(tab, simulate.p.value = TRUE)
  
  tibble(
    gene        = gene_name,
    n_K1        = sum(samples$K_group == "K1"),
    n_Kgt1      = sum(samples$K_group == "K_gt1"),
    n_gene_K1   = sum(samples$K_group == "K1"   & samples$gene_present),
    n_gene_Kgt1 = sum(samples$K_group == "K_gt1" & samples$gene_present),
    OR          = ft$estimate,
    p_value     = ft$p.value
  )
}

results <- df_filt |>
  group_by(IntoGen_cancer_type, class) |>
  group_modify(~ {
    genes_in_stratum <- unique(.x$gene)
    map_dfr(genes_in_stratum, \(g) run_fisher(.x, g))
  }) |>
  ungroup()

# Step 4: FDR correction within each stratum
results <- results |>
  group_by(IntoGen_cancer_type, class) |>
  mutate(q_value = p.adjust(p_value, method = "BH")) |>
  ungroup() |>
  arrange(IntoGen_cancer_type, class, p_value)