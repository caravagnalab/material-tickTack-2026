.libPaths("~/R/orfeo_R_4.4/")
library(dplyr)
base_path <- "/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack_2026/simulations/"
data_dir = "/orfeo/cephfs/scratch/cdslab/scocomello/material_tickTack_2026/simulations/data/"

selected_scenarios_expanded <- readRDS(paste0(data_dir,"selected_scenarios.rds"))

selected_scenarios_expanded %>%
  select(
    clusters     = n_clusters,
    segments     = n_seg,
    purity       = purity_sim,
    coverage     = coverage_sim,
    density      = median_density_mut
  ) %>%
  write.table(
    file      = paste0(base_path, "parameter_combinations.txt"),
    sep       = ",",
    row.names = FALSE,
    col.names = FALSE,
    quote     = FALSE
  )
