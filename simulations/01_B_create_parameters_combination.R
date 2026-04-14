base_path = "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/simulations/"

# Define parameter values
clusters <- c(1, 2, 3, 4, 5)
segments <- c(5, 20, 40)
purity <- c(0.4, 0.8)
coverage <- c(40, 80)
mutation_density <- c(2e-07,1e-06,2e-06,3e-06)

# Create full grid of combinations
grid <- expand.grid(
  clusters = clusters,
  segments = segments,
  purity = purity,
  coverage = coverage,
  density = mutation_density
)

# Reorder columns to match the input order in the next script 
grid <- grid[, c("clusters", "segments", "purity", "coverage", "density")]


# Write to comma-separated txt file
write.table(
  grid,
  file = paste0(base_path,"parameter_combinations.txt"),
  sep = ",",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
