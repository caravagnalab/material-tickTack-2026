rm(list=ls())
.libPaths("~/R/orfeo_R_4.4/")
library(dplyr)
library(ggplot2)
library(patchwork)

FIT_DIR <-"/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit/inference_results_5ncomponents/"
PLOT_DIR <-"/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit_preparation_for_analysis//plot/"
RES_FINAL_DIR <- "/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Fit_preparation_for_analysis/results/"

################################################################
################ check robustness with standardized data #######

info_HM <- readRDS(paste0(RES_FINAL_DIR, "00_HM_max_clusters.rds"))
Segments <- readRDS("/orfeo/cephfs/scratch/cdslab/scocomello/material-tickTack-2026/PCAWG/Data/Segments.rds")


info_HM_maxgroup = info_HM %>% 
  group_by(sample) %>% 
  filter(n_cna_events_per_timing_group == max(n_cna_events_per_timing_group)) %>%
  ungroup()

df <- info_HM_maxgroup 
# %>%
#   select(n_chr_affected_per_timing_group,
#           frac_genome_affected_per_timing_group,
#           n_cna_events_per_timing_group,
#           ploidy)


vars <- df %>%
  dplyr::select(n_chr_affected_per_timing_group,
         frac_genome_affected_per_timing_group,
         n_cna_events_per_timing_group,
         ploidy)

X_scaled <- scale(vars)

set.seed(123)

kmeans_fit <- kmeans(X_scaled, centers = 2, nstart = 50)

df$cluster <- kmeans_fit$cluster

#### BOOTSTRAP CLUSTER STABILITY - ROBUSTNESS ASSESMENT ###########

B <- 200
ari_values <- numeric(B)

for (b in 1:B) {
  
  idx <- sample(1:nrow(X_scaled), replace = TRUE)
  
  X_boot <- X_scaled[idx, ]
  
  km_boot <- kmeans(X_boot, centers = 2, nstart = 20)
  
  ari_values[b] <- adjustedRandIndex(
    df$cluster[idx],
    km_boot$cluster
  )
}

mean(ari_values)
sd(ari_values)

### clusters separation - SILHOUETTE #######
sil <- silhouette(df$cluster, dist(X_scaled))

mean(sil[, 3])

### between distance of centroids #########
centroids <- kmeans_fit$centers
centroid_distance <- dist(centroids)

centroid_distance
###########################################



########################################### TRAIN CLASSIFIER - logistic regression ###########################################
df$cluster_factor <- factor(df$cluster)

library(caret)
library(pROC)

set.seed(123)

control <- caret::trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = caret::twoClassSummary,
  savePredictions = TRUE
)

# Rename classes to valid names
levels(df$cluster_factor) <- c("HM", "Classic")


# df_scaled <- as.data.frame(scale(df[, c(
#   "n_cna_events_per_timing_group",
#   "frac_genome_affected_per_timing_group",
#   "ploidy",
#   "n_chr_affected_per_timing_group"
# )]))

# df_scaled$cluster_factor <- df$cluster_factor
df$cluster_factor <- relevel(df$cluster_factor, ref = "Classic")

model <- caret::train(
  cluster_factor ~ n_chr_affected_per_timing_group +
    frac_genome_affected_per_timing_group +
    n_cna_events_per_timing_group,
  data = df,
  method = "glm",
  family = binomial,
  trControl = control,
  metric = "ROC"
)

model
saveRDS(model, file = paste0(RES_FINAL_DIR,"00_HM_logistic_model.rds"))

######### extract decision rule #########
### we are modelling prob of being HM ###
summary(model$finalModel)

model$results$ROC
confusionMatrix(model)

################################################################