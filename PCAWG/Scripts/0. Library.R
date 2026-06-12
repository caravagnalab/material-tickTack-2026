library(dplyr)
library(ggplot2)
library(tidyverse)
library(vegan)
library(ggrepel)
library(rstatix)    # for dunn_test post-hoc
library(patchwork)
library(factoextra)   # fviz_* helpers
library(cluster) 

plot_cnv_signatures_profiles = function(signatures_profiles){
  signatures = signatures_profiles %>% reshape2::melt() %>% pull(variable) %>% unique() %>% as.character() %>% sort()
  my_cols <- setNames(
    hcl.colors(length(signatures), "Dynamic"),
    levels(signatures)
  )
  names(my_cols) = signatures
  
  ## Signatures profile
  signatures_profiles %>% 
    reshape2::melt() %>%
    tidyr::separate(MutationType, into = c("Copy_number", "Type", "Size"), sep = ":") %>%
    mutate(
      Copy_number = factor(Copy_number, levels = c("0","1","2","3-4","5-8","9+")), 
      Type = factor(Type, levels = c("homdel", "LOH", "het")),
      Size = factor(Size, levels = c("0-100kb","100kb-1Mb", ">1Mb", "1Mb-10Mb","10Mb-40Mb",">40Mb"))
    ) %>%
    ggplot() +
    geom_col(aes(x = Size, y = value, fill = variable)) +
    facet_grid(variable~Type+Copy_number, scales = "free_y")+
    scale_fill_manual(values = my_cols) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.position = "none",
      panel.grid = element_blank()
    ) + xlab("") + ylab("")
}
plot_activities = function(df, group_by_class = T, relative_exp = F){
  
  df = df %>% reshape2::melt() %>% 
    group_by(sample) %>%
    mutate(max = sum(value))
  if (relative_exp) df = df %>% mutate(value = value/max) 
  
  signatures = df %>% pull(variable) %>% unique() %>% as.character() %>% sort()
  
  my_cols <- setNames(
    hcl.colors(length(signatures), "Dynamic"),
    levels(signatures)
  )
  names(my_cols) = signatures
  
  df %>% ungroup() %>%
    # rowwise() %>% mutate(sample = strsplit(sample, "-")[[1]][1]) %>% 
    # as_tibble() %>%
    ggplot(aes(x = forcats::fct_reorder(sample, max, .desc = TRUE), 
               y = value, fill = variable)) +
    geom_col() + 
    scale_fill_manual(values = my_cols) +
    facet_wrap(~class, scales = "free_x")+
    theme_bw() +
    theme(
      axis.text.x = element_blank(),#element_text(angle = 45, hjust = 1, size = 0),
      axis.ticks.x = element_blank(), 
      legend.position = "none",
      panel.grid = element_blank(),
    ) + xlab("") + ylab("") + labs(fill = NULL)
  
}
pca_signatures = function(df, samples_tt, class_colors){
  mat <- df %>% select(!c("class", "sample")) %>% as.matrix()
  rownames(mat) <- df$sample
  # Row-normalise so exposures are compositional (sum to 1 per sample)
  mat_norm <- mat / rowSums(mat)
  
  pca <- prcomp(mat_norm, scale. = F)
  pca_df   <- as.data.frame(pca$x) %>%
    rownames_to_column("sample") %>%
    left_join(samples_tt, by = "sample")
  var_exp  <- summary(pca)$importance[2, ] * 100   # % variance per PC
  # Loadings (which signatures drive each PC)
  loadings_df <- as.data.frame(pca$rotation) %>%
    rownames_to_column("signature")
  
  ggplot(pca_df, aes(PC1, PC2, color = class)) +
    geom_point(size = 3, alpha = 0.8) +
    # Confidence ellipses per class
    #stat_ellipse(aes(fill = class), geom = "polygon",
    #alpha = 0.1, level = 0.95, show.legend = FALSE) +
    # Loading arrows
    geom_segment(data = loadings_df,
                 aes(x = 0, y = 0,
                     xend = PC1 * 3, yend = PC2 * 3),   # scale factor for visibility
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "grey30", linewidth = 0.5, inherit.aes = FALSE) +
    geom_text_repel(data = loadings_df,
                    aes(x = PC1 * 3, y = PC2 * 3, label = signature),
                    color = "grey20", size = 3, inherit.aes = FALSE) +
    scale_color_manual(values = class_colors) +
    scale_fill_manual(values = class_colors) +
    labs(
      title = "PCA of CN signature profiles",
      x = sprintf("PC1 (%.1f%%)", var_exp[1]),
      y = sprintf("PC2 (%.1f%%)", var_exp[2]),
      color = "Class"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")
}
clustering = function(df, class_colors){
  # 1. Choose the number of k for k-means
  mat <- df %>% select(!c("class", "sample")) %>% as.matrix()
  rownames(mat) <- df$sample
  mat_norm <- mat / rowSums(mat)
  
  # 1a. Elbow (within-cluster sum of squares)
  # set.seed(59)
  k_opt_vector = c()
  wss_df = data.frame(
    "k" = numeric(),
    "value" = numeric(),
    "rep" = numeric()
  )
  for (i in 1:50){
    # p_elbow <- fviz_nbclust(mat_norm, kmeans, method = "wss", k.max = 10) +
      # labs(title = "Elbow method")
    res_elbow <- fviz_nbclust(mat_norm, kmeans, method = "wss", k.max = 10)
    wss <- res_elbow$data$y
    # wss_df[[as.character(i)]] = wss
    wss_df = rbind(
      wss_df,
      data.frame("k" = 1:10, "value" = wss, "rep" = rep(i, 10))
    )
    distances = diff(diff(wss))
    k_opt = which.max(distances[2:length(distances)]) + 1
    k_opt_vector[i] = k_opt
  }
  k_opt = round(mean(k_opt_vector))
  
  elbow_plot = wss_df %>% ggplot(aes(x = as.factor(k), y = value)) + 
    stat_summary(fun = median, geom = "point", size = 3) +
    stat_summary(fun.data = median_hilow, geom = "errorbar", 
                 width = 0.2, fun.args = list(conf.int = 0.95))+
    theme_bw() + xlab("K") + ylab("Within cluster sum of squares") #+ ggtitle("Within cluster sum of squares")
  
  # # 1b. Average silhouette width (higher = better separation)
  # p_sil <- fviz_nbclust(mat_norm, kmeans, method = "silhouette", k.max = 10) +
  #   labs(title = "Silhouette")
  # res <- fviz_nbclust(mat_norm, kmeans, method = "silhouette", k.max = 10)
  # # The optimal k is the one with maximum average silhouette width
  # k_opt <- as.integer(res$data$clusters[which.max(res$data$y)])
  
  # 2. Perform k-means
  km <- kmeans(mat_norm, centers = k_opt, nstart = 25, iter.max = 100)
  
  df_clusters <- data.frame(
    sample  = rownames(mat),      # or seq_len(nrow(mat)) if no rownames
    cluster = factor(km$cluster)
  )
  
  # 3. Visualize 
  signatures = colnames(df) %>% unique() %>% as.character() %>% sort()
  my_cols <- setNames(
    hcl.colors(length(signatures), "Dynamic"),
    levels(signatures)
  )
  names(my_cols) = signatures
  
  cluster_plot = df %>% left_join(df_clusters, by = "sample") %>%
    reshape2::melt() %>% 
    group_by(sample) %>%
    mutate(max = sum(value)) %>% 
    mutate(value = value/max) %>% ungroup() %>%
    ggplot(aes(x = forcats::fct_reorder(sample, max, .desc = TRUE), 
               y = value, fill = variable, color = class)) +
    geom_col() + 
    scale_fill_manual(values = my_cols) +
    scale_color_manual(values = class_colors)+
    facet_wrap(~cluster, scales = "free_x")+
    theme_bw() +
    theme(
      axis.text.x = element_blank(),#element_text(angle = 45, hjust = 1, size = 0),
      axis.ticks.x = element_blank(), 
      legend.position = "none",
      panel.grid = element_blank(),
    ) + xlab("") + ylab("") + labs(fill = NULL)
  
  elbow_plot + cluster_plot
}
Kruskal_Wallis_test = function(df, class_colors){
  mat <- df %>% select(!c("class", "sample")) %>% as.matrix()
  rownames(mat) <- df$sample
  mat_norm <- mat / rowSums(mat)
  
  sig_long <- as.data.frame(mat_norm) %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "signature", values_to = "exposure") %>%
    left_join(samples_tt, by = "sample")
  
  # Kruskal-Wallis per signature
  kw_res <- sig_long %>%
    group_by(signature) %>%
    summarise(
      statistic = kruskal.test(exposure ~ class)$statistic,
      p_kruskal = kruskal.test(exposure ~ class)$p.value,
      .groups = "drop"
    ) %>%
    mutate(p_adj = p.adjust(p_kruskal, method = "BH")) %>%
    arrange(p_adj)
  
  # print(kw_res)
  
  # Dunn post-hoc for each signature (pairwise class comparisons)
  dunn_res <- sig_long %>%
    group_by(signature) %>%
    dunn_test(exposure ~ class, p.adjust.method = "BH") %>%
    ungroup()
  
  # Violin plot per signature, annotated with KW significance
  
  # Significance labels
  sig_labels <- kw_res %>%
    mutate(sig_label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    ))
  
  sig_long %>%
    left_join(sig_labels %>% select(signature, sig_label, p_adj), by = "signature") %>%
    mutate(
      signature = fct_reorder(signature, p_adj)   # order by significance
    ) %>%
    ggplot(aes(x = class, y = exposure, fill = class)) +
    geom_violin(alpha = 1, trim = TRUE) +
    # geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.8) +
    geom_text(
      data = sig_labels %>% mutate(signature = fct_reorder(signature, p_adj)),
      aes(x = 2, y = Inf, label = sig_label),   # centred above
      vjust = 1.5, size = 5, inherit.aes = FALSE, color = "grey20"
    ) +
    facet_wrap(~ signature, scales = "free_y", nrow = 2) +
    scale_fill_manual(values = class_colors) +
    labs(
      title = "CN signature exposure by class",
      subtitle = "KW adjusted p-value: *** <0.001, ** <0.01, * <0.05, ns",
      x = NULL, y = "Normalised exposure",
      fill = "Class"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position    = "bottom",
      strip.background   = element_rect(fill = "grey92"),
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank()
    )
}