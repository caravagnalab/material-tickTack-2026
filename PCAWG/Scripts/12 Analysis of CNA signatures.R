setwd("~/dati_Orfeo/antonelloa/material-tickTack-2026/PCAWG/Scripts/")
source("0. Library.R")

samples = readRDS("~/dati_Orfeo/antonelloa/material-tickTack-2026/PCAWG/Data/Samples.rds")
ttypes = list.files("~/dati_Orfeo/antonelloa/material-tickTack-2026/PCAWG/Analysis_results/CNA_Signatures/Sigprofiler") #samples %>% filter(!is.na(IntoGen_cancer_type)) %>% pull(IntoGen_cancer_type) %>% unique()
# t = ttypes[2]
for (t in ttypes){
  print(t)
  sigprofiler_dir = "~/dati_Orfeo/antonelloa/material-tickTack-2026/PCAWG/Analysis_results/CNA_Signatures/Sigprofiler/"
  output = path.expand(paste0(sigprofiler_dir, t, "/CNV48/Suggested_Solution/COSMIC_CNV48_Decomposed_Solution"))
  signatures_profiles = read.table(paste0(output, "/Signatures/COSMIC_CNV48_Signatures.txt"), header = 1)
  activities = read.table(paste0(output, "/Activities/COSMIC_CNV48_Activities.txt"), header = 1)
  class_colors <- c("WGD" = "firebrick", "Classic" = "darkgrey", "HM" = "steelblue")
  
  cnsig_profiles = plot_cnv_signatures_profiles(signatures_profiles)
  samples_tt = samples %>% filter(IntoGen_cancer_type == t) %>% select(sample, class) %>% distinct() 
  activities = activities %>% rename(sample = Samples)
  df = samples_tt %>% left_join(activities, by = "sample")
  abs_exp = plot_activities(df, relative_exp = F) 
  rel_exp = plot_activities(df, relative_exp = T) 
  
  final_plot1 = patchwork::wrap_plots(cnsig_profiles, abs_exp, rel_exp,
                        design = "AAAAABBB
                                  AAAAABBB
                                  AAAAACCC
                                  AAAAACCC")
  
  # ── 1. PCA ────────────────────────────────────────────────────────────────────
  final_plot2 = pca_signatures(df, samples_tt, class_colors)
  
  # ── 2. Unsupervised Clustering ────────────────────────────────────────────────
  final_plot3 = clustering(df, class_colors) + ggtitle("K-means clustering")
  
  # ── 3. Kruskal-Wallis + Dunn test ─────────────────────────────────────────────
  final_plot4 = Kruskal_Wallis_test(df, class_colors)
  
  # Save final Plot
  out_dir_plots = "~/dati_Orfeo/antonelloa/material-tickTack-2026/PCAWG/Analysis_results/CNA_Signatures/Plots/"
  patchwork::wrap_plots(
    final_plot1, final_plot3, final_plot2, final_plot4,
    design = "AA
              AA
              BB
              C#
              DD") %>% ggsave(filename = paste0(out_dir_plots, t, ".pdf"), height = 20, width = 12)
}
