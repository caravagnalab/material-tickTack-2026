## PCAWG summary plots 
library(dplyr)
library(ggplot2)
library(tidyr)
dir = "~/Docs/GitHub/material-tickTack-2026/PCAWG/Data/"
segments = readRDS(paste0(dir, "Segments.rds"))

# ## 1. Frequency of arm-level events by tumour type
# plot1= segments %>%
#   filter(arm_level_class!="focal",
#          !(class %in% c("WGD"))) %>%
#   mutate(arm_level_class = case_when(arm_level_class=="p_arm"~"p", arm_level_class=="q_arm"~"q",.default="")) %>%
#   mutate(event = paste0(chr, arm_level_class)) %>%
#   select(sample, event, clock_mean, class, IntoGen_cancer_type, best_K) %>%
#   
#   filter(!is.na(IntoGen_cancer_type)) %>%
#   group_by(IntoGen_cancer_type) %>% mutate(n_samples = n_distinct(sample)) %>%
#   filter(n_samples >=5) %>%
#   #group_by(IntoGen_cancer_type, sample) %>% mutate(n_ttype=n()) %>% filter(n_ttype > 1) %>%
#   ungroup() %>%
#   mutate(IntoGen_cancer_type = paste0(IntoGen_cancer_type, " (",n_samples,")")) %>%
#   group_by(IntoGen_cancer_type, event) %>% summarise(n_events = n()) %>%
#   dplyr::arrange(desc(n_events)) %>%
#   filter(n_events > 3) %>%
#   ggplot(aes(x = reorder(event, n_events), y = n_events, fill = event)) + geom_col()+
#   facet_wrap(~IntoGen_cancer_type, scales = "free") + theme_minimal() + xlab("") + ylab("")+
#   theme(legend.position = "none", 
#         axis.text.x = element_text(angle=90, hjust = 1)) + 
#   ggtitle("Most frequent arm-level and chromosome events in PCAWG (samples per ttype > 5, events > 3)") 
# 
# ## 2. Rank of arm-level events by ttype
# plot2=segments %>%
#   filter(arm_level_class!="focal",
#          !(class %in% c("WGD"))) %>%
#   mutate(arm_level_class = case_when(arm_level_class=="p_arm"~"p", arm_level_class=="q_arm"~"q",.default="")) %>%
#   mutate(event = paste0(chr, arm_level_class)) %>%
#   select(sample, event, clock_mean, class, class, IntoGen_cancer_type, clock_rank) %>%
#   group_by(sample) %>% mutate(max_clock_rank=max(clock_rank)) %>% ungroup() %>%
#   filter(!is.na(IntoGen_cancer_type)) %>%
#   group_by(IntoGen_cancer_type) %>% mutate(n_samples = n_distinct(sample)) %>%
#   filter(n_samples >=5) %>%
#   #group_by(IntoGen_cancer_type, sample) %>% mutate(n_ttype=n()) %>% filter(n_ttype > 1) %>%
#   ungroup() %>%
#   mutate(IntoGen_cancer_type = paste0(IntoGen_cancer_type, " (",n_samples,")")) %>%
#   group_by(IntoGen_cancer_type, event) %>% mutate(n_events = n()) %>%
#   dplyr::arrange(desc(n_events)) %>%
#   filter(n_events > 3) %>%
#   ggplot(aes(x = reorder(event, clock_rank), y = clock_rank, color = event, fill= event, size)) + 
#   geom_jitter(size=1, alpha=.5)+
#   #geom_violin()+
#   facet_grid(max_clock_rank~IntoGen_cancer_type, scales = "free_x") + theme_bw() + xlab("") + ylab("")+
#   theme(legend.position = "none", 
#         axis.text.x = element_text(angle=90, hjust = 1),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor = element_blank())


library(tidyverse)

df_agg <- segments %>%
  filter(arm_level_class!="focal",
         !(class %in% c("WGD"))) %>%
  mutate(arm_level_class = case_when(arm_level_class=="p_arm"~"p", arm_level_class=="q_arm"~"q",.default="")) %>%
  mutate(event = paste0(chr, arm_level_class)) %>%
  select(sample, event, clock_mean, class, IntoGen_cancer_type, clock_rank) %>%
  group_by(sample) %>% mutate(max_clock_rank=max(clock_rank)) %>% ungroup() %>%
  filter(!is.na(IntoGen_cancer_type)) %>%
  #group_by(IntoGen_cancer_type) %>% mutate(n_samples = n_distinct(sample)) %>%
  group_by(event, IntoGen_cancer_type) %>%
  summarise(
    n_samples  = n_distinct(sample),
    clock_med  = median(clock_mean),
    modal_rank = as.integer(names(sort(table(clock_rank), decreasing=TRUE))[1]),
    modal_class = names(sort(table(class), decreasing=TRUE))[1],
    .groups = "drop"
  ) %>%
  filter(n_samples >= 3)   # tune threshold

# rank_colors <- c("1"="#e8853a","2"="#4aaf82","3"="#3b78c0","4"="#9044c5","5"="#c04040")

summary_arm_level_plot = ggplot(df_agg, aes(
  x     = IntoGen_cancer_type,
  y     = fct_reorder(event, clock_med),   # sort arms by pseudotime
  size  = n_samples,
  color = clock_med,
  shape = factor(modal_rank)  # secondary encoding for rank, avoids color-only
)) +
  geom_point(alpha = 0.85) +
  scale_color_gradientn(
    colours = c("#EF9F27","#D05538","#3266ad","#185FA5"),
    name    = "Pseudotime\n(median)"
  ) +
  scale_size_area(max_size = 14, name = "Samples") +
  scale_shape_manual(values = c(19,17,15,18,8), name = "Clock rank") +
  facet_wrap(~ modal_class, nrow = 1) +   # split panels by Bayesian class
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(linewidth = 0.3),
    legend.position = "bottom"
  ) +
  labs(x = NULL, y = NULL,
       title = "CNA event timing by tumour type",
       subtitle = "Size = sample frequency · colour = pseudotime · shape = clock rank") 

ggsave(plot = summary_arm_level_plot, filename = "~/Docs/GitHub/material-tickTack-2026/PCAWG/Plot/Summary_arm_level_events.png", height = 15, width = 15)





