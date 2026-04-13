## PCAWG summary plots 
library(dplyr)
library(ggplot2)
library(tidyr)
dir = "~/Dropbox/TickTack 2026"
drivers = readRDS(paste0(dir, "/Tables/PCAWG/Drivers.rds")) %>% select(
  "segment_id","gene","karyotype","sample_id",
  "mutatation_status","mult_estimate",
  "timed", "sample" ,"clock_mean" , "clock_rank" ,"class",
  "best_K", "IntoGen_cancer_type") 

## 1. Frequency of amplified drivers per ttype (colored if mutated)
drivers %>%
  filter(timed == "TRUE", mutatation_status!="CI_M", class == "Classic") %>%
  mutate(status = case_when(mutatation_status == "M" & mult_estimate == 2 ~ "Mut+Amp",
                            mutatation_status == "M" & mult_estimate == 1 ~ "Mut",
                            .default="Amp"
                            )) %>%
  group_by(IntoGen_cancer_type) %>% 
  mutate(n_samples = n_distinct(sample)) %>% 
  mutate(IntoGen_cancer_type = paste0(IntoGen_cancer_type, "(",n_samples,")")) %>% ungroup() %>%
  group_by(IntoGen_cancer_type, gene) %>% mutate(n_gene = n()) %>% ungroup() %>%
  filter(n_samples > 5, n_gene > 5) %>%
  #group_by(IntoGen_cancer_type) %>% mutate(n_by_ttype = n_distinct()) %>% 
  #ungroup() %>% filter(n_by_ttype > 5) %>%
  #ungroup() %>% group_by(IntoGen_cancer_type) %>% mutate(f_by_ttype = n_by_ttype/n()) %>%
  #ungroup() %>% filter(f_by_ttype > .05) %>% 
  ggplot(aes(x = reorder(gene, n_gene), fill = status)) + geom_bar() + 
  facet_wrap(~IntoGen_cancer_type, scales = "free")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Frequency of Amplified drivers")+xlabs("")+ylab("")

drivers %>%
  filter(timed == "TRUE", mutatation_status!="CI_M", class == "HM") %>%
  mutate(status = case_when(mutatation_status == "M" & mult_estimate == 2 ~ "Mut+Amp",
                            mutatation_status == "M" & mult_estimate == 1 ~ "Mut",
                            .default="Amp"
  )) %>%
  group_by(IntoGen_cancer_type) %>% 
  mutate(n_samples = n_distinct(sample)) %>% 
  mutate(IntoGen_cancer_type = paste0(IntoGen_cancer_type, "(",n_samples,")")) %>% ungroup() %>%
  group_by(IntoGen_cancer_type, gene) %>% mutate(n_gene = n()) %>% ungroup() %>%
  #filter(n_samples > 5, n_gene > 5) %>%
  ggplot(aes(x = reorder(gene, n_gene), fill = status)) + geom_bar() + 
  facet_wrap(~IntoGen_cancer_type, scales = "free")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Frequency of Amplified drivers (HM)")+xlabs("")+ylab("")

## 3. Timing of amplified drivers (coloured if mutated)
