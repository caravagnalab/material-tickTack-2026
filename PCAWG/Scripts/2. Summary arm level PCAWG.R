## PCAWG summary plots 
library(dplyr)
library(ggplot2)
library(tidyr)
dir = "~/Dropbox/TickTack 2026"
segments = readRDS(paste0(dir, "/Tables/PCAWG/Segments.rds"))

## 1. Frequency of arm-level events by tumour type
plot1= segments %>%
  filter(cna_class!="focal",
         !(class %in% c("WGD"))) %>%
  mutate(cna_class = case_when(cna_class=="p_arm"~"p", cna_class=="q_arm"~"q",.default="")) %>%
  mutate(event = paste0(chr, cna_class)) %>%
  select(sample, event, clock_mean, class, Intogen_ttype, clock_clock_rank) %>%
  
  filter(!is.na(Intogen_ttype)) %>%
  group_by(Intogen_ttype) %>% mutate(n_samples = n_distinct(sample)) %>%
  filter(n_samples >=10) %>%
  #group_by(Intogen_ttype, sample) %>% mutate(n_ttype=n()) %>% filter(n_ttype > 1) %>%
  ungroup() %>%
  mutate(Intogen_ttype = paste0(Intogen_ttype, " (",n_samples,")")) %>%
  group_by(Intogen_ttype, event) %>% summarise(n_events = n()) %>%
  dplyr::arrange(desc(n_events)) %>%
  filter(n_events > 3) %>%
  ggplot(aes(x = reorder(event, n_events), y = n_events, fill = event)) + geom_col()+
  facet_wrap(~Intogen_ttype, scales = "free_x") + theme_minimal() + xlab("") + ylab("")+
  theme(legend.position = "none", 
        axis.text.x = element_text(angle=90, hjust = 1)) + 
  ggtitle("Most frequent arm-level and chromosome events in PCAWG (samples per ttype > 5, events > 3)") 

## 2. Rank of arm-level events by ttype
plot2=segments %>%
  filter(cna_class!="focal",
         !(class %in% c("WGD"))) %>%
  mutate(cna_class = case_when(cna_class=="p_arm"~"p", cna_class=="q_arm"~"q",.default="")) %>%
  mutate(event = paste0(chr, cna_class)) %>%
  select(sample, event, clock_mean, class, class, Intogen_ttype, clock_rank) %>%
  group_by(sample) %>% mutate(max_clock_rank=max(clock_rank)) %>% ungroup() %>%
  filter(!is.na(Intogen_ttype)) %>%
  group_by(Intogen_ttype) %>% mutate(n_samples = n_distinct(sample)) %>%
  filter(n_samples >=10) %>%
  #group_by(Intogen_ttype, sample) %>% mutate(n_ttype=n()) %>% filter(n_ttype > 1) %>%
  ungroup() %>%
  mutate(Intogen_ttype = paste0(Intogen_ttype, " (",n_samples,")")) %>%
  group_by(Intogen_ttype, event) %>% mutate(n_events = n()) %>%
  dplyr::arrange(desc(n_events)) %>%
  filter(n_events > 3) %>%
  ggplot(aes(x = reorder(event, clock_rank), y = clock_rank, color = event, fill= event, size)) + 
  geom_jitter(size=1, alpha=.5)+
  #geom_violin()+
  facet_grid(max_clock_rank~Intogen_ttype, scales = "free_x") + theme_bw() + xlab("") + ylab("")+
  theme(legend.position = "none", 
        axis.text.x = element_text(angle=90, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())





