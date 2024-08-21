# ______________________________________________________________________________
# Plot mutations in methylated motifs
# DRACH motifs that have been shown to be methylated
# ______________________________________________________________________________

# clear workspace --------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
options(scipen = 999)
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, readxl, ggsci, cowplot)

# load data --------------------------------------------------------------------
load("20230322_publication_figs/processed_data/001_tcga_mc3_df_GRCh38.RData")

# wrangle ----------------------------------------------------------------------
# filter for DRACH motifs
only_DRACH_motifs <- tcga_mc3_df_GRCh38 %>% filter(in_DRACH == 1)

# remove 'TCGA-' prefix
only_DRACH_motifs <- only_DRACH_motifs %>% 
  mutate(cancer_type = str_remove(cancer_type, "TCGA-"))

# absolute number of mutated methylated DRACH motifs per cancer type
abs_mut_DRACH <- only_DRACH_motifs %>% 
  group_by(cancer_type, in_methylated_motif) %>% 
  tally() %>% 
  drop_na()

# proportion of all mutations in DRACH sites per cancer type
# NOT all mutations per cancer type
rel_mut_DRACH <- table(tcga_mc3_df_GRCh38$cancer_type, tcga_mc3_df_GRCh38$in_methylated_motif) %>% 
  prop.table(1) %>% 
  as.data.frame() %>%
  mutate(Freq = round(Freq * 100, digits = 1)) %>% 
  setNames(c("cancer_type", "in_methylated_motif", "percent")) %>% 
  filter(in_methylated_motif == 1)

df_temp <- abs_mut_DRACH %>% filter(in_methylated_motif == 1)
summary(df_temp$n)


# plot -------------------------------------------------------------------------
# absolute number
p_abs <- abs_mut_DRACH %>% 
  filter(in_methylated_motif == 1) %>% 
  drop_na() %>% 
  ggplot(aes(x = reorder(cancer_type, n),
             y = n)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8, fill = "grey") +
  xlab("") +
  ylab("Number of mutations in methylated DRACH motifs") +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
    #    axis.text = element_text(size = 6),
    plot.margin = margin(5.5, 0, 5.5, 5.5)) +
  scale_fill_jco()

# relative number
p_rel <- rel_mut_DRACH %>% 
  ggplot(aes(x = cancer_type, abs_mut_DRACH$n[abs_mut_DRACH$in_methylated_motif == 1],
             y = percent)) +
  geom_hline(yintercept = mean(rel_mut_DRACH$percent), linetype = "dashed") +
  geom_bar(stat = "identity", color = "black", alpha = 0.8, fill = "grey") +
  xlab("") +
  ylab("% of methylated DRACHs") +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(5.5, 5.5, 5.5, 0)) +
  scale_fill_jco()

# combine plots
# p_004 <- plot_grid(p_abs, 
#                    p_rel, nrow = 1, rel_widths = c(4,1))

p_004 <- p_abs

p_004

# export plot
ggsave(filename = "20230322_publication_figs/plots/004_plot_mutation_in_methylated_motifs.png", height = 6, width = 12)
