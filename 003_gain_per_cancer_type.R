## -----------------------------------------------------------------------------
## Purpose of script: Plot the number of gained motifs per cancer type
##
## Author: Oliver Artz
## Date Created:
## Email: artzo@mskcc.org
##
## -----------------------------------------------------------------------------

# clear work space -------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, data.table, ggpubr, ggsci)

# load data --------------------------------------------------------------------
# annotated TCGA mutations
load("20240403_gained_motifs/processed_data/001_tcga_mc3_motifs_anno.RData")

# set parameters ---------------------------------------------------------------

# load data --------------------------------------------------------------------

# wrangle ----------------------------------------------------------------------
# filter for gained motifs
gained_motifs <- tcga_mc3_df_GRCh38 %>% 
  filter(motif_class == "gain")

# analysis ---------------------------------------------------------------------
# get stats
summary_stats <- gained_motifs %>% 
  group_by(cancer_type) %>%
  summarise(n = n()) %>% 
  mutate(cancer_type = str_remove(cancer_type, "TCGA-"))

# stats for publication text
# median(summary_stats$n)


# plot -------------------------------------------------------------------------
p003_cancer_types <- summary_stats %>% 
  drop_na() %>% 
  ggplot(aes(x = reorder(cancer_type, n),
             y = n)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8, fill = "grey") +
  xlab("") +
  ylab("Number of gained motifs") +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        plot.margin = margin(5.5, 0, 5.5, 5.5)) +
  scale_fill_jco()

p003_cancer_types

# export -----------------------------------------------------------------------
ggsave("20240403_gained_motifs/plots/003_gain_cancer_type.png", p003_cancer_types, width = 8, height = 6, units = "in")

# cleanup ----------------------------------------------------------------------
rm(gained_motifs, summary_stats)