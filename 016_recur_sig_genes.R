## -----------------------------------------------------------------------------
## Purpose of script: Plot number of recurrences for significant genes
##
## Author: Oliver Artz
## Date Created: Jun 24, 2024
## Email: artzo@mskcc.org
##
## -----------------------------------------------------------------------------

# clear work space -------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, data.table)

# set parameters ---------------------------------------------------------------

# load data --------------------------------------------------------------------
mutation_info_list <- fread("20240403_gained_motifs/processed_data/015_mutation_info_list.csv")

# wrangle ----------------------------------------------------------------------
mutation_info_list_mod <- mutation_info_list %>% 
  mutate(cancer_type = cancer_type %>% str_remove("TCGA-"))

# make columns for axis label
mutation_info_list_mod <- mutation_info_list_mod %>% 
  mutate(mut_label = paste0(Hugo_Symbol, " - ", cancer_type))

# analysis ---------------------------------------------------------------------

# plot -------------------------------------------------------------------------
p_016_recur_sig_genes <- mutation_info_list_mod %>% 
  ggplot(aes(x = total_mutations, y = reorder(mut_label, total_mutations), fill = cancer_type)) +
  geom_col(alpha = 0.8, color = "black", position = position_dodge()) +
  geom_text(aes(label = total_mutations), position = position_dodge(width = 0.9), hjust = -0.2, size = 3) +
  theme_bw() +
  theme(text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text.y = element_text(size = 6)) +
  labs(title = "Number of recurrences for significant genes",
       x = "Number of recurrences",
       y = "",
       fill = "Cancer type")

p_016_recur_sig_genes

# export -----------------------------------------------------------------------
ggsave("20240403_gained_motifs/plots/016_recur_sig_genes.png", width = 12, height = 6)

# cleanup ----------------------------------------------------------------------
rm(mutation_info_list, mutation_info_list_mod)
