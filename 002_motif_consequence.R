## -----------------------------------------------------------------------------
## Purpose of script: Plot motif consequences
##
## Author: Oliver Artz
## Date Created: 04-26-2024
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

# wrangle ----------------------------------------------------------------------
tcga_mc3_mod <- tcga_mc3_df_GRCh38

# analysis ---------------------------------------------------------------------
## plot all motif classes ----
# get stats
summary_stats <- tcga_mc3_mod %>% 
  group_by(motif_class) %>%
  summarise(n = n()/1000) %>%
  arrange(desc(n))

# in percent
summary_stats <- summary_stats %>% 
  mutate(percent = round(n/sum(n)*100, 1))

# plot -------------------------------------------------------------------------
p002_all_classes <- summary_stats %>% 
  ggplot(aes(x = reorder(motif_class, -n), y = n, fill = motif_class)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  geom_text(aes(label = paste0(round(n,0), " (", percent, "%)")), vjust = -0.5, size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black")) +
  labs(title = "Consequence of mutation on DRACH motifs",
       x = "",
       y = "Count (in thousands)") +
  scale_fill_jco() +
  ylim(0, max(summary_stats$n) * 1.1)

p002_all_classes

# export -----------------------------------------------------------------------
ggsave("20240403_gained_motifs/plots/002_all_classes.png", p002_all_classes, width = 8, height = 6, units = "in")

# cleanup ----------------------------------------------------------------------
rm(summary_stats, tcga_mc3_df_GRCh38)