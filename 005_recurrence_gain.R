## -----------------------------------------------------------------------------
## Purpose of script: Plot recurrence of gained motifs
##
## Author: Oliver Artz
## Date Created: May 24, 2024
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
# annotated TCGA mutations
load("20240403_gained_motifs/processed_data/001_tcga_mc3_motifs_anno.RData")

# wrangle ----------------------------------------------------------------------
# filter for gained motifs
gained_motifs <- tcga_mc3_df_GRCh38 %>% 
  filter(motif_class == "gain")

# analysis ---------------------------------------------------------------------
# get stats
summary_stats <- gained_motifs %>% 
  filter(Variant_Classification %in% c("3'UTR", "Missense_Mutation", "Silent")) %>% 
  group_by(mut_motif_id, Variant_Classification) %>% 
  summarise(n = n()) %>%
  arrange(desc(n)) %>% 
  mutate(rec_category = case_when(n >= 4 ~ "4+",
                                  n == 3 ~ "3",
                                  n == 2 ~ "2",
                                  n == 1 ~ "1")) %>% 
  group_by(Variant_Classification, rec_category) %>%
  summarise(n = n())

# make df for plotting
df_plot <- summary_stats

# plot -------------------------------------------------------------------------
p005_recurrence <- df_plot %>% 
  ggplot(aes(x = reorder(rec_category, -n),
             y = n)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  geom_text(aes(label = n), vjust = -0.5, size = 3) +
  xlab("") +
  ylab("Number of gained motifs") +
  theme_bw() +
  theme(text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black")) +
  scale_fill_jco() +
  ylim(0, max(df_plot$n) * 1.1) +
  facet_wrap(~Variant_Classification) +
  labs(title = "Recurrence of gained motifs",
       x = "Number of recurrences",
       y = "Number of gained motifs")

p005_recurrence

# export -----------------------------------------------------------------------
ggsave("20240403_gained_motifs/plots/005_gain_recurrence.png", p005_recurrence, width = 8, height = 6, units = "in")

# cleanup ----------------------------------------------------------------------
rm(summary_stats, df_plot, tcga_mc3_df_GRCh38, gained_motifs)
