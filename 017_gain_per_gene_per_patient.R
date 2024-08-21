## -----------------------------------------------------------------------------
## Purpose of script: Check number of gains per gene per patient
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
## motif data
load("20240403_gained_motifs/processed_data/001_tcga_mc3_motifs_anno.RData")

# wrangle ----------------------------------------------------------------------

# analysis ---------------------------------------------------------------------
## number of gains per gene per patient
n_gain_gene_patient <- tcga_mc3_df_GRCh38 %>% 
  filter(motif_class == "gain") %>% 
  group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>%
  summarise(n_gain = n())

## count 
n_gain <- n_gain_gene_patient %>% 
  group_by(n_gain) %>% 
  summarize(n_patients = n()) %>% 
  drop_na() %>% 
  mutate(total = sum(n_patients),
         percentage = (n_patients / total) * 100)


# plot -------------------------------------------------------------------------
p_017_gain_per_gene_per_patient <- n_gain %>% 
  ggplot(aes(x = n_gain, y = n_patients)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "black") +
  geom_text(aes(label = n_patients), vjust = -0.5, size = 3) +
  labs(
    title = "Gains per gene per patient",
    x = "Number of gains",
    y = "Number of patients"
  ) +
  theme_bw() +
  theme(text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black")) +
  scale_x_continuous(breaks = seq(0, 10, 1))

# export -----------------------------------------------------------------------
ggsave("20240403_gained_motifs/plots/017_gain_per_gene_per_patient.png", width = 12, height = 6)

# cleanup ----------------------------------------------------------------------
rm(n_gain, n_gain_gene_patient, tcga_mc3_df_GRCh38)

