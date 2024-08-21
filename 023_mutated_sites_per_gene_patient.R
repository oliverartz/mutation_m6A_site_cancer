## -----------------------------------------------------------------------------
## Purpose of script: Plotting number of mutated m6A sites for each gene per patient
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

p_load(tidyverse, ggpubr, data.table, ggsci, ggpmisc)

# load data --------------------------------------------------------------------
load("20230322_publication_figs/processed_data/001_tcga_mc3_df_GRCh38.RData")

# analysis ---------------------------------------------------------------------
# filter mutational database
tcga_mc3_df_GRCh38_filtered <- tcga_mc3_df_GRCh38 %>%
  filter(Variant_Classification %in% c("3'UTR", "Missense_Mutation", "Silent"))

# count m6A mutations
count_m6A_mut_gene_patient <- tcga_mc3_df_GRCh38_filtered %>% 
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>% 
  summarize(meth_mut = sum(in_methylated_motif == 1),
            total_mut = n()) %>% 
  mutate(methylated_total = meth_mut / total_mut) %>% 
  drop_na() %>% 
  filter(meth_mut > 0)

# prop table
prop_n_mut <- table(count_m6A_mut_gene_patient$meth_mut) %>% 
  prop.table() %>% 
  as.data.frame() %>% 
  mutate(Freq = round(Freq, 4) * 100)

# plot -------------------------------------------------------------------------
p_023 <- prop_n_mut %>% 
  ggplot(aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "black") +
  geom_text(aes(label = Freq), vjust = -0.5) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(color = "black"),
              axis.text = element_text(color = "black"),
              axis.title = element_text(color = "black")) +
  labs(
    x = "Mutations per gene per patient",
    y = "% genes with mutated m6A site") +
  scale_fill_jco() +
  ylim(0, 105)

p_023

# export -----------------------------------------------------------------------
ggsave(filename = "20230322_publication_figs/plots/023_mutated_sites_per_gene_patient.png",
      height = 6)
      
    
# cleanup ----------------------------------------------------------------------
rm(count_m6A_mut_gene_patient, prop_n_mut, tcga_mc3_df_GRCh38, tcga_mc3_df_GRCh38_filtered)