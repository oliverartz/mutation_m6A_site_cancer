# ______________________________________________________________________________
# Number of methylated genes per cancer type above threshold with loss
# Genes with X number of mutated m6A sites leading to loss
# Count number of genes per cancer type
# ______________________________________________________________________________

# clear workspace --------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, ggsci, ggpubr, cowplot)

# load data --------------------------------------------------------------------
mutations_in_m6A_motifs <- readxl::read_xlsx("20230322_publication_figs/processed_data/001.1_mutations_in_m6A_motifs.xlsx")

# analysis ---------------------------------------------------------------------
# filter for mutations in DRACH that lead to loss 
mutations_in_m6A_motifs_filtered <- mutations_in_m6A_motifs %>%
  filter(in_DRACH == 1) %>% 
  filter(motif_consequence == "loss") %>% 
  filter(Variant_Classification %in% c("3'UTR", "Missense_Mutation", "Silent"))


# count number of genes per cancer type
count_genes <- mutations_in_m6A_motifs_filtered %>% 
  group_by(Hugo_Symbol, cancer_type) %>% 
  tally()

# define threshold for recurrence
threshold_recurrence <- 2

# filter above threshold recurrences and count number of genes per cancer type
count_genes_filtered <- count_genes %>% 
  filter(n >= threshold_recurrence) %>% 
  group_by(cancer_type) %>% 
  tally() %>% 
  na.omit()

# remove 'TCGA-' prefix
count_genes_filtered$cancer_type <- str_remove(count_genes_filtered$cancer_type, "TCGA-")

summary(count_genes_filtered$n)

# plot -------------------------------------------------------------------------
p_010 <- count_genes_filtered %>% 
  ggplot(aes(
    x = n,
    y = reorder(cancer_type, n, sum))) +
  geom_bar(stat = "identity", alpha = 0.8, color = "black") +
  theme_bw() +
  scale_fill_jco() +
  xlab("Number of genes") +
  ylab("") +
  labs(fill = "Variant type") +
  theme(text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"))

p_010

# export -----------------------------------------------------------------------
ggsave(filename = "20230322_publication_figs/plots/010_methylated_genes_per_cancer_type_above_threshold_with_loss.png", width = 12, height = 6)

# cleanup ----------------------------------------------------------------------
rm(count_genes, count_genes_filtered, mutations_in_m6A_motifs, mutations_in_m6A_motifs_filtered, threshold_recurrence)
