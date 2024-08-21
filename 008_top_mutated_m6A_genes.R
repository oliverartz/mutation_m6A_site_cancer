# ______________________________________________________________________________
# Plot genes and motifs with most mutations in m6A sites
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
# count mutations per gene per Variant_Classification
mutations_in_m6A_motifs_filtered <- mutations_in_m6A_motifs %>%
  filter(in_DRACH == 1) %>% 
  filter(motif_consequence == "loss") %>% 
  filter(Variant_Classification %in% c("3'UTR", "Missense_Mutation", "Silent"))

counts_gene_classif <- mutations_in_m6A_motifs_filtered %>%
  group_by(Hugo_Symbol, Variant_Classification) %>% 
  tally() %>% 
  arrange(desc(n)) %>% 
  dplyr::rename(n_gene_classif = n)

# count mutations per gene
counts_gene <- mutations_in_m6A_motifs_filtered %>%
  group_by(Hugo_Symbol) %>% 
  tally() %>% 
  arrange(desc(n)) %>% 
  dplyr::rename(n_gene = n)

# merge data
counts <- left_join(counts_gene_classif, counts_gene, by = "Hugo_Symbol")

# cutoff for counts
counts_plot <- counts %>% 
  filter(n_gene > 9)
  

# plot -------------------------------------------------------------------------
p_008 <- counts_plot %>% 
  ggplot(aes(x = reorder(Hugo_Symbol, -n_gene), y = n_gene_classif, fill = Variant_Classification)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  geom_text(
    data = . %>% distinct(Hugo_Symbol, .keep_all = TRUE),
    aes(label = Hugo_Symbol),
    size = 3,
    angle = 90,
    y = max(counts_plot$n_gene),
    hjust = -0.1
  ) +
  theme_bw() +
  scale_fill_jco() +
  labs(x = "Gene",
       y = "Number of mutations",
       fill = "Variant Classification") +
  theme(axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(),
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black")) +
  ylim(0, max(counts_plot$n_gene * 1.3))

p_008

# export -----------------------------------------------------------------------
ggsave(filename = "20230322_publication_figs/plots/008_top_mutated_m6A_genes.png", width = 12, height = 6)

# cleanup
rm(counts, counts_gene, counts_gene_classif, counts_plot, mutations_in_m6A_motifs, mutations_in_m6A_motifs_filtered)
