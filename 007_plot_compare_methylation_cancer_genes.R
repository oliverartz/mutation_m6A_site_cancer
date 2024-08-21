# ______________________________________________________________________________
# Plot comparison of methylated vs non-methylated cancer genes
# Comparing mutations in DRACH motifs that are either methylated or not methylated
# ______________________________________________________________________________

# clear workspace --------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, ggsci, ggpubr, cowplot)

# load data --------------------------------------------------------------------
load("20230322_publication_figs/processed_data/001_tcga_mc3_df_GRCh38.RData")
mutations_in_m6A_motifs <- readxl::read_xlsx("20230322_publication_figs/processed_data/001.1_mutations_in_m6A_motifs.xlsx")


# wrangle ----------------------------------------------------------------------
# define new labels for facets
facet_labels <- c("3'UTR", "missense", "silent")
names(facet_labels) <- c("3'UTR", "Missense_Mutation", "Silent")

# reorder x-axis for plot
cancer_gene_order <- c("OG", "TS", "OG & TS", "n.d.")

# refactor 'not in OncoKB'
tcga_mc3_df_GRCh38$cancer_gene[tcga_mc3_df_GRCh38$cancer_gene == "not in OncoKB"] <- "n.d."
mutations_in_m6A_motifs$cancer_gene[mutations_in_m6A_motifs$cancer_gene == "not in OncoKB"] <- "n.d."

# refactor in_methylated_motif
tcga_mc3_df_GRCh38$in_methylated_motif[tcga_mc3_df_GRCh38$in_methylated_motif == 0] <- "no m6A"
tcga_mc3_df_GRCh38$in_methylated_motif[tcga_mc3_df_GRCh38$in_methylated_motif == 1] <- "m6A"

mutations_in_m6A_motifs$in_methylated_motif[mutations_in_m6A_motifs$in_methylated_motif == 0] <- "no m6A"
mutations_in_m6A_motifs$in_methylated_motif[mutations_in_m6A_motifs$in_methylated_motif == 1] <- "m6A"

## gene_level comparison ----
gene_level_recurrence <- tcga_mc3_df_GRCh38 %>% 
  filter(Variant_Classification %in% c("Missense_Mutation", "Silent", "3'UTR")) %>% 
  filter(in_DRACH == 1) %>% 
  filter(cancer_gene %in% c("n.d.", "TS", "OG", "OG & TS")) %>% 
  group_by(Hugo_Symbol, in_methylated_motif, cancer_gene, Variant_Classification) %>% 
  tally() %>% 
  mutate(in_methylated_motif = as.factor(in_methylated_motif))

## motif_level comparison ----
motif_level_recurrence <- tcga_mc3_df_GRCh38 %>% 
  filter(Variant_Classification %in% c("Missense_Mutation", "Silent", "3'UTR")) %>% 
  filter(in_DRACH == 1) %>% 
  filter(cancer_gene %in% c("n.d.", "TS", "OG", "OG & TS")) %>% 
  group_by(motif_id, in_methylated_motif, cancer_gene, Variant_Classification) %>%
  tally() %>% 
  mutate(in_methylated_motif = as.factor(in_methylated_motif))

## base-level comparison ----
base_level_recurrence <- tcga_mc3_df_GRCh38 %>% 
  filter(Variant_Classification %in% c("Missense_Mutation", "Silent", "3'UTR")) %>% 
  filter(in_DRACH == 1) %>% 
  filter(cancer_gene %in% c("n.d.", "TS", "OG", "OG & TS")) %>%
  mutate(mutation_id = paste0(seqnames, ":", start, "_", STRAND)) %>% 
  group_by(mutation_id, in_methylated_motif, cancer_gene, Variant_Classification) %>% 
  tally() %>% 
  mutate(in_methylated_motif = as.factor(in_methylated_motif))

# plot -------------------------------------------------------------------------

# define plotting function
plot_recurrence_meth <- function(dataset, title){
  dataset %>%
  ggplot(aes(x = fct_relevel(cancer_gene, rev(cancer_gene_order)),
             y = n, 
             fill = in_methylated_motif,
             colour = in_methylated_motif
  )) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.3),
    alpha = 0.3,
    shape = 21,
    size = 1
  ) +
  theme_bw() +
  theme(text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black")) +
  xlab("") +
  ylab("Recurrence") +
  scale_fill_jco() +
  scale_color_jco() +
  guides(fill = guide_legend(title = "DRACH motif with", reverse = TRUE),
         color = "none") +
  stat_compare_means(aes(label = after_stat(p.signif)), vjust = 0.5, method = "wilcox.test", size = 3) +
  coord_flip() +
  ggtitle(title) +
  facet_wrap(Variant_Classification ~ ., ncol = 1, strip.position = "right",  labeller = labeller(Variant_Classification = facet_labels), scales = "free")
}

# call plotting function
p1 <- plot_recurrence_meth(gene_level_recurrence, "Gene-level")
p2 <- plot_recurrence_meth(motif_level_recurrence, "Motif-level")
p3 <- plot_recurrence_meth(base_level_recurrence, "Base-level")

# extract legend
# legend <- get_legend(p1 + theme(legend.box.margin = margin(0, 0, 0, 0)))


# remove legends
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "right")

# make composite plot
pg_1 <- plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(1,1,1.4))
p_007 <- pg_1

p_007

# p_007 <- plot_grid(pg_1, legend, rel_widths = c(5, 1), nrow = 1) +
#   theme(plot.background = element_rect(fill = "white"))



# export
ggsave(filename = "20230322_publication_figs/plots/007_plot_compare_methylation_cancer_genes.png", width = 12, height = 6)

# cleanup
rm(base_level_recurrence, gene_level_recurrence, motif_level_recurrence, mutations_in_m6A_motifs, p1, p2, p3, pg_1, tcga_mc3_df_GRCh38, cancer_gene_order, facet_labels, plot_recurrence_meth)
