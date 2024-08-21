# ______________________________________________________________________________
# Plot how often motifs are mutated
# ______________________________________________________________________________

# clear workspace --------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, ggsci, readxl)

# load data --------------------------------------------------------------------
mutations_in_m6A_motifs <- read_excel("20230322_publication_figs/processed_data/001.1_mutations_in_m6A_motifs.xlsx")
load("20230322_publication_figs/processed_data/001_tcga_mc3_df_GRCh38.RData")

# wrangle ----------------------------------------------------------------------
# count how often specific motifs are mutated
count_recurrence_motif_mut <- mutations_in_m6A_motifs %>% 
  group_by(motif_position, Variant_Classification, Hugo_Symbol) %>%
  tally() 

# refactor n
count_recurrence_motif_mut$n[count_recurrence_motif_mut$n >= 4] <- "4+"

# count frequency of recurrences
count_frequency <- count_recurrence_motif_mut %>% 
  group_by(Variant_Classification, n) %>% 
  summarize(n_frequency = n())

# plot -------------------------------------------------------------------------
df_plot <- count_frequency %>% 
  filter(Variant_Classification %in% c("3'UTR", "Missense_Mutation", "Silent"))

# define ylims
ylim_df <- count_frequency %>% 
  group_by(Variant_Classification) %>% 
  slice_max(n_frequency, n = 1) %>% 
  mutate(max_freq = n_frequency * 1.1)

# plot
p_006 <- df_plot %>% 
  ggplot(aes(x = n, y = n_frequency)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  geom_text(aes(label = n_frequency), vjust = -0.25, size = 3) +
  theme_bw() +
  theme(text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black")) +
  xlab("Number of recurrences") +
  ylab("Number m6A motifs") +
  scale_fill_jco() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  facet_wrap(~ Variant_Classification, 
             scales = "free")

p_006

# export
ggsave(filename = "20230322_publication_figs/plots/006_plot_recurrence_mutated_motifs.png", width = 12, height = 6)

# cleanup
rm(count_frequency, count_recurrence_motif_mut, df_plot, mutations_in_m6A_motifs, tcga_mc3_df_GRCh38, ylim_df)
