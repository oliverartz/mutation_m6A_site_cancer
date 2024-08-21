# ______________________________________________________________________________
# Plot total mutations vs mutations in m6A sites
# ______________________________________________________________________________

# clear workspace --------------------------------------------------------------
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
count_m6A_mut <- tcga_mc3_df_GRCh38_filtered %>% 
  group_by(Tumor_Sample_Barcode, cancer_type) %>% 
  summarize(meth_mut = sum(in_methylated_motif == 1),
            total_mut = n()) %>% 
  mutate(methylated_total = meth_mut / total_mut) %>% 
  drop_na()

# plot -------------------------------------------------------------------------
df_plot <- count_m6A_mut

p_015 <- df_plot %>%
  ggplot(aes(x = total_mut, y = meth_mut)) +
  geom_point(color = "black", shape = 21, fill = "grey", alpha = 0.3) +
  geom_smooth(method = 'lm', formula = y ~ x, linetype = "dashed", color = "black", se = FALSE, linewidth = 0.3) +
  stat_cor(method = "spearman", p.accuracy = 0.001) +
  theme_bw() +
  theme(text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black")) +
  xlab("Total mutations") +
  ylab("Mutations in m6A sites")

p_015

# export -----------------------------------------------------------------------
# define plot name as name of script
plot_name <- rstudioapi::getSourceEditorContext()$path %>%
  basename() %>% 
  str_replace(".R", ".png")

ggsave(filename = paste0("20230322_publication_figs/plots/", plot_name), width = 12, height = 6)

# cleanup ----------------------------------------------------------------------
rm(count_m6A_mut, df_plot, tcga_mc3_df_GRCh38, tcga_mc3_df_GRCh38_filtered, plot_name)
