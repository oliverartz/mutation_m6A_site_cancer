# ______________________________________________________________________________
# Plot fraction of total mutations that are in methylated sites
# for each cancer type
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

# filter patients with more than 10 total mutations
 count_m6A_mut <- count_m6A_mut %>% 
   filter(total_mut > 5)

# filter patients with more than 0 mutation in m6A sites
# count_m6A_mut <- count_m6A_mut %>% 
#   filter(meth_mut > 0)

# Calculate quartiles for each cancer_type
quantiles <- count_m6A_mut %>%
  group_by(cancer_type) %>%
  summarize(q1 = quantile(methylated_total, 0.25),
            median = quantile(methylated_total, 0.5),
            q3 = quantile(methylated_total, 0.75))

# Add quantile information to count_m6A_mut
# below 1st quartile, above third quartile
count_m6A_mut <- count_m6A_mut %>%
  left_join(quantiles, by = "cancer_type") %>%
  mutate(quantile_group = case_when(
    methylated_total <= q1 ~ "Below 1st Quartile",
    methylated_total >= q3 ~ "Above 3rd Quartile",
    TRUE ~ "Within 1st and 3rd Quartile"
  ))
 
# plot -------------------------------------------------------------------------
df_plot <- count_m6A_mut

p_016 <- df_plot %>%
  mutate(cancer_type = str_remove(cancer_type, "TCGA-")) %>% 
#  filter(total_mut > 10) %>% 
  ggplot(aes(x = cancer_type, y = methylated_total)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  geom_jitter(aes(fill = quantile_group), shape = 21, color = "black", alpha = 0.2, width = 0.2) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(color = "black"),
              axis.text = element_text(color = "black"),
              axis.title = element_text(color = "black")) +
  labs(title = "m6A mutation load",
       subtitle = "Patients with >5 total mutations",
       x = "",
       y = "Mut. m6A site / mut. total",
       fill = "Quartile group") +
  scale_fill_jco()

p_016

# export -----------------------------------------------------------------------
# define plot name as name of script
plot_name <- "016_fraction_mutation_in_meth_site_cancer_type.png"

ggsave(filename = paste0("20230322_publication_figs/plots/", plot_name), width = 12, height = 6)

# export count and quartile data
# write_csv(count_m6A_mut, file = "20230322_publication_figs/processed_data/m6A_load_quartiles.csv")
write_csv(count_m6A_mut, file = "20230322_publication_figs/processed_data/016_m6A_load_quartiles.csv")

# cleanup ----------------------------------------------------------------------
rm(count_m6A_mut, df_plot, quantiles, tcga_mc3_df_GRCh38, tcga_mc3_df_GRCh38_filtered, plot_name)
