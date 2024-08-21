# ______________________________________________________________________________
# Plot difference in total mutations between quartiles for significant cancer types
# ______________________________________________________________________________

# clear workspace --------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, data.table, RTCGA.clinical, survival, survminer)

# load data --------------------------------------------------------------------
# m6A mutational load with quartiles
count_m6A_mut <- read_csv("20230322_publication_figs/processed_data/017_m6A_load_quartiles_survival.csv", show_col_types = FALSE)

# analysis ---------------------------------------------------------------------
# define data set
cancer_data_set <- "TCGA-THCA"

# make data frame for plotting
df_plot <- count_m6A_mut %>% 
  filter(cancer_type == cancer_data_set) %>% 
  filter(quantile_group %in% c("Below 1st Quartile", "Above 3rd Quartile"))


# plot -------------------------------------------------------------------------
# p_019 <- df_plot %>% 
#   ggplot(aes(x = quantile_group, y = total_mut, fill = quantile_group)) +
#   geom_boxplot(alpha = 0.8, outlier.shape = NA) +
#   geom_jitter(shape = 21, color = "black", alpha = 0.3, width = 0.3) +
#   stat_compare_means(
#     aes(label = paste0("p = ", format(after_stat(p.value), scientific = TRUE))),
#     method = "wilcox.test"
#   ) +
#   theme_bw() +
#   theme(legend.position = "none") +
#   labs(title = "Comparison of total mutational load",
#        x = "",
#        y = "total mutations") +
#   scale_fill_jco()

# Manual calculation of Wilcoxon test
test_result <- wilcox.test(total_mut ~ quantile_group, data = df_plot)
p_val <- format(signif(test_result$p.value, digits = 2), scientific = TRUE)

# Build the plot
p_019 <- df_plot %>%
  ggplot(aes(x = quantile_group, y = total_mut, fill = quantile_group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(shape = 21, color = "black", alpha = 0.3, width = 0.3) +
  annotate("text", x = 1, y = max(df_plot$total_mut, na.rm = TRUE), label = paste0("p = ", p_val)) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(color = "black"),
              axis.text = element_text(color = "black"),
              axis.title = element_text(color = "black")) +
  labs(title = "Mutational load THCA",
       x = "",
       y = "Total mutations") +
  scale_fill_jco() +
  scale_y_log10()

p_019

# export -----------------------------------------------------------------------
# define plot name as name of script
plot_name <- rstudioapi::getSourceEditorContext()$path %>%
  basename() %>% 
  str_replace(".R", "")

ggsave(filename = paste0("20230322_publication_figs/plots/", plot_name, "_", cancer_data_set, ".png"), width = 6, height = 6)

# cleanup ----------------------------------------------------------------------
rm(count_m6A_mut, df_plot, cancer_data_set, plot_name)

