# ______________________________________________________________________________
# Plot DRACH motif from HOMER analysis
# ______________________________________________________________________________

# clear workspace --------------------------------------------------------------
#rm(list = ls())
#gc()

# load packages ----------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, ggseqlogo, data.table)

# load data --------------------------------------------------------------------
motif <- fread("/Users/artzo/Documents/4_projects/20220623_synonymous_m6A/20230810_HOMER/HOMER_output/homerResults/motif1.motif") %>% 
  t() %>% 
  as.data.frame() %>% 
  drop_na()

# plot -------------------------------------------------------------------------
rownames(motif) <- c("A", "C", "G", "U")
colnames(motif) <- seq(1, ncol(motif), 1)

# select DRACH
motif <- motif[, 3:7]

p_020 <- ggplot() +
  annotate('rect', xmin = 2.5, xmax = 3.5, ymin = -0.05, ymax = 1, alpha = .1, color = 'black', fill = '#EFC000FF') +
  geom_logo(as.matrix(motif), method = "prob") +
  theme_logo() +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("-2", "-1", "m6A", "+1", "+2")) +
  theme(axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black")) +
  labs(subtitle = paste0("p-value: 1e-7847"))

# export -----------------------------------------------------------------------
# define plot name as name of script
plot_name <- rstudioapi::getSourceEditorContext()$path %>%
  basename() %>% 
  str_replace(".R", ".png")

ggsave(filename = paste0("20230322_publication_figs/plots/", plot_name), height = 6)

# cleanup
rm(motif, plot_name)
