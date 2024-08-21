# ______________________________________________________________________________
# Plot number of m6A sites per data set
# ______________________________________________________________________________

# clear workspace ----
#rm(list = ls())
#gc()

# load libraries ----
options(repos = "https://cran.rstudio.com")
if (!require("pacman")) install.packages("pacman", dependencies = TRUE)
library(pacman)

p_load(tidyverse, readxl, ggsci)

# load data ----
total_m6A_sites <- read_xlsx("20230322_publication_figs/processed_data/001_total_m6A_sites.xlsx")

# analysis ----
# count total unique sites
total_m6A_unique <- total_m6A_sites %>%
  distinct(`Mutation Genome Position`) %>%
  nrow()

# count sites per data set
m6A_sites_per_data_set <- total_m6A_sites %>%
  group_by(data_set) %>%
  summarize(n = n()) %>%
  add_row(data_set = "total distinct sites", n = total_m6A_unique)


# set colors for bars
m6A_sites_per_data_set$color <- "grey"
m6A_sites_per_data_set$color[m6A_sites_per_data_set$data_set == "total distinct sites"] <- "#EFC000FF"

# plot
p_002 <- m6A_sites_per_data_set %>%
  ggplot(aes(x = n, y = reorder(data_set, n))) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8, fill = m6A_sites_per_data_set$color) +
  scale_fill_jco() +
  theme_bw() +
  xlab("Number of m6A sites") +
  ylab("") +
  scale_x_continuous(labels = scales::comma) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

p_002

# cleanup
rm(m6A_sites_per_data_set, total_m6A_sites, total_m6A_unique)

# export plot ----
ggsave(filename = "20230322_publication_figs/plots/002_plot_n_m6A_sites.png", height = 6, width = 12)
