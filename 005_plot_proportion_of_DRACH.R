# ______________________________________________________________________________
# Plot proportion of m6A deposition in DRACH vs outside DRACH
# independent of mutations, just all m6A sites from data set
# ______________________________________________________________________________

# clear workspace --------------------------------------------------------------
#rm(list = ls())
#gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, readxl, ggsci, cowplot)

# load data --------------------------------------------------------------------
total_m6A_sites <- read_xlsx("20230322_publication_figs/processed_data/001_total_m6A_sites.xlsx")
total_m6A_unique <- read_xlsx("20230322_publication_figs/processed_data/001_total_m6A_unique.xlsx")

# wrangle ----------------------------------------------------------------------
# count per data set
DRACH_count <- total_m6A_sites %>% 
  group_by(data_set, DRACH) %>% 
  tally()

# count for unique sites
unique_DRACH_count <- total_m6A_unique %>% 
  group_by(DRACH) %>% 
  tally() %>% 
  mutate(data_set = "total distinct sites") %>% 
  dplyr::select(data_set, DRACH, n)

# add to df
DRACH_count <- rbind(DRACH_count, unique_DRACH_count)

# refactor fill
DRACH_count <- DRACH_count %>% 
  mutate(DRACH = case_when(DRACH == 0 ~ "outside DRACH",
                           DRACH == 1 ~ "in DRACH"))

# plot -------------------------------------------------------------------------
# determine order of y axis
order_y <- total_m6A_sites %>% group_by(data_set) %>% 
  tally() %>% 
  ungroup() %>% 
  add_row(data_set = "total distinct sites", n = nrow(total_m6A_unique)) %>% 
  arrange(n) %>% 
  mutate(data_set = as.factor(data_set))

DRACH_count$data_set <- factor(DRACH_count$data_set, levels = order_y$data_set)
DRACH_count$DRACH <- factor(DRACH_count$DRACH)

# plot
p_005 <- DRACH_count %>%
  ggplot(aes(x = n, y = data_set, fill = DRACH)) +
  geom_bar(stat = "identity", position = "fill", alpha = 0.8, color = "black") +
  scale_fill_jco() +
  theme_bw() +
  theme(text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black")) +
  ylab("") +
  xlab("Proportion of m6A sites")

p_005

# export plot
ggsave(filename = "20230322_publication_figs/plots/005_plot_proportion_of_DRACH.png", height = 6, width = 12)

# cleanup
rm(DRACH_count, order_y, total_m6A_sites, total_m6A_unique, unique_DRACH_count)

