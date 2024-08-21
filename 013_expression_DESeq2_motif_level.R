# ______________________________________________________________________________
# Plot DESeq2 results for genes that have loss of m6A motif
# Motif level
# ______________________________________________________________________________

# clear workspace --------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, ggsci, ggpubr, cowplot, ggrepel, data.table)

# load data --------------------------------------------------------------------
# DESeq2 analysis was peformed on HPC
files <- list.files(path = "20230519_expression_HPC/processed_data/motiflevel", pattern = "*_expression_DEseq2.csv", full.names = TRUE)

# function to read in a file and add a column with the file number
read_file <- function(file) {
  data <- read_csv(file,show_col_types = FALSE)
  file_num <- str_extract(basename(file), "\\d+") %>% as.numeric()
  data <- data %>% mutate(gene_row = file_num)
  return(data)
}

# read all files into a list of data frames
data_list <- map(files, read_file)

# bind all data frames together
deseq2_results_df <- bind_rows(data_list) %>% arrange(gene_row)

# recurrent motifs DESeq2
recurrent_motifs <- fread("20230519_expression_HPC/data/recurrent_motifs.csv")

# analysis ---------------------------------------------------------------------


# plot -------------------------------------------------------------------------
# add Hugo Symbol to DESeq2 results
idx <- match(deseq2_results_df$Ensembl_ID, recurrent_motifs$ensembl_id)
deseq2_results_df$Hugo_Symbol <- recurrent_motifs$Hugo_Symbol[idx] 

# define thresholds
psig_threshold <- 0.05
fc_threshold <- 2

p_013 <- deseq2_results_df %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(psig_threshold), linetype = "dashed", alpha = 0.5) +
  geom_point(shape = 21, color = "black", fill = "grey", alpha = 0.5) +
  geom_point(data = filter(deseq2_results_df, -log10(padj) >= -log10(psig_threshold) & abs(log2FoldChange) >= fc_threshold),
             aes(fill = cancer_type), 
             shape = 21, 
             color = "black", alpha = 1) +
  geom_text_repel(data = filter(deseq2_results_df, -log10(padj) >= -log10(psig_threshold) & abs(log2FoldChange) >= fc_threshold),
                  aes(label = Hugo_Symbol),
                  size = 3,
                  max.overlaps = Inf,
                  alpha = 0.8) +
  theme_bw() +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  ) +
  scale_fill_jco() +
  labs(title = "DGE mutated vs wildtype",
       subtitle = "Motif-level",
       fill = "Cancer type")
  
p_013

# export -----------------------------------------------------------------------
# define plot name as name of script
plot_name <- rstudioapi::getSourceEditorContext()$path %>%
  basename() %>% 
  str_replace(".R", ".png")

ggsave(filename = paste0("20230322_publication_figs/plots/", plot_name), width = 12, height = 6)

# cleanup ----------------------------------------------------------------------
rm(data_list, deseq2_results_df, recurrent_motifs, idx, fc_threshold, files, plot_name, psig_threshold, read_file)
