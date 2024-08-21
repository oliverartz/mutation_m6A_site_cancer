# ______________________________________________________________________________
# Plot DESeq2 results for genes that have loss of m6A motif
# Gene level
# ______________________________________________________________________________

# clear workspace --------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, ggsci, ggpubr, cowplot, ggrepel)

# set parameters ---------------------------------------------------------------
threshold_recurrence <- 2

# load data --------------------------------------------------------------------
# DESeq2 analysis was peformed on HPC
# Results were stored in respective folder
files <- list.files(path = "20230519_expression_HPC/processed_data/genelevel", pattern = "*_expression_DEseq2.csv", full.names = TRUE)

# Function to read in a file and add a column with the file number
read_file <- function(file) {
  data <- read_csv(file,show_col_types = FALSE)
  file_num <- str_extract(basename(file), "\\d+") %>% as.numeric()
  data <- data %>% mutate(gene_row = file_num)
  return(data)
}

# Read all files into a list of data frames
data_list <- map(files, read_file)

# Bind all data frames together
deseq2_results_df <- bind_rows(data_list) %>% arrange(gene_row)

# mutations in m6A motifs ------------------------------------------------------
mutations_in_m6A_motifs <- readxl::read_xlsx("20230322_publication_figs/processed_data/001.1_mutations_in_m6A_motifs.xlsx")

# wrangle ----------------------------------------------------------------------
# find recurring mutations
recurring_genes <- mutations_in_m6A_motifs %>% 
  filter(Variant_Classification %in% c("Missense_Mutation", "Silent", "3'UTR")) %>% 
  filter(motif_consequence == "loss") %>% 
  filter(in_methylated_motif == 1) %>% 
  filter(in_DRACH == 1) %>%
  group_by(cancer_type, Hugo_Symbol) %>% 
  summarize(n = n()) %>% 
  filter(n >= threshold_recurrence) %>% 
  arrange(desc(n)) %>% 
  mutate(cancer_type_gene = paste0(cancer_type, "_", Hugo_Symbol))

# filter results for recurring threshold
deseq2_results_df <- deseq2_results_df %>% 
  mutate(cancer_type_gene = paste0(cancer_type, "_", Hugo_Symbol)) %>% 
  filter(cancer_type_gene %in% recurring_genes$cancer_type_gene)

# plot -------------------------------------------------------------------------
# define thresholds
psig_threshold <- 0.05
fc_threshold <- 2

# make volcano plot
p_012 <- deseq2_results_df %>% 
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
  theme(text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black")) +
  scale_fill_jco() +
  labs(title = "DGE mutated vs wildtype",
       subtitle = "Gene-level",
       fill = "Cancer type")

p_012

# export -----------------------------------------------------------------------
# define plot name as name of script
plot_name <- rstudioapi::getSourceEditorContext()$path %>%
  basename() %>% 
  str_replace(".R", ".png")

ggsave(filename = paste0("20230322_publication_figs/plots/", plot_name), width = 12, height = 6)

# cleanup ----------------------------------------------------------------------
rm(plot_name, data_list, deseq2_results_df, fc_threshold, files, plot_name, psig_threshold, read_file)
