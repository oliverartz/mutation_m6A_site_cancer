## -----------------------------------------------------------------------------
## Purpose of script: Analyze motif-based DGE results
##
## Author: Oliver Artz
## Date Created: Jun 18, 2024
## Email: artzo@mskcc.org
##
## -----------------------------------------------------------------------------

# clear work space -------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, data.table, ggrepel)

# set parameters ---------------------------------------------------------------
padj_val_threshold <- 0.05
logFC_threshold <- 2

# load data --------------------------------------------------------------------
## motif data
load("20240403_gained_motifs/processed_data/001_tcga_mc3_motifs_anno.RData")

## input data frame for DESeq2
load("20240403_gained_motifs/processed_data/008_patient_data_df_motif_level.RData")

## recurring gained motifs (patient_data_motiflevel_df)
load("20240403_gained_motifs/processed_data/008_patient_data_df_motif_level_all_data.RData")

## load motif-based DGE results
filenames <- list.files("20240403_gained_motifs/processed_data/motif_based", pattern = "*.csv", full.names = TRUE)

read_file <- function(file) {
  data <- fread(file)
  file_num <- str_extract(basename(file), "\\d+") %>% as.numeric()
  data <- data %>% 
    mutate(motif_num = file_num,
           motif = as.list(motif))
  return(data)
}

# Read all files into a list of data frames
data_list <- map(filenames, read_file)

# Bind all data frames together
deseq2_results_df <- bind_rows(data_list) %>% arrange(motif_num)

# wrangle ----------------------------------------------------------------------
## merge with patient data
deseq2_results_patients <- patient_data_df_temp %>% 
  left_join(deseq2_results_df %>% 
              dplyr::select(-Hugo_Symbol), 
            by = c("motif", "cancer_type"))

# QC ---------------------------------------------------------------------------
## Check if all recurring motifs have expression data
df_qc <- deseq2_results_patients %>% filter(is.na(log2FoldChange))

df_qc

# plot -------------------------------------------------------------------------
n_significant_motifs <- deseq2_results_patients %>% 
  filter(padj < padj_val_threshold & abs(log2FoldChange) > logFC_threshold) %>% 
  nrow()

n_total_motifs <- deseq2_results_patients %>% 
  nrow()

## make a volcano plot
p_011_volcano_motif <- deseq2_results_patients %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(shape = 21, alpha = 0.8, fill = "grey") +
  geom_point(data = deseq2_results_patients %>% 
               filter(padj < padj_val_threshold & abs(log2FoldChange) > logFC_threshold), 
             aes(x = log2FoldChange, 
                 y = -log10(padj),
                 fill = cancer_type),
             shape = 21) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(data = deseq2_results_patients %>% 
                    filter(padj < padj_val_threshold & abs(log2FoldChange) > logFC_threshold),  
                   aes(label = paste0(Hugo_Symbol, "\n", motif)),
                   box.padding = 0.5,
                   point.padding = 0.5,
                   segment.color = "grey50",
                   segment.size = 0.5,
                   segment.alpha = 0.5,
                   nudge_x = 0.5,
                   nudge_y = 0.5) +
  theme_bw() +
  theme(text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black")) +
  labs(title = "Volcano plot of motif-based DGE results",
       subtitle = paste0("padj < ", padj_val_threshold, " & abs(log2FC) > ", logFC_threshold, " (n = ", n_significant_motifs,"/", n_total_motifs, ")"),
       x = "higher in mut <- log2 fold change -> higher in WT", 
       y = "-log10 adjusted p-value", 
       fill = "Cancer type")

p_011_volcano_motif

# export -----------------------------------------------------------------------
ggsave("20240403_gained_motifs/plots/011_volcano_plot_motif.png", width = 12, height = 6)

# cleanup ----------------------------------------------------------------------
rm(data_list, deseq2_results_df, deseq2_results_patients, df_qc, patient_data_df_temp, patient_data_motiflevel_df, tcga_mc3_df_GRCh38,
   filenames, logFC_threshold, n_significant_motifs, n_total_motifs, padj_val_threshold, read_file)
