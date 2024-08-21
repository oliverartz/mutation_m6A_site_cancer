## -----------------------------------------------------------------------------
## Purpose of script: Analyze motif-based DGE results
##
## Author: Oliver Artz
## Date Created: Jun 20, 2024
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
load("20240403_gained_motifs/processed_data/012_patient_data_df_gene_level.RData")

## recurring gained motifs (patient_data_motiflevel_df)
load("20240403_gained_motifs/processed_data/012_patient_data_df_gene_level_all_data.RData")

## load motif-based DGE results
filenames <- list.files("20240403_gained_motifs/processed_data/gene_based", pattern = "*.csv", full.names = TRUE)

read_file <- function(file) {
  data <- fread(file)
  file_num <- str_extract(basename(file), "\\d+") %>% as.numeric()
  data <- data %>% 
    mutate(motif_num = file_num)
  return(data)
}

# Read all files into a list of data frames
data_list <- map(filenames, read_file)

# Bind all data frames together
deseq2_results_df <- bind_rows(data_list) %>% 
  arrange(motif_num)

# wrangle ----------------------------------------------------------------------
## merge with patient data
deseq2_results_patients <- patient_data_df_temp %>% 
  mutate(Hugo_Symbol = as.character(Hugo_Symbol)) %>% 
  left_join(deseq2_results_df,
            by = c("Hugo_Symbol", "cancer_type"))

# remove 'TCGA' from cancer type
deseq2_results_patients$cancer_type <- deseq2_results_patients$cancer_type %>% 
  str_remove("TCGA-")

# QC ---------------------------------------------------------------------------
## Check if all recurring motifs have expression data
df_qc <- deseq2_results_patients %>% filter(is.na(log2FoldChange))

df_qc$motif_num

# plot -------------------------------------------------------------------------
n_significant_genes <- deseq2_results_patients %>% 
  filter(padj < padj_val_threshold & abs(log2FoldChange) > logFC_threshold) %>% 
  nrow()

n_total_genes <- deseq2_results_patients %>% 
  nrow()

## make a volcano plot
p015_deseq2_gene <- deseq2_results_patients %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(shape = 21, alpha = 0.5, fill = "grey") +
  geom_point(data = deseq2_results_patients %>% 
               filter(padj < padj_val_threshold & abs(log2FoldChange) > logFC_threshold), 
             aes(x = log2FoldChange, 
                 y = -log10(padj),
                 fill = cancer_type),
             shape = 21) +
  geom_hline(yintercept = -log10(padj_val_threshold), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(data = deseq2_results_patients %>% 
                    filter(padj < padj_val_threshold & abs(log2FoldChange) > logFC_threshold),  
                   aes(label = paste0(Hugo_Symbol)),
                   box.padding = 0.5,
                   point.padding = 0.5,
                   segment.color = "grey50",
                   segment.size = 0.5,
                   segment.alpha = 0.5,
                   nudge_x = 0.5,
                   nudge_y = 0.5,
                  size = 3) +
  theme_bw() +
  theme(text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black")) +
  labs(title = "Volcano plot of gene-based DGE results",
       subtitle = paste0("padj < ", padj_val_threshold, " & abs(log2FC) > ", logFC_threshold, " (n = ", n_significant_genes,"/", n_total_genes, ")"),
       x = "higher in mut <- log2 fold change -> higher in WT", 
       y = "-log10 adjusted p-value", 
       fill = "Cancer type")

p015_deseq2_gene

# stats for publication text
n_significant_genes <- deseq2_results_patients %>% 
  filter(padj < padj_val_threshold & abs(log2FoldChange) > logFC_threshold) %>% 
  nrow()


# number of genes that are downregulated after gain in mut
deseq2_results_patients %>% filter(log2FoldChange > logFC_threshold,
                                   padj < padj_val_threshold) %>% 
  mutate(id = paste0(cancer_type, "-", Hugo_Symbol)) %>% 
  pull(id) %>% length()

# number of genes that are upregulated after gain in mut
deseq2_results_patients %>% filter(log2FoldChange < -logFC_threshold,
                                   padj < padj_val_threshold) %>% 
  mutate(id = paste0(cancer_type, "-", Hugo_Symbol)) %>% 
  pull(id)



# all significant genes
deseq2_results_patients %>% filter(abs(log2FoldChange) > logFC_threshold,
                                   padj < padj_val_threshold) %>% 
  mutate(id = paste0(cancer_type, "-", Hugo_Symbol)) %>% 
  pull(id)

# export -----------------------------------------------------------------------
ggsave("20240403_gained_motifs/plots/015_volcano_gene_motif.png", width = 12, height = 6)

# analysis ---------------------------------------------------------------------
## get significant genes
significant_genes <- deseq2_results_patients %>% 
  filter(padj < padj_val_threshold & abs(log2FoldChange) > logFC_threshold) %>% 
  select(Hugo_Symbol, cancer_type, log2FoldChange, padj)

## get mutation info for significant genes
get_mutation_info <- function(gene_count){
df_temp <- tcga_mc3_df_GRCh38 %>% 
  filter(motif_class == "gain") %>%
  filter(Hugo_Symbol %in% significant_genes$Hugo_Symbol[gene_count]) %>% 
  filter(cancer_type %in% paste0("TCGA-", significant_genes$cancer_type[gene_count]))

dfSNP_rs <- df_temp$dbSNP_RS %>% 
  unique() %>% 
  paste(collapse = ", ")

variant_classification <- df_temp %>% 
  group_by(Variant_Classification) %>% 
  summarize(n = n()) %>%
  mutate(Hugo_Symbol = significant_genes$Hugo_Symbol[gene_count],
         cancer_type = paste0("TCGA-", significant_genes$cancer_type[gene_count])) %>% 
  pivot_wider(names_from = Variant_Classification, values_from = n, values_fill = 0) %>% 
  mutate(dbSNP_RS = dfSNP_rs)


return(variant_classification)
}

# get_mutation_info for all significant genes
mutation_info_list <- map(1:nrow(significant_genes), get_mutation_info) %>% 
  bind_rows()

# make NA to 0
mutation_info_list[is.na(mutation_info_list)] <- 0

# calculate mutation stats
mutation_info_list <- mutation_info_list %>% 
  mutate(total_mutations = Missense_Mutation + Silent + `3'UTR`) %>% 
  relocate(dbSNP_RS, .after = last_col())

# export
write_csv(mutation_info_list, "20240403_gained_motifs/processed_data/015_mutation_info_list.csv")

# cleanup ----------------------------------------------------------------------
rm(data_list, deseq2_results_df, deseq2_results_patients, df_qc,
  mutation_info_list, patient_data_df_temp, patient_data_genelevel_df,
  significant_genes, tcga_mc3_df_GRCh38, filenames, logFC_threshold, 
  n_significant_genes, n_total_genes, padj_val_threshold, get_mutation_info, 
  read_file)
