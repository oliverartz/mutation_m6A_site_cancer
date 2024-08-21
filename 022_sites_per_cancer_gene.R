## -----------------------------------------------------------------------------
## Purpose of script: Plot number of m6A sites per cancer gene (OG, TS)
##
## Author: Oliver Artz
## Date Created:
## Email: artzo@mskcc.org
##
## -----------------------------------------------------------------------------

# clear work space -------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, data.table, readxl, ChIPseeker, TxDb.Hsapiens.UCSC.hg38.knownGene, cowplot, ggsci, ggpubr, biomaRt)

# load data --------------------------------------------------------------------
# m6A sites
total_m6A_sites <- read_xlsx("20230322_publication_figs/processed_data/001_total_m6A_sites.xlsx")

# mutations
load("20230322_publication_figs/processed_data/001_tcga_mc3_df_GRCh38.RData")

# oncoKB oncogenes / tumor suppressors
oncokb_gene_list <- read_delim("20221025_coexpression_machinery/data/cancerGeneList.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)

# wrangle ---------------------------------------------------------------------
# annotate loci
# load txdb
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# add 'chr' to seqnames
total_m6A_sites_mod <- total_m6A_sites %>% 
  mutate(chr = paste0("chr", chr))

# make GRanges object from m6A sites
total_m6A_sites_gr <- makeGRangesFromDataFrame(total_m6A_sites_mod,
                                               keep.extra.columns = TRUE,
                                               ignore.strand = FALSE,
                                               seqinfo = NULL,
                                               start.field = "Start",
                                               end.field = "End",
                                               strand.field = "Strand",
                                               starts.in.df.are.0based = FALSE)

peakAnno <- annotatePeak(total_m6A_sites_gr,
                         TxDb = txdb, 
                         annoDb = "org.Hs.eg.db",
                         genomicAnnotationPriority = c("3UTR", "Exon", "5UTR", "Intron", "Promoter", "Downstream", "Intergenic"))

total_m6A_sites_annotated <- peakAnno %>% as.data.frame()

# find distinct sites
distinct_m6A_sites_anno <- total_m6A_sites_annotated %>% 
  mutate(locus = paste0(seqnames,":", start)) %>% 
  distinct(locus, .keep_all = TRUE)

## annotate cancer genes ----
## add entrezgene_id information to OncoKB genes -------------------------------
# extract oncogene and tumor suppressor information
oncokb_df <- oncokb_gene_list %>%
  dplyr::select(`Hugo Symbol`, `Is Oncogene`, `Is Tumor Suppressor Gene`) %>%
  mutate(cancer_gene = ifelse(`Is Oncogene` == "Yes" & `Is Tumor Suppressor Gene` == "Yes", "OG & TS",
                              ifelse(`Is Oncogene` == "No" & `Is Tumor Suppressor Gene` == "No", "OnkoKB, neither",
                                     ifelse(`Is Oncogene` == "Yes", "OG",
                                            ifelse(`Is Tumor Suppressor Gene` == "Yes", "TS", "not defined")
                                     ))))

idx <- match(distinct_m6A_sites_anno$SYMBOL, oncokb_df$`Hugo Symbol`)
distinct_m6A_sites_anno$cancer_gene <- oncokb_df$cancer_gene[idx]
distinct_m6A_sites_anno$cancer_gene[is.na(distinct_m6A_sites_anno$cancer_gene)] <- "not in OncoKB"

# analysis ---------------------------------------------------------------------
# count genes
count_sites_per_gene <- distinct_m6A_sites_anno %>% 
  group_by(SYMBOL, cancer_gene) %>% 
  tally() %>% 
  arrange(desc(n))

count_sites_per_gene %>% 
  group_by(cancer_gene) %>% 
  summarize(median_n = median(n))

# plot -------------------------------------------------------------------------
df_plot <- count_sites_per_gene %>%
  filter(cancer_gene %in% c("OG", "TS", "OG & TS")) 

# calculate mean for each cancer_gene category
summary_stats <- df_plot %>% 
  group_by(cancer_gene) %>% 
  summarize(avg_number = mean(n))

# Specify the desired order of levels for cancer_gene
gene_order <- c("OG", "TS", "OG & TS")

# Reorder the cancer_gene column as a factor with the desired order
df_plot$cancer_gene <- factor(df_plot$cancer_gene, levels = gene_order)

# histogram
p_hist <- df_plot %>% 
  ggplot(aes(y = n)) +
  geom_histogram(binwidth = 10, aes(fill = cancer_gene), alpha = 0.8, color = "black") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Number of m6A sites per gene",
       #       subtitle = "",
       x = "Number of genes",
       y = "Number of m6A sites") +
  scale_fill_jco() +
  ylim(0, 400) +
  facet_wrap(~ cancer_gene) +
  theme(text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"))

# specify comparisons for boxplot
my_comparisons <- list( c("OG", "TS"), c("TS", "OG & TS"), c("OG", "OG & TS") )

# boxplot
p_box <- df_plot %>% 
  ggplot(aes(x = cancer_gene, y = n, fill = cancer_gene)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(shape = 21, width = 0.3, alpha = 0.5) +
  stat_compare_means(comparisons = my_comparisons) +
  scale_fill_jco() +
  theme_bw() +
  theme(legend.position = "none") +
  ylim(0, 400) +
  labs(title = "",
       x = "",
       y = "Number of m6A sites") +
  theme(text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"))

# combine plots
p_022 <- plot_grid(p_hist, p_box, ncol = 2, rel_widths = c(1,1))

p_022

# export -----------------------------------------------------------------------
ggsave(filename = paste0("20230322_publication_figs/plots/022_sites_per_cancer_gene.png"), 
      height = 6)

# cleanup ----------------------------------------------------------------------
rm(txdb, total_m6A_sites_mod, total_m6A_sites_gr, peakAnno, oncokb_gene_list, 
   idx, oncokb_df, count_sites_per_gene, df_plot, distinct_m6A_sites_anno, my_comparisons, 
   p_box, p_hist, summary_stats, tcga_mc3_df_GRCh38, total_m6A_sites, total_m6A_sites_annotated,
   ChIPseekerEnv, gene_order)
