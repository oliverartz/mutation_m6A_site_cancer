## -----------------------------------------------------------------------------
## Purpose of script: Plot m6A sites per gene for entire data set
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
if (!require("pacman")) install.packages("pacman", dependencies = TRUE)
library(pacman)

p_load(tidyverse, readxl, TxDb.Hsapiens.UCSC.hg38.knownGene, ChIPseeker, cowplot, org.Hs.eg.db)

# load data --------------------------------------------------------------------
# m6A sites
total_m6A_sites <- read_xlsx("20230322_publication_figs/processed_data/001_total_m6A_sites.xlsx")

# mutations
load("20230322_publication_figs/processed_data/001_tcga_mc3_df_GRCh38.RData")

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
  starts.in.df.are.0based = FALSE
)

peakAnno <- annotatePeak(total_m6A_sites_gr,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db",
  genomicAnnotationPriority = c("3UTR", "Exon", "5UTR", "Intron", "Promoter", "Downstream", "Intergenic")
)

total_m6A_sites_annotated <- peakAnno %>% as.data.frame()

# analysis ---------------------------------------------------------------------
# find distinct sites
distinct_m6A_sites_anno <- total_m6A_sites_annotated %>%
  mutate(locus = paste0(seqnames, ":", start)) %>%
  distinct(locus, .keep_all = TRUE)

# count genes
count_sites_per_gene <- distinct_m6A_sites_anno %>%
  group_by(SYMBOL) %>%
  tally() %>%
  arrange(desc(n))

median(count_sites_per_gene$n)

# count number of total genes
total_genes <- distinct_m6A_sites_anno$SYMBOL %>%
  unique() %>%
  length()

# summary stats
summary_stats <- summary(count_sites_per_gene$n)

# plot -------------------------------------------------------------------------
p_hist <- count_sites_per_gene %>%
  ggplot(aes(y = n)) +
  geom_histogram(binwidth = 5) +
  theme_bw() +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")) +
  labs(
#    title = "",
#    subtitle = paste0("Total number of genes: ", total_genes, ", mean m6A sites: ", summary_stats[4] %>% round(0)),
    x = "Number of genes",
    y = "Number of m6A sites"
  )

p_box <- count_sites_per_gene %>%
  ggplot(aes(x = 1, y = n)) +
  geom_boxplot() +
  theme_bw() +
  theme(
    axis.text.x = element_text(color = "white"),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(
    x = " ",
    y = "Number of m6A sites"
  )

# combine plots
p_021 <- plot_grid(p_hist, p_box, ncol = 2, rel_widths = c(3, 1))

# export -----------------------------------------------------------------------
# define plot name as name of script
plot_name <- rstudioapi::getSourceEditorContext()$path %>%
  basename() %>%
  str_replace(".R", ".png")

ggsave(
  filename = paste0("20230322_publication_figs/plots/", plot_name),
  height = 6
)

# cleanup ----------------------------------------------------------------------
rm(total_m6A_sites, txdb, total_m6A_sites_mod, total_m6A_sites_gr, peakAnno, total_m6A_sites_annotated, distinct_m6A_sites_anno, count_sites_per_gene, total_genes, summary_stats, p_box, p_hist, plot_name, tcga_mc3_df_GRCh38, ChIPseekerEnv)
