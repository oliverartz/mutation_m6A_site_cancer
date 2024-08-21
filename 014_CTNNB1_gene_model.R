# ______________________________________________________________________________
# Plot CTNNB1 gene model with mutations
# ______________________________________________________________________________

# clear workspace --------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, ggpubr, cowplot, data.table, TxDb.Hsapiens.UCSC.hg38.knownGene, BSgenome.Hsapiens.UCSC.hg38, org.Hs.eg.db, ggsci, trackViewer)

# load data --------------------------------------------------------------------
mutations_in_m6A_motifs <- readxl::read_xlsx("20230322_publication_figs/processed_data/001.1_mutations_in_m6A_motifs.xlsx")
load("20230322_publication_figs/processed_data/001_tcga_mc3_df_GRCh38.RData")

# analysis ---------------------------------------------------------------------


# plot -------------------------------------------------------------------------
# plot gene model ----
# define ENTREZ ID
entrez_id <- 1499

# load transcript
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# load transcript features
transcript_features <- genes(txdb)
transcript_features <- transcript_features[transcript_features$gene_id == entrez_id, ]

trs <- geneModelFromTxdb(txdb,
                         org.Hs.eg.db,
                         chrom = seqnames(transcript_features),
                         start = start(transcript_features),
                         end = end(transcript_features),
                         strand = NULL)


features <- paste0(seqnames(trs[[2]]$dat),":", start(trs[[2]]$dat),"-", end(trs[[2]]$dat))

features_GRanges <- GRanges(features)

features_GRanges <- c(range(trs[[2]]$dat))
features_GRanges <- trs[[2]]$dat
names(features_GRanges) <- trs[[2]]$dat %>% as.data.frame() %>% pull(feature)

features_GRanges_mod <- features_GRanges
names(features_GRanges_mod) <- seq(1, length(features_GRanges_mod), 1)

features_df <- features_GRanges_mod %>% as.data.frame()

# define xlims
min_x <- min(features_df$start)
max_x <- max(features_df$end)

gene_model <- ggplot() +
  geom_hline(yintercept = 20, linewidth = 4, color = "grey") +
  geom_rect(data = features_df, 
            aes(xmin = start, xmax = end, ymin = 0, ymax = 40, fill = feature)) +
  xlim(min_x, max_x) +
  # guides(fill="none") +
  ylab(" ") +
  xlab("Position on chr3") +
  labs(subtitle = "Gene model",
       fill = "Feature") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(color = "white"),
        axis.text.x = element_text(color = "black"),
        panel.ontop = element_blank(),
        legend.position = "bottom",
        panel.border = element_blank(),
        axis.line.x.bottom = element_line(color = 'black'),
        plot.margin = unit(c(0, 5.5, 5.5, 5.5), "pt")) +
  scale_fill_manual(values = c(pal_jco("default")(10)[3], 
                               pal_jco("default")(10)[1], 
                               pal_jco("default")(10)[2]), 
                    breaks = c("utr5", "CDS", "utr3"), 
                    labels = c("5'UTR", "CDS", "3'UTR"))

# plot mutation pileup ----
df_plot <- tcga_mc3_df_GRCh38 %>% 
  filter(Variant_Classification %in% c("3'UTR", "Missense_Mutation", "Silent")) %>% 
  filter(Hugo_Symbol == "CTNNB1") %>% 
  mutate(in_methylated_motif = case_when(in_methylated_motif == 0 ~ "no",
                                         in_methylated_motif == 1 ~ "yes"))

p_mutations <- df_plot %>% 
  ggplot(aes(x = start, color = in_methylated_motif, fill = in_methylated_motif)) +
  geom_histogram(position = "identity", binwidth = 1) +
  xlim(min_x, max_x) +
  # guides(color = "none") +
  xlab("") +
  theme_bw() +
  labs(title = "CTNNB1", 
       subtitle = "Mutations",
       fill = "Site methylated",
       color = "Site methylated",
       y = "Count") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.15,0.7),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt"),
        text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black")) +
  scale_color_jco()  +
  scale_fill_jco()

# combine plots ----
p_014 <- plot_grid(p_mutations, gene_model, ncol = 1,
                 rel_heights = c(1.5, 1)
)

p_014

# export -----------------------------------------------------------------------
# define plot name as name of script
plot_name <- rstudioapi::getSourceEditorContext()$path %>%
  basename() %>% 
  str_replace(".R", ".png")

ggsave(filename = paste0("20230322_publication_figs/plots/", plot_name), width = 12, height = 6)

# cleanup ----------------------------------------------------------------------
rm(df_plot, features_df, features_GRanges, features_GRanges_mod, gene_model, 
   mutations_in_m6A_motifs, p_mutations, tcga_mc3_df_GRCh38, transcript_features, 
   trs, entrez_id, max_x, min_x, plot_name, txdb, features)
