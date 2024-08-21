# ______________________________________________________________________________
# Annotate motifs
# ______________________________________________________________________________

# clear work space -------------------------------------------------------------
rm(list = ls())
gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, readxl, writexl, BSgenome, BSgenome.Hsapiens.NCBI.GRCh38)

# load data --------------------------------------------------------------------
total_m6A_unique <- read_xlsx("20230322_publication_figs/processed_data/001_total_m6A_unique.xlsx")
load("20230322_publication_figs/processed_data/001_tcga_mc3_df_GRCh38.RData")

# Mutations in m6A motifs ------------------------------------------------------
## The same mutation might affect multiple motifs. This is represented here.
### make GRanges objects with motif positions as ranges
total_m6A_unique_GRanges <- makeGRangesFromDataFrame(total_m6A_unique,
  keep.extra.columns = TRUE,
  ignore.strand = FALSE,
  seqinfo = NULL,
  start.field = "start_motif",
  end.field = "end_motif",
  strand.field = "Strand",
  starts.in.df.are.0based = FALSE
)

tcga_mc3_GRanges <- makeGRangesFromDataFrame(tcga_mc3_df_GRCh38,
  keep.extra.columns = TRUE,
  ignore.strand = FALSE,
  seqinfo = NULL,
  start.field = "start",
  end.field = "end",
  strand.field = "STRAND",
  starts.in.df.are.0based = FALSE
)

# find overlaps
idx <- findOverlaps(tcga_mc3_GRanges, total_m6A_unique_GRanges, ignore.strand = FALSE) %>% as.data.frame()

# annotate mutations with specific motif they disrupt
# if a mutation is in multiple motifs, the row will be duplicated and annotated
motif_positions <- paste0(total_m6A_unique$chr, ":", total_m6A_unique$start_motif, total_m6A_unique$Strand)
mutations_in_m6A_motifs <- tcga_mc3_df_GRCh38[idx$queryHits, ] %>% 
  mutate(motif_position = motif_positions[idx$subjectHits],
         start_motif = total_m6A_unique$start_motif[idx$subjectHits]) 

# How many motifs are disrupted by a given mutation?
n_overlaps <- countOverlaps(tcga_mc3_GRanges, total_m6A_unique_GRanges) %>% as.data.frame()
# There are 405 mutations that overlap with multiple motifs


### Annotate effect of mutation on DRACH motif
# add WT motif sequence
mutations_in_m6A_motifs <- mutations_in_m6A_motifs %>%
  mutate(motif_wt = getSeq(BSgenome.Hsapiens.NCBI.GRCh38,
    seqnames,
    start = start_motif,
    end = start_motif + 4,
    strand = STRAND,
    as.character = TRUE
  ))


# position in motif that is mutated
# Important: Mutated position has to be calculated differently for + and - Strands
mutations_in_m6A_motifs <- mutations_in_m6A_motifs %>% 
  mutate(mut_pos = case_when(mutations_in_m6A_motifs$STRAND == "+" ~ mutations_in_m6A_motifs$start - mutations_in_m6A_motifs$start_motif + 1,
                             mutations_in_m6A_motifs$STRAND == "-" ~ 5 - (mutations_in_m6A_motifs$start - mutations_in_m6A_motifs$start_motif)))

# add MUT motif sequence
alt_allele <- str_sub(mutations_in_m6A_motifs$HGVSc, -1)
mutations_in_m6A_motifs <- mutations_in_m6A_motifs %>% mutate(motif_mut = motif_wt)

substr(mutations_in_m6A_motifs$motif_mut, mutations_in_m6A_motifs$mut_pos, mutations_in_m6A_motifs$mut_pos) <- alt_allele

# find DRACH in WT
## define DRACH
DRACH <- "(A|G|T)(A|G)AC(A|C|T)"

mutations_in_m6A_motifs <- mutations_in_m6A_motifs %>% 
  mutate(motif_wt_DRACH = if_else(grepl(DRACH, motif_wt), 1, 0))

# Careful: There are cases in which a mutation overlaps with a methylated site and a DRACH motif, but the methylated site is not within a DRACH motif

# find DRACH in MUT
mutations_in_m6A_motifs <- mutations_in_m6A_motifs %>% 
  mutate(motif_mut_DRACH = if_else(grepl(DRACH, motif_mut), 1, 0))

# annotate effect of mutation
mutations_in_m6A_motifs <- mutations_in_m6A_motifs %>% 
  mutate(motif_consequence = case_when(
    motif_wt_DRACH == 1 & motif_mut_DRACH == 1 ~ "neutral",
    motif_wt_DRACH == 0 & motif_mut_DRACH == 0 ~ "neutral",
    motif_wt_DRACH == 1 & motif_mut_DRACH == 0 ~ "loss",
    motif_wt_DRACH == 0 & motif_mut_DRACH == 1 ~ "gain"
  ))

# export data ------------------------------------------------------------------
write_xlsx(mutations_in_m6A_motifs,"20230322_publication_figs/processed_data/001.1_mutations_in_m6A_motifs.xlsx")

m6A_motifs <- read_xlsx("20230322_publication_figs/processed_data/001.1_mutations_in_m6A_motifs.xlsx")
