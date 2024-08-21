## -----------------------------------------------------------------------------
## Purpose of script: Find mutations that lead to gain of DRACH
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

p_load(tidyverse, data.table, readxl, maftools, BSgenome, BSgenome.Hsapiens.NCBI.GRCh38, Biostrings)

# load data --------------------------------------------------------------------
# TCGA cancer types ----
# load cancer types
sample_cancer_type <- read_delim("~/Documents/4_projects/20220623_synonymous_m6A/20220623_gathering_m6A_sites/syn_mut_data/TCGA_James/20220801_data/id2type.txt", delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE)

# rename columns
colnames(sample_cancer_type) <- c("patient", "data_set")

# TCGA MC3 ----
# define which mutation types to load
# mutation_types <- c("3'Flank", "3'UTR", "5'Flank", "5'UTR", "Frame_Shift_Del", 
#                     "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Intron", 
#                     "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", 
#                     "RNA", "Silent", "Splice_Site", "Translation_Start_Site")

mutation_types <- c("3'UTR", "Missense_Mutation", "Silent")



# takes about 5-10 min to run
tcga_mc3 <- read.maf(
  maf = "~/Documents/4_projects/20220623_synonymous_m6A/20220623_gathering_m6A_sites/syn_mut_data/TCGA/mc3.v0.2.8.PUBLIC.maf.gz",
  # define which variants should be included
  vc_nonSyn = mutation_types
)

# make df
tcga_mc3_df <- as.data.frame(tcga_mc3@data)

# wrangle ----------------------------------------------------------------------
## add cancer type ----
patients <- tcga_mc3_df$Tumor_Sample_Barcode %>% str_sub(1, 12)
idx <- match(patients, sample_cancer_type$patient)
tcga_mc3_df$cancer_type <- sample_cancer_type$data_set[idx]

## remove MSI & POLE patients ----
tcga_mc3_df <- tcga_mc3_df %>% filter(!cancer_type %in% c("TCGA-MSI", "TCGA-POLE"))

## liftOver ----
# import chain file
chainObject <- import.chain("20230322_publication_figs/data/liftOver/hg19ToHg38.over.chain")

# prep data frame for liftOver
tcga_mc3_df <- tcga_mc3_df %>% mutate(Chromosome = paste0("chr", Chromosome))

# make GRanges object for LiftOver
tcga_mc3_df_GRanges <- makeGRangesFromDataFrame(tcga_mc3_df,
                                                keep.extra.columns = TRUE,
                                                ignore.strand = TRUE,
                                                seqinfo = NULL,
                                                seqnames.field = "Chromosome",
                                                start.field = "Start_Position",
                                                end.field = "End_Position",
                                                starts.in.df.are.0based = FALSE
)

# perform liftover
tcga_mc3_df_GRCh38 <- as.data.frame(liftOver(tcga_mc3_df_GRanges, chainObject))

# filter for SNPs 
tcga_mc3_df_GRCh38 <- tcga_mc3_df_GRCh38 %>% 
  filter(Variant_Type == "SNP")

# convert strand column from c(1,-1) to c("+","-")
tcga_mc3_df_GRCh38$STRAND[tcga_mc3_df_GRCh38$STRAND == 1] <- "+"
tcga_mc3_df_GRCh38$STRAND[tcga_mc3_df_GRCh38$STRAND == -1] <- "-"

## add sequence context ----
# wild-type sequence
tcga_mc3_df_GRCh38 <- tcga_mc3_df_GRCh38 %>%
  as.data.frame() %>%
  mutate(seqnames = gsub("^chr","", seqnames)) %>% 
  mutate(rna_seq_wt = getSeq(BSgenome.Hsapiens.NCBI.GRCh38,
                      seqnames, 
                      start = start - 4, 
                      end = end + 4,
                      strand = STRAND,
                      as.character = TRUE))


# Reverse complement sequences from minus strand -------------------------------

# Tumor_Seq_Allele2
# The tumor allele is sequences from DNA, so it needs to be reverse transcribed
tcga_mc3_df_GRCh38$Tumor_Seq_Allele2 <- DNAStringSet(tcga_mc3_df_GRCh38$Tumor_Seq_Allele2)

tcga_mc3_df_GRCh38$Tumor_Seq_Allele2_onRNA <- ifelse(
  tcga_mc3_df_GRCh38$STRAND == "-",
  as.character(reverseComplement(tcga_mc3_df_GRCh38$Tumor_Seq_Allele2)),
  as.character(tcga_mc3_df_GRCh38$Tumor_Seq_Allele2)
)

# generate mutated sequence
tcga_mc3_df_GRCh38$rna_seq_mut <- tcga_mc3_df_GRCh38$rna_seq_wt
substr(tcga_mc3_df_GRCh38$rna_seq_mut, 5, 5) <- tcga_mc3_df_GRCh38$Tumor_Seq_Allele2_onRNA %>% str_sub(-1, -1)

# revert back to character
tcga_mc3_df_GRCh38$Tumor_Seq_Allele2 <- as.character(tcga_mc3_df_GRCh38$Tumor_Seq_Allele2)

# find DRACH motifs ------------------------------------------------------------
# define DRACH motif
drach_motif <- "(A|G|T)(A|G)AC(A|C|T)"

# find DRACH motifs in wild-type sequence
# locate motif occurrences
DRACH_coord_wt <- str_locate(tcga_mc3_df_GRCh38$rna_seq_wt, drach_motif) %>% 
  as.data.frame()

DRACH_coord_mut <- str_locate(tcga_mc3_df_GRCh38$rna_seq_mut, drach_motif) %>% 
  as.data.frame()

# add motif location to df
tcga_mc3_df_GRCh38 <- tcga_mc3_df_GRCh38 %>% 
  mutate(motif_wt_coord = paste0(DRACH_coord_wt$start, "-", DRACH_coord_wt$end),
         motif_mut_coord = paste0(DRACH_coord_mut$start, "-", DRACH_coord_mut$end))

# add motif sequence to df
tcga_mc3_df_GRCh38 <- tcga_mc3_df_GRCh38 %>% 
   mutate(motif_wt = substr(tcga_mc3_df_GRCh38$rna_seq_wt, DRACH_coord_wt$start, DRACH_coord_wt$end),
          motif_mut = substr(tcga_mc3_df_GRCh38$rna_seq_mut, DRACH_coord_mut$start, DRACH_coord_mut$end))

# classify the motifs
tcga_mc3_df_GRCh38 <- tcga_mc3_df_GRCh38 %>% 
  mutate(motif_class = case_when(motif_wt_coord == "NA-NA" & motif_mut_coord == "NA-NA" ~ "no motif",
                                 motif_wt_coord == motif_mut_coord ~ "same position",
                                 motif_wt_coord == "NA-NA" & motif_mut_coord != "NA-NA" ~ "gain",
                                 motif_wt_coord != "NA-NA" & motif_mut_coord == "NA-NA" ~ "loss",
                                 motif_wt_coord != motif_mut_coord & motif_wt_coord != "NA-NA" & motif_mut_coord != "NA-NA" ~ "shifted position"))

# remove columns that only contain "."
tcga_mc3_df_GRCh38 <- tcga_mc3_df_GRCh38 %>%
  select_if(~ !all(. == "."))

# add motif ids
tcga_mc3_df_GRCh38 <- tcga_mc3_df_GRCh38 %>%
  mutate(gen_mut_motif_start = if_else(
    str_detect(str_sub(motif_mut_coord, 1, 1), "^[0-9]+$"),
    start - (5 - as.numeric(str_sub(motif_mut_coord, 1, 1))),
    as.numeric(NA))) %>% 
  mutate(mut_motif_id = paste0(seqnames, ":", gen_mut_motif_start))

# filter mutations that recur because they were observed in multiple samples of
# the same patient
tcga_mc3_df_GRCh38 <- tcga_mc3_df_GRCh38 %>% 
  mutate(bcr_patient_barcode = str_sub(Tumor_Sample_Barcode, 1, 12)) %>% 
  mutate(mutation_sample_id = paste0(bcr_patient_barcode, Hugo_Symbol, HGVSc, HGVSp_Short)) %>% 
  distinct(mutation_sample_id, .keep_all = TRUE)

# export
save(tcga_mc3_df_GRCh38, file = "20240403_gained_motifs/processed_data/001_tcga_mc3_motifs_anno.RData")

# cleanup ----------------------------------------------------------------------
# cleanup
rm(chainObject, DRACH_coord_mut, DRACH_coord_wt, sample_cancer_type, tcga_mc3_df, tcga_mc3_df_GRanges, drach_motif, idx, mutation_types, patients, tcga_mc3)
