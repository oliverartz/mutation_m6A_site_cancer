# ______________________________________________________________________________
# Make m6A data set 
# Homogenize m6A data sets from the literature and make one singular data base
# ______________________________________________________________________________

# clear workspace ----
rm(list = ls())
gc()

# load libraries ----
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, readxl, BSgenome, BSgenome.Hsapiens.NCBI.GRCh38, writexl, maftools)

# m6A data sets ----
## Hu 2022 ----
# load data
# Hu 2022 reports data for different cell lines in different sheets of the excel file
# do not load ribo- data, to make sure we only get messenger RNAs and exclude miRNAs, SNORs, etc.
sheets_to_import <- c(
  "HeLa polyA", 
  "HEK293 polyA", 
  "HEK293 ribo-", 
  "HepG2 polyA", 
  "HSPC ribo- d0", 
  "HSPC ribo- d3", 
  "HSPC ribo- d6", 
  "HSPC ribo- d9",
  NULL
)

Hu_2022 <- lapply(sheets_to_import, 
                  function(x) read_excel("~/Documents/4_projects/20220623_synonymous_m6A/20220623_gathering_m6A_sites/m6A_data/01_41587_2022_1243_MOESM3_ESM.xlsx", 
                                         sheet = x))

for (i in seq(1:length(Hu_2022))) {
  Hu_2022[[i]]$experiment <- sheets_to_import[i]
}

Hu_2022 <- Hu_2022 %>% do.call(rbind, .)

# change colnames
colnames(Hu_2022) <- c("chr", "Start", "Strand", "m6A_frac", "5-mer", "experiment")

Hu_2022 <- Hu_2022 %>% 
  # add end position for coordinate
  mutate(End = Start) %>% 
  # select important columns
  dplyr::select(chr, Start, End, Strand) %>% 
  # re-format data frame for getSeq
  mutate(chr = gsub("^chr","", chr),
         chr = gsub("23", "X", chr),
         chr = gsub("24", "Y", chr),
         chr = gsub("25", "MT", chr)) %>%
  # remove mitochondrial hits
  filter(chr != "M") %>% 
  # add sequence context
  mutate(seq = getSeq(BSgenome.Hsapiens.NCBI.GRCh38,
                      chr, 
                      start = Start - 2, 
                      end = End + 2, 
                      strand = Strand, 
                      as.character = TRUE)) %>% 
  # add columns for later identification of data set
  mutate(data_set = "Hu 2022",
         category = "m6A")



# cleanup
rm(i, sheets_to_import)

## Tegowski bulk ----
# load data
Tegowski_2022_bulk <- read_excel("20230322_publication_figs/data/m6A_data/02_1-s2.0-S1097276521011436-mmc2.xlsx", 
                                 skip = 1,
                                 col_types = c("text","numeric", "numeric", "text", "numeric", "text"))

# change colnames
colnames(Tegowski_2022_bulk) <- c("chr","Start","End","Gene","Edit_Ratio","Strand")

Tegowski_2022_bulk <- Tegowski_2022_bulk %>%
  # select important columns
  select(chr, Start, End, Strand) %>% 
  # make end coordinate = start coordinate
  mutate(Start = case_when(Strand == "-" ~ Start + 2,
                           Strand == "+" ~ Start),
         End = Start) %>% 
  # remove NAs
  drop_na() %>% 
  # add sequence context
  mutate(seq = getSeq(BSgenome.Hsapiens.NCBI.GRCh38,
                      chr, 
                      start = Start - 2, 
                      end = End + 2, 
                      strand = Strand, 
                      as.character = TRUE)) %>% 
  # add columns for later identification of data set
  mutate(data_set = "Tegowski 2022 bulk",
         category = "m6A")

## Tegowski sc ----
# load data
Tegowski_2022_sc <- read_excel("20230322_publication_figs/data/m6A_data/03_1-s2.0-S1097276521011436-mmc3.xlsx", 
                               skip = 1, 
                               sheet = "10xGenomicsHEK293T_m6Asites",
                               col_types = c("text","numeric", "numeric", "text", "text", "text"))

# change colnames
colnames(Tegowski_2022_sc) <- c("chr", "Start", "End", "Gene", "Motif", "Strand")

Tegowski_2022_sc <- Tegowski_2022_sc %>%
  # select important columns
  select(chr, Start, End, Strand) %>% 
  # re-format data frame for getSeq
  mutate(chr = gsub("^chr","", chr)) %>%
  # make end coordinate = start coordinate
  mutate(Start = case_when(Strand == "-" ~ Start + 2,
                           Strand == "+" ~ Start),
         End = Start) %>% 
  # remove NAs
  drop_na() %>% 
  # add sequence context
  mutate(seq = getSeq(BSgenome.Hsapiens.NCBI.GRCh38,
                      chr, 
                      start = Start - 2, 
                      end = End + 2, 
                      strand = Strand, 
                      as.character = TRUE)) %>% 
  # add columns for later identification of data set
  mutate(data_set = "Tegowski 2022 sc",
         category = "m6A")

## Tegowski SMARTseq2 ----
# load data
Tegowski_2022_SMARTseq2 <- read_excel("20230322_publication_figs/data/m6A_data/05_1-s2.0-S1097276521011436-mmc4.xlsx", 
                                      sheet = "SMARTseq2_HighConfidenceSites", 
                                      skip = 1,
                                      col_types = c("text","numeric", "numeric", "text", "text", "text"))

# change colnames
colnames(Tegowski_2022_SMARTseq2) <- c("chr", "Start", "End", "Gene", "Motif", "Strand")

Tegowski_2022_SMARTseq2 <- Tegowski_2022_SMARTseq2 %>%
  # select important columns
  select(chr, Start, End, Strand) %>% 
  # make end coordinate = start coordinate
  mutate(Start = case_when(Strand == "-" ~ Start + 2,
                           Strand == "+" ~ Start),
         End = Start) %>% 
  # remove NAs
  drop_na() %>% 
  # add sequence context
  mutate(seq = getSeq(BSgenome.Hsapiens.NCBI.GRCh38,
                      chr, 
                      start = Start - 2, 
                      end = End + 2, 
                      strand = Strand, 
                      as.character = TRUE)) %>% 
  # add columns for later identification of data set
  mutate(data_set = "Tegowski 2022 SMARTseq2",
         category = "m6A")

## Garcia 2019 ----
Garcia_2019 <- read_excel("20230322_publication_figs/data/m6A_data/04_1-s2.0-S0092867419306762-mmc6.xlsx") %>% 
  filter(prevKnown_ConfidenceGroups %in% c("high", "highest"))


Garcia_2019 <- Garcia_2019 %>% 
  # select important columns
  select(chr, start, end, strand)

# liftOver from hg19 to GRCh38
## import chain file
chainObject <- import.chain("20230322_publication_figs/data/liftOver/hg19ToHg38.over.chain")

## make gRanges object for liftover
Garcia_2019_gRanges <- makeGRangesFromDataFrame(Garcia_2019,
                                                keep.extra.columns = FALSE,
                                                ignore.strand = FALSE,
                                                seqinfo = NULL,
                                                start.field = "start",
                                                end.field = "end",
                                                strand.field = "strand",
                                                starts.in.df.are.0based = FALSE)

# perform liftover
Garcia_2019_GRCh38 <- as.data.frame(liftOver(Garcia_2019_gRanges, 
                                             chainObject)) %>%
  # add chr column
  mutate(chr = gsub("^chr","", seqnames)) %>% 
  # select important columns
  select(chr, start, end, strand) %>% 
  # add sequence context
  mutate(seq = getSeq(BSgenome.Hsapiens.NCBI.GRCh38,
                      chr, 
                      start = start - 2, 
                      end = end + 2, 
                      strand = strand, 
                      as.character = TRUE))
# change colnames
colnames(Garcia_2019_GRCh38) <- c("chr","Start","End","Strand","seq")


# add columns for later identification of data set
Garcia_2019_GRCh38 <- Garcia_2019_GRCh38 %>% 
  mutate(data_set = "Garcia 2019",
         category = "m6A")

# cleanup
Garcia_2019 <- Garcia_2019_GRCh38
rm(Garcia_2019_gRanges, Garcia_2019_GRCh38)

## Meyer 2019 ----
Meyer_2019 <- read_excel("20230322_publication_figs/data/m6A_data/06_41592_2019_570_MOESM4_ESM.xlsx",
                         sheet = "C2U sites high stringency")

# Meyer 2019 report the C2U site, which is immediately downstream of the methylated A as the end position
# Analogous to the other studies I will report the site of the A
Meyer_2019 <- Meyer_2019 %>%
  mutate(start = `C2U start`,
         end = `C2U start`) %>%
  # select important columns
  select(Chr, start, end, Strand)

# liftOver from hg19 to GRCh38
## make gRanges object for liftover
Meyer_2019_gRanges <- makeGRangesFromDataFrame(Meyer_2019,
                                               keep.extra.columns = FALSE,
                                               ignore.strand = FALSE,
                                               seqinfo = NULL,
                                               start.field = "start",
                                               end.field = "end",
                                               strand.field = "Strand",
                                               starts.in.df.are.0based = FALSE)

# perform liftover
Meyer_2019_GRCh38 <- as.data.frame(liftOver(Meyer_2019_gRanges, 
                                            chainObject)) %>%
  # add chr column
  mutate(chr = gsub("^chr","", seqnames)) %>% 
  mutate(start = case_when(strand == "-" ~ start + 2,
                           strand == "+" ~ start),
         end = start) %>% 
  # select important columns
  select(chr, start, end, strand) %>% 
  # add sequence context
  mutate(seq = getSeq(BSgenome.Hsapiens.NCBI.GRCh38,
                      chr, 
                      start = start - 2, 
                      end = end + 2, 
                      strand = strand, 
                      as.character = TRUE))
# rename columns
colnames(Meyer_2019_GRCh38) <- c("chr","Start", "End", "Strand", "seq")

# add columns for later identification of data set
Meyer_2019 <- Meyer_2019_GRCh38 %>% 
  mutate(data_set = "Meyer 2019",
         category = "m6A")

# cleanup
rm(Meyer_2019_GRCh38, Meyer_2019_gRanges)

## Linder 2015 ----
Linder_2015a <- read_excel("20230322_publication_figs/data/m6A_data/07_41592_2015_BFnmeth3453_MOESM77_ESM.xlsx") %>% 
  select(-"PeakID (cmd=annotatePeaks.pl hek293.abcam.CIMS.m6A.9536.bed hg19)")

Linder_2015b <- read_excel("20230322_publication_figs/data/m6A_data/08_41592_2015_BFnmeth3453_MOESM78_ESM.xlsx") %>% 
  select(-"PeakID (cmd=annotatePeaks.pl hek293.sysy.CITS.m6A.mRNA.6543.bed hg19)")

# merge both data sets
Linder_2015 <- rbind(Linder_2015a, Linder_2015b) %>% 
  select(Chr, Start, End, Strand)

# liftOver from hg19 to GRCh38
## make gRanges object for liftover
Linder_2015_gRanges <- makeGRangesFromDataFrame(Linder_2015,
                                                keep.extra.columns = FALSE,
                                                ignore.strand = FALSE,
                                                seqinfo = NULL,
                                                start.field = "Start",
                                                end.field = "End",
                                                strand.field = "Strand",
                                                starts.in.df.are.0based = FALSE)

# perform liftover
Linder_2015_GRCh38 <- as.data.frame(liftOver(Linder_2015_gRanges, 
                                             chainObject)) %>%
  # add chr column
  mutate(chr = gsub("^chr","", seqnames)) %>% 
  # select important columns
  select(chr, start, end, strand) %>%
  # remove mitochondrian reads
  filter(chr != "M") %>% 
  # add sequence context
  mutate(seq = getSeq(BSgenome.Hsapiens.NCBI.GRCh38,
                      chr, 
                      start = start - 2, 
                      end = end + 2, 
                      strand = strand, 
                      as.character = TRUE))

# rename columns
colnames(Linder_2015_GRCh38) <- c("chr","Start", "End", "Strand", "seq")

# add columns for later identification of data set
Linder_2015 <- Linder_2015_GRCh38 %>% 
  mutate(data_set = "Linder 2019",
         category = "m6A")

# cleanup
rm(Linder_2015a, Linder_2015b, Linder_2015_GRCh38, Linder_2015_gRanges)

## Liu 2023 ----
Liu_2023 <- read_excel("20230322_publication_figs/data/m6A_data/41587_2022_1487_MOESM3_ESM.xlsx") %>% 
  select("Chr", "Sites", "Strand")

# change colnames
colnames(Liu_2023) <- c("chr", "Start", "Strand")

Liu_2023 <- Liu_2023 %>% 
  # add end position for coordinate
  mutate(End = Start) %>% 
  # select important columns
  select(chr, Start, End, Strand) %>% 
  # re-format data frame for getSeq
  mutate(chr = gsub("^chr","", chr)) %>%
  # remove mitochondrial hits
  filter(chr != "M") %>% 
  # add sequence context
  mutate(seq = getSeq(BSgenome.Hsapiens.NCBI.GRCh38,
                      chr, 
                      start = Start - 2, 
                      end = End + 2, 
                      strand = Strand, 
                      as.character = TRUE)) %>% 
  # add columns for later identification of data set
  mutate(data_set = "Liu 2023",
         category = "m6A")

## merge data sets ----
total_m6A_sites <- rbind(Garcia_2019, 
                         Hu_2022,
                         Tegowski_2022_SMARTseq2,
                         Tegowski_2022_bulk,
                         Tegowski_2022_sc,
                         Meyer_2019,
                         Linder_2015,
                         Liu_2023)

# cleanup
rm(Garcia_2019, 
   Hu_2022,
   Tegowski_2022_SMARTseq2,
   Tegowski_2022_bulk,
   Tegowski_2022_sc,
   Meyer_2019,
   Linder_2015,
   Liu_2023)

## add additional info columns ----
### Annotate DRACH motifs ----
# define DRACH
DRACH <- "(A|G|T)(A|G)AC(A|C|T)"

# add DRACH information to data frame
total_m6A_sites <- total_m6A_sites %>% mutate(DRACH = 0)
total_m6A_sites$DRACH[grepl(DRACH, total_m6A_sites$seq)] <- 1

# add motif coordinate
total_m6A_sites <- total_m6A_sites %>% 
  mutate(start_motif = Start - 2, 
         end_motif = End + 2,
         `Mutation Genome Position` = paste0(chr, ":", Start, Strand))

## make data set for unique sites ----
total_m6A_unique <- total_m6A_sites %>% distinct(`Mutation Genome Position`, .keep_all = TRUE) %>% select(-data_set)

### TESTING m6A DATA SET #######################################################
# Are the DRACH sites for each data set balanced for each strand?
test <- total_m6A_sites %>% 
  group_by(data_set, DRACH, Strand) %>% 
  count() %>% 
  group_by(data_set, Strand) %>% 
  mutate(prop_count = n/sum(n))

rm(test)

## export data sets ----
write_xlsx(total_m6A_sites,"20230322_publication_figs/processed_data/001_total_m6A_sites.xlsx")
write_xlsx(total_m6A_unique,"20230322_publication_figs/processed_data/001_total_m6A_unique.xlsx")

# TCGA cancer types ----
# load cancer types
sample_cancer_type <- read_delim("~/Documents/4_projects/20220623_synonymous_m6A/20220623_gathering_m6A_sites/syn_mut_data/TCGA_James/20220801_data/id2type.txt", delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE)

# rename columns
colnames(sample_cancer_type) <- c("patient", "data_set")

# TCGA MC3 ----
# define which mutation types to load
mutation_types <- c("3'Flank", "3'UTR", "5'Flank", "5'UTR", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Intron", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "RNA", "Silent", "Splice_Site", "Translation_Start_Site")

# takes about 5-10 min to run
tcga_mc3 <- read.maf(maf = "~/Documents/4_projects/20220623_synonymous_m6A/20220623_gathering_m6A_sites/syn_mut_data/TCGA/mc3.v0.2.8.PUBLIC.maf.gz",
                     # define which variants should be included
                     vc_nonSyn = mutation_types)
rm(mutation_types)

# make df
tcga_mc3_df <- as.data.frame(tcga_mc3@data)

#rm(tcga_mc3)

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

# cleanup
rm(tcga_mc3_df, tcga_mc3_df_GRanges, chainObject, patients)

# select important columns
tcga_mc3_df_GRCh38 <- tcga_mc3_df_GRCh38 %>%
  dplyr::select(
    seqnames,
    start,
    end,
    STRAND,
    Hugo_Symbol,
    HGNC_ID,
    cancer_type,
    Tumor_Sample_Barcode,
    Variant_Classification,
    Variant_Type,
    HGVSc,
    HGVSp_Short
  )

# filter for SNPs 
tcga_mc3_df_GRCh38 <- tcga_mc3_df_GRCh38 %>% 
  filter(Variant_Type == "SNP")

# convert strand column from c(1,-1) to c("+","-")
tcga_mc3_df_GRCh38$STRAND[tcga_mc3_df_GRCh38$STRAND == 1] <- "+"
tcga_mc3_df_GRCh38$STRAND[tcga_mc3_df_GRCh38$STRAND == -1] <- "-"

## annotate m6A motifs GRanges -------------------------------------------------
# make GRanges objects
# all TCGA mutations
tcga_mc3_GRanges <- makeGRangesFromDataFrame(tcga_mc3_df_GRCh38,
                                             keep.extra.columns = FALSE,
                                             ignore.strand = FALSE,
                                             seqinfo = NULL,
                                             start.field = "start",
                                             end.field = "end",
                                             strand.field = "STRAND",
                                             starts.in.df.are.0based = FALSE)

# motif around m6A site
total_m6A_unique_GRanges <- makeGRangesFromDataFrame(total_m6A_unique %>% mutate(chr = paste0("chr", chr)),
                                                     keep.extra.columns = FALSE,
                                                     ignore.strand = FALSE,
                                                     seqinfo = NULL,
                                                     start.field = "start_motif",
                                                     end.field = "end_motif",
                                                     strand.field = "Strand",
                                                     starts.in.df.are.0based = FALSE)

# Are there overlapping m6A motifs?
n_overlaps <- countOverlaps(total_m6A_unique_GRanges, total_m6A_unique_GRanges)
find_overlaps <- findOverlaps(total_m6A_unique_GRanges, total_m6A_unique_GRanges)

total_m6A_unique[c(4, 286298), ]

n_overlaps_mut <- countOverlaps(tcga_mc3_GRanges, total_m6A_unique_GRanges)

n_overlaps_mut %>% table()

# find overlaps
idx_methylated <- overlapsAny(tcga_mc3_GRanges, total_m6A_unique_GRanges, ignore.strand = FALSE)

# annotate mutations
tcga_mc3_df_GRCh38$in_methylated_motif <- ifelse(idx_methylated == TRUE, 1, 0)

## add sequence context ----
tcga_mc3_df_GRCh38 <- tcga_mc3_df_GRCh38 %>%
  as.data.frame() %>%
  mutate(seqnames = gsub("^chr","", seqnames)) %>% 
  mutate(seq = getSeq(BSgenome.Hsapiens.NCBI.GRCh38,
                      seqnames, 
                      start = start - 4, 
                      end = end + 4,
                      strand = STRAND,
                      as.character = TRUE))

## annotate DRACH motif in sequence
## Is the mutation within the context of a DRACH (regardless of methylation status)
tcga_mc3_df_GRCh38 <- tcga_mc3_df_GRCh38 %>% 
  mutate(in_DRACH = if_else(grepl(DRACH, seq), 1, 0))

table(tcga_mc3_df_GRCh38$in_methylated_motif, tcga_mc3_df_GRCh38$in_DRACH)

tcga_mc3_df_GRCh38 %>% 
  group_by(in_methylated_motif, in_DRACH) %>% 
  count()

## annotate DRACH motif position -----------------------------------------------
# get motif position within extracted sequence
tcga_mc3_df_GRCh38$mc3_DRACH_start <- str_locate(tcga_mc3_df_GRCh38$seq, 
                                                 pattern = DRACH) %>% 
  as.data.frame() %>% 
  pull(start)

# calculate motif start position
tcga_mc3_df_GRCh38 <- tcga_mc3_df_GRCh38 %>% 
  mutate(
    motif_start = case_when(
      tcga_mc3_df_GRCh38$STRAND == "+" ~ tcga_mc3_df_GRCh38$start - tcga_mc3_df_GRCh38$mc3_DRACH_start + 1,
      tcga_mc3_df_GRCh38$STRAND == "-" ~ 5 - (tcga_mc3_df_GRCh38$mc3_DRACH_start - tcga_mc3_df_GRCh38$start))) %>% 
  select(-mc3_DRACH_start)

# add motif ID to find recurrence of specific motifs
tcga_mc3_df_GRCh38$motif_id <- paste0(tcga_mc3_df_GRCh38$seqnames,":", as.character(tcga_mc3_df_GRCh38$motif_start), "_", tcga_mc3_df_GRCh38$STRAND)

# substitute WT bases with mutated bases in motif sequences
#alt_allele <- str_sub(tcga_mc3_df_GRCh38$HGVSc, -1)
#tcga_mc3_df_GRCh38 <- tcga_mc3_df_GRCh38 %>% mutate(seq_mutated = seq)
#substr(tcga_mc3_df_GRCh38$seq_mutated, 5, 5) <- alt_allele

# cleanup
#rm(alt_allele)

## add recurrence of mutation to data frame ----
# make id for each mutation
tcga_mc3_df_GRCh38 <- tcga_mc3_df_GRCh38 %>% 
  mutate(mutation_id = paste0(seqnames, start, end, STRAND, Hugo_Symbol, HGVSc))

# count recurrence of specific mutations
mutation_recur <- tcga_mc3_df_GRCh38$mutation_id %>% table() %>% as.data.frame()

# add recurrence information to df
idx <- match(tcga_mc3_df_GRCh38$mutation_id, mutation_recur$.)
tcga_mc3_df_GRCh38$recurrence_total <- mutation_recur$Freq[idx]

## add recurrence of mutation per cancer type ----
# count recurrence
mutation_recur <- tcga_mc3_df_GRCh38 %>% 
  group_by(cancer_type, mutation_id) %>% 
  summarize(recurrence_cancer_type = n())

tcga_mc3_df_GRCh38 <- left_join(tcga_mc3_df_GRCh38, mutation_recur, by = c("cancer_type", "mutation_id"))

# cleanup
rm(mutation_recur)

## annotate cancer genes ----
oncokb_gene_list <- read_delim("20221025_coexpression_machinery/data/cancerGeneList.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)

## add entrezgene_id information to OncoKB genes -------------------------------
# extract oncogene and tumor suppressor information
oncokb_df <- oncokb_gene_list %>%
  dplyr::select(`Hugo Symbol`, `Is Oncogene`, `Is Tumor Suppressor Gene`) %>%
  mutate(cancer_gene = ifelse(`Is Oncogene` == "Yes" & `Is Tumor Suppressor Gene` == "Yes", "OG & TS",
                              ifelse(`Is Oncogene` == "No" & `Is Tumor Suppressor Gene` == "No", "OnkoKB, neither",
                                     ifelse(`Is Oncogene` == "Yes", "OG",
                                            ifelse(`Is Tumor Suppressor Gene` == "Yes", "TS", "not defined")
                                     ))))

idx <- match(tcga_mc3_df_GRCh38$Hugo_Symbol, oncokb_df$`Hugo Symbol`)
tcga_mc3_df_GRCh38$cancer_gene <- oncokb_df$cancer_gene[idx]
tcga_mc3_df_GRCh38$cancer_gene[is.na(tcga_mc3_df_GRCh38$cancer_gene)] <- "not in OncoKB"

rm(oncokb_gene_list, idx, oncokb_df)

## add DRACH motif sequence ----------------------------------------------------
# define DRACH motif
# pattern <- "(A|G|T)(A|G)AC(A|C|T)"
# 
# # locate motif occurrences
# DRACH_coord_wt <- str_locate(tcga_mc3_df_GRCh38$seq, pattern) %>% as.data.frame()
# DRACH_coord_mut <- str_locate(tcga_mc3_df_GRCh38$seq_mutated, pattern) %>% as.data.frame()
# 
# # extract wt and mutated motif sequences and add to data frame
# tcga_mc3_df_GRCh38 <- tcga_mc3_df_GRCh38 %>% 
#   mutate(motif_wt = substr(tcga_mc3_df_GRCh38$seq, DRACH_coord_wt$start, DRACH_coord_wt$end),
#          motif_mut = substr(tcga_mc3_df_GRCh38$seq_mutated, DRACH_coord_mut$start, DRACH_coord_mut$end))
# 
# # add information about gained/lost motifs
# tcga_mc3_df_GRCh38$motif_consequence <- NA
# tcga_mc3_df_GRCh38$motif_consequence[is.na(DRACH_coord_wt$start) & is.na(DRACH_coord_mut$start)] <- "no motif"
# tcga_mc3_df_GRCh38$motif_consequence[is.na(DRACH_coord_wt$start) == FALSE & is.na(DRACH_coord_mut$start) == FALSE] <- "neutral"
# tcga_mc3_df_GRCh38$motif_consequence[is.na(DRACH_coord_wt$start) == TRUE & is.na(DRACH_coord_mut$start) == FALSE] <- "gain"
# tcga_mc3_df_GRCh38$motif_consequence[is.na(DRACH_coord_wt$start) == FALSE & is.na(DRACH_coord_mut$start) == TRUE] <- "loss"
# 
# # cleanup
# rm(pattern, DRACH_coord_mut, DRACH_coord_wt)

# ## add coordinates of motif ----------------------------------------------------
# # locate motif within sequence
# loc_wt <- str_locate(tcga_mc3_df_GRCh38$seq, tcga_mc3_df_GRCh38$motif_wt) %>% as.data.frame()
# loc_mut <- str_locate(tcga_mc3_df_GRCh38$seq_mut, tcga_mc3_df_GRCh38$motif_mut) %>% as.data.frame()
# 
# # consolidate start coordinate within sequence for wt and mut
# loc_start_both <- coalesce(loc_wt$start, loc_mut$start)
# loc_end_both <- coalesce(loc_wt$end, loc_mut$end)
# 
# # determine genomic locus from SNV coordinate and add to data frame
# # tcga_mc3_df_GRCh38$motif_start <- tcga_mc3_df_GRCh38$start - 5 + loc_start_both
# 
# tcga_mc3_df_GRCh38$loc_start_both <- coalesce(loc_wt$start, loc_mut$start)
# tcga_mc3_df_GRCh38$loc_end_both <- coalesce(loc_wt$end, loc_mut$end)
# 
# tcga_mc3_df_GRCh38 <- tcga_mc3_df_GRCh38 %>% mutate(motif_start = case_when(STRAND == "+" ~ start - 5 + loc_start_both,
#                                                                             STRAND == "-" ~ start + 5 - loc_end_both))
# 
# tcga_mc3_df_GRCh38 <- tcga_mc3_df_GRCh38 %>% mutate(motif_start = ifelse(motif_consequence == "no motif", "no motif", paste0("chr", seqnames, ":", motif_start)))

# cleanup
#rm(loc_wt, loc_mut, loc_start_both)

## annotate m6A motifs ---------------------------------------------------------
# DRACH_meth describes whether a given DRACH motif is methylated
# checks whether any mutation lies in DRACH motif
# checks whether given DRACH motif is methylated
# we focus only on mutations in DRACH motifs because we will be able to predict the effect on methylation
# If there is an m6A site, which is not within DRACH, the column value will be "no"

# total_m6A_unique$motif_id <- paste0("chr", total_m6A_unique$chr, ":", total_m6A_unique$start_motif,":", total_m6A_unique$Strand)
# tcga_mc3_df_GRCh38$motif_id <- paste0(tcga_mc3_df_GRCh38$motif_start, ":", tcga_mc3_df_GRCh38$STRAND)
# 
# idx <- match(tcga_mc3_df_GRCh38$motif_id, total_m6A_unique$motif_id)
# 
# tcga_mc3_df_GRCh38$DRACH_meth <- "yes"
# tcga_mc3_df_GRCh38$DRACH_meth[is.na(idx)] <- "no"
# 
# # cleanup
# rm(idx)

# filter mutations that recur because they were observed in multiple samples of
# the same patient
tcga_mc3_df_GRCh38 <- tcga_mc3_df_GRCh38 %>% 
  mutate(bcr_patient_barcode = str_sub(Tumor_Sample_Barcode, 1, 12)) %>% 
  mutate(mutation_sample_id = paste0(bcr_patient_barcode, Hugo_Symbol, HGVSc, HGVSp_Short)) %>% 
  distinct(mutation_sample_id, .keep_all = TRUE)

# export annotated mutation data
# too big to export to csv
save(tcga_mc3_df_GRCh38, file = "20230322_publication_figs/processed_data/001_tcga_mc3_df_GRCh38.RData")
