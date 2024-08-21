## -----------------------------------------------------------------------------
## Purpose of script: Plot KM-curves of interesting mutations
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

p_load(tidyverse, data.table, survival, survminer, cowplot)

# set parameters ---------------------------------------------------------------

# load data --------------------------------------------------------------------
load("20240403_gained_motifs/processed_data/006_cox_results.RData")
load("20240403_gained_motifs/processed_data/006_tcga_mc3_mod.RData")

load("20240403_gained_motifs/processed_data/001_tcga_mc3_motifs_anno.RData")


# wrangle ----------------------------------------------------------------------
# take only significant comparisons
results_mod <- results %>% 
  filter(p_value < 0.05)

# analysis ---------------------------------------------------------------------
# DEBUG 
# row_num <- 1
# END DEBUG

# Kaplan Meier analysis --------------------------------------------------------
draw_KM <- function(row_num){
cancer_data_set <- results_mod$cancer_type[row_num]
stage_of_interest <- results_mod$stage_simple[row_num]
mut_motif_id_of_interest <- results_mod$mut_motif_id[row_num]
hugo_symbol_of_interest <- results_mod$Hugo_Symbol[row_num]
hgvsc_of_interest <- results_mod$HGVSc[row_num]

df_temp <- tcga_mc3_mod %>% 
  filter(cancer_type == cancer_data_set) %>%
  filter(stage_simple == stage_of_interest)

# annotate mutation
df_temp <- df_temp %>% 
  mutate(has_mutation = ifelse(mut_motif_id == mut_motif_id_of_interest, "Yes", "No"))

# filter for patients with and without gain
patients_with_gain <- df_temp %>% 
  filter(has_mutation == "Yes") %>% 
  distinct(Tumor_Sample_Barcode, .keep_all = TRUE)

patients_without_gain <- df_temp %>% 
  filter(has_mutation == "No") %>% 
  distinct(Tumor_Sample_Barcode, .keep_all = TRUE)

# combine datasets
df_temp <- rbind(patients_with_gain, patients_without_gain)

# count patients with or without mutation
n_mut <- df_temp %>% group_by(has_mutation) %>% summarise(n = n())

# define survival object 
surv_object <- Surv(time = df_temp$times, event = df_temp$patient.vital_status)

# fit model
fit <- surv_fit(surv_object ~ has_mutation, data = df_temp)

# plot 
p <- ggsurvplot(
  fit = fit,
  data = df_temp,
  pval = TRUE,
  conf.int = TRUE,
  palette = "jco",
  xlim = c(0, 2200),
  title = paste0(cancer_data_set %>% str_remove("TCGA-"), ", ", 
                 stage_of_interest, ", ",
                 hugo_symbol_of_interest, " ",
                 hgvsc_of_interest),
  legend.title = "",
  xscale = "d_y",
  xlab = "Time (years)",
  break.time.by = 365.25,
  surv.median.line = "hv",
  risk.table = FALSE,
  legend = "top",
  ggtheme = theme_bw() + 
    theme(plot.title = element_text(size = 9, color = "black"),
          axis.title = element_text(size = 9, color = "black"),
          axis.text = element_text(size = 9, color = "black"),
          text = element_text(color = "black"),
          legend.text = element_text(size = 5),
          legend.key.size = unit(0.5, "cm")),
  pval.size = 3,
  pval.coord = c(1500, 0.8),
  legend.labs = c(paste0("wildtype (n = ", n_mut$n[1], ")"),
                  paste0("mutated (n = ", n_mut$n[2], ")"))
)

# extract only the plots
p$plot
}

# make plots
p_KM <- lapply(1:nrow(results_mod), draw_KM)

# combine plots ----------------------------------------------------------------
# no panel labels
p_007_KM_curves <- plot_grid(plotlist = p_KM, nrow = 1)
p_007_KM_curves

# export -----------------------------------------------------------------------
ggsave("20240403_gained_motifs/plots/007_KM_curves.png", p_007_KM_curves, width = 10, height = 9, units = "in")

# cleanup ----------------------------------------------------------------------
rm(p_KM, results, results_mod, tcga_mc3_mod, draw_KM)
