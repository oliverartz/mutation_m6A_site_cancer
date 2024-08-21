# ______________________________________________________________________________
# Plot survival curves for significant cancer types
# ______________________________________________________________________________

# clear workspace --------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, data.table, RTCGA.clinical, survival, survminer)

# load data --------------------------------------------------------------------
# m6A mutational load with quartiles
count_m6A_mut <- read_csv("20230322_publication_figs/processed_data/017_m6A_load_quartiles_survival.csv", show_col_types = FALSE)

# analysis ---------------------------------------------------------------------
# define data set
# cancer_data_set <- "TCGA-LUAD"

# define function to plot and export KM curves
make_KM_curve <- function(cancer_data_set){

# Kaplan Meier analysis
df_temp <- count_m6A_mut %>% 
  filter(cancer_type == cancer_data_set) %>% 
  filter(quantile_group %in% c("Below 1st Quartile", "Above 3rd Quartile"))

# define survival object 
surv_object <- Surv(time = df_temp$times, event = df_temp$patient.vital_status)

# fit model
fit <- surv_fit(surv_object ~ quantile_group, data = df_temp)

# plot -------------------------------------------------------------------------
p <- ggsurvplot(
  fit = fit,
  data = df_temp,
  pval = TRUE,
  conf.int = TRUE,
  palette = "jco",
  xlim = c(0, 3650),
  title = cancer_data_set,
  legend.title = "",
  xscale = "d_y",
  xlab = "Time (years)",
  break.time.by = 365.25,
  surv.median.line = "hv",
  risk.table = FALSE,
  #  ggtheme = theme_bw()
    legend.labs = c("Above 3rd Quartile", "Below 1st Quartile")
)

# export -----------------------------------------------------------------------
# define plot name as name of script
plot_name <- rstudioapi::getSourceEditorContext()$path %>%
  basename() %>% 
  str_replace(".R", "")

ggsave(filename = paste0("20230322_publication_figs/plots/", plot_name, "_", cancer_data_set, ".png"), width = 6, height = 6)

p
}

#cancer_types_plotting <- c("TCGA-THCA")
#sapply(cancer_types_plotting, make_KM_curve)

p_018 <- make_KM_curve("TCGA-THCA")

# cleanup ----------------------------------------------------------------------
rm(make_KM_curve, count_m6A_mut)

