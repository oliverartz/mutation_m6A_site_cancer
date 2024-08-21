# ______________________________________________________________________________
# Plot survival stats for patients with low m6A mutational load
# for each cancer type
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
count_m6A_mut <- read_csv("20230322_publication_figs/processed_data/m6A_load_quartiles.csv", show_col_types = FALSE)


# TCGA survival data ----
all_survival_data <- survivalTCGA(
  ACC.clinical, BLCA.clinical, BRCA.clinical,
  CESC.clinical, CHOL.clinical, COAD.clinical,
  COADREAD.clinical, DLBC.clinical, ESCA.clinical,
  FPPP.clinical, GBM.clinical, GBMLGG.clinical,
  HNSC.clinical, KICH.clinical, KIPAN.clinical,
  KIRC.clinical, KIRP.clinical, LAML.clinical,
  LGG.clinical, LIHC.clinical, LUAD.clinical,
  LUSC.clinical, MESO.clinical, OV.clinical,
  PAAD.clinical, PCPG.clinical, PRAD.clinical,
  READ.clinical, SARC.clinical, SKCM.clinical,
  STAD.clinical, STES.clinical, TGCT.clinical,
  THCA.clinical, THYM.clinical, UCEC.clinical,
  UCS.clinical, UVM.clinical
)

# analysis ---------------------------------------------------------------------
# add survival data to mutations
patients <- str_sub(count_m6A_mut$Tumor_Sample_Barcode, 1, 12)
idx <- match(patients, all_survival_data$bcr_patient_barcode)

count_m6A_mut$times <- all_survival_data$times[idx]
count_m6A_mut$patient.vital_status <- all_survival_data$patient.vital_status[idx]

# DEBUG
# cancer_data_set <- "TCGA-BRCA"
# cancer_data_set <- "TCGA-LUAD"

# define function to get survival stats for each cancer type
get_surv_pval <- function(cancer_data_set){
  
  # make data frame for cancer type
  df_temp <- count_m6A_mut %>% 
    filter(cancer_type == cancer_data_set)
  
  # filter for quantiles
  df_temp <- df_temp %>% filter(quantile_group %in% c("Below 1st Quartile", "Above 3rd Quartile"))
  
  # define survival object 
  surv_object <- Surv(time = df_temp$times, event = df_temp$patient.vital_status)
  
  # fit model
  fit <- surv_fit(surv_object ~ quantile_group, data = df_temp)
  
  # get p-value of difference in total mutations
  wilcox_p <- wilcox.test(df_temp$total_mut[df_temp$quantile_group == "Below 1st Quartile"], 
                          df_temp$total_mut[df_temp$quantile_group == "Above 3rd Quartile"])$p.value
  
  # print results
  data.frame(cancer_type = cancer_data_set,
             p_val = surv_pvalue(fit)$pval,
             wilcox_p_total_mut = wilcox_p)
}

# call function on all cancer types to get all p_values
p_val <- lapply(unique(count_m6A_mut$cancer_type), 
                FUN = get_surv_pval) %>% 
  do.call(rbind, .)

# plot -------------------------------------------------------------------------
p_val$cancer_type <- str_remove(p_val$cancer_type, "TCGA-")

p_017 <- p_val %>% 
  ggplot(aes(x = reorder(cancer_type, -p_val), y = p_val)) +
  geom_bar(stat = "identity", color = "black", fill = "grey", alpha = 0.8) +
  geom_bar(data = p_val[p_val$p_val < 0.05, ], stat = "identity", color = "black", fill = "orange", alpha = 0.8) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  theme_bw() +
  scale_x_discrete(limits = rev) +
  labs(title = "Kaplan-Meier p-values",
#       subtitle = "Stratification: Below 1st quartile vs. above 3rd quartile",
       x = "",
       y = "p-value") +
#  theme(axis.text = element_text(size = 6)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(color = "black"),
              axis.text = element_text(color = "black"),
              axis.title = element_text(color = "black"))

p_017

# export -----------------------------------------------------------------------
# define plot name as name of script
plot_name <- rstudioapi::getSourceEditorContext()$path %>%
  basename() %>% 
  str_replace(".R", ".png")

ggsave(filename = paste0("20230322_publication_figs/plots/", plot_name), width = 12, height = 6)

# export count and quartile data plus survival
write_csv(p_val, file = "20230322_publication_figs/processed_data/017_m6A_load_quartiles_survival_pval.csv")
write_csv(count_m6A_mut, file = "20230322_publication_figs/processed_data/017_m6A_load_quartiles_survival.csv")

# cleanup ----------------------------------------------------------------------
rm(all_survival_data, count_m6A_mut, p_val, idx, patients, plot_name, get_surv_pval)
