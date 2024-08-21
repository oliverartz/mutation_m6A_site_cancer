## -----------------------------------------------------------------------------
## Purpose of script: Survival analysis for recurring DRACH gains
##
## Author: Oliver Artz
## Date Created: April 29, 2024
## Email: artzo@mskcc.org
##
## -----------------------------------------------------------------------------

# clear work space -------------------------------------------------------------
# rm(list = ls())
# gc()

# load libraries ---------------------------------------------------------------
options(repos = "https://cran.rstudio.com")
if (!require('pacman')) install.packages('pacman', dependencies = TRUE); library(pacman)

p_load(tidyverse, data.table, RTCGA.clinical, survival, survminer, gtsummary, broom, ggrepel)

# define parameters ------------------------------------------------------------
recurrence_threshold <- 3

# load data --------------------------------------------------------------------
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

# stage data
load("20230711_survival_stratify_stage_age_sex/processed_data/TCGA_stage_info.RData")

# annotated TCGA mutations
load("20240403_gained_motifs/processed_data/001_tcga_mc3_motifs_anno.RData")

# wrangle ----------------------------------------------------------------------
# add survival data
tcga_mc3_mod <- tcga_mc3_df_GRCh38 %>% 
  select(Tumor_Sample_Barcode, mut_motif_id, motif_class, Hugo_Symbol, HGVSc, HGVSp_Short) %>%
  mutate(bcr_patient_barcode = str_sub(Tumor_Sample_Barcode, 1, 12))

# tcga_mc3_mod <- left_join(tcga_mc3_mod, 
#                           all_survival_data, 
#                           by = c("bcr_patient_barcode" = "bcr_patient_barcode"),
#                           relationship = "many-to-many")

idx <- match(tcga_mc3_mod$bcr_patient_barcode, all_survival_data$bcr_patient_barcode)
tcga_mc3_mod$times <- all_survival_data$times[idx]
tcga_mc3_mod$patient.vital_status <- all_survival_data$patient.vital_status[idx]

# add stage data
tcga_mc3_mod <- left_join(tcga_mc3_mod, 
                          all_cancer_stages_mod, 
                          by = c("bcr_patient_barcode" = "bcr_patient_barcode"),
                          relationship = "many-to-many")

# remove mutations that were picked up in multiple samples per patient
tcga_mc3_mod <- tcga_mc3_mod %>% 
  mutate(mutation_sample_id = paste0(bcr_patient_barcode, Hugo_Symbol, HGVSc, HGVSp_Short)) %>% 
  distinct(mutation_sample_id, .keep_all = TRUE)

# analysis ---------------------------------------------------------------------
# count recurrence of DRACH gains per cancer type
gained_recurrence <- tcga_mc3_mod %>% 
  filter(motif_class == "gain") %>%
  group_by(mut_motif_id, cancer_type, stage_simple, Hugo_Symbol, HGVSc, HGVSp_Short) %>% 
  summarize(n = n()) %>% 
  filter(n >= recurrence_threshold)

# add list of bcr_patient_barcodes for each gain to the data frame in column "patients"
# make sure it is specific for cancer_type and stage_simple
gained_recurrence$patients <- sapply(1:nrow(gained_recurrence), function(i) {
  tcga_mc3_mod %>% 
    filter(motif_class == "gain",
           mut_motif_id == gained_recurrence$mut_motif_id[i], 
           cancer_type == gained_recurrence$cancer_type[i], 
           stage_simple == gained_recurrence$stage_simple[i], 
           Hugo_Symbol == gained_recurrence$Hugo_Symbol[i], 
           HGVSc == gained_recurrence$HGVSc[i], 
           HGVSp_Short == gained_recurrence$HGVSp_Short[i]) %>% 
    pull(bcr_patient_barcode)
})

# check whether all gained_recurrence$patients are distinct
# all(sapply(gained_recurrence$patients, function(x) length(x) == length(unique(x))))


# survival analysis ------------------------------------------------------------
# make data frame for analysis

# DEBUG ------------------------------------------------------------------------
# stage_of_interest <- gained_recurrence$stage_simple[1]
# mut_motif_of_interest <- gained_recurrence$mut_motif_id[1]
# cancer_type_of_interest <- gained_recurrence$cancer_type[1]
# END DEBUG --------------------------------------------------------------------

perform_cox_analysis <- function(stage_of_interest, mut_motif_of_interest, cancer_type_of_interest) {
  survival_data <- tcga_mc3_mod %>% 
  filter(stage_simple == stage_of_interest) %>% 
  filter(cancer_type == cancer_type_of_interest) %>% 
  mutate(has_gain = ifelse(mut_motif_id == mut_motif_of_interest, "Yes", "No"))
  
  patients_with_gain <- survival_data %>% 
    filter(has_gain == "Yes") %>% 
    distinct(Tumor_Sample_Barcode, .keep_all = TRUE)
  
  patients_without_gain <- survival_data %>% 
    filter(has_gain == "No") %>% 
    distinct(Tumor_Sample_Barcode, .keep_all = TRUE)
  
  survival_data <- rbind(patients_with_gain, patients_without_gain)

# Kaplan-Meier survival analysis
surv_object <- Surv(time = survival_data$times, event = survival_data$patient.vital_status)

# fit model
coxph_results <- coxph(surv_object ~ has_gain, data = survival_data) 
coxph_summary <- summary(coxph_results)
cox_p <- coxph_summary$coefficients[, "Pr(>|z|)"]
cox_hazard_ratio <- coxph_summary$coefficients[, "exp(coef)"]

# results
data.frame(
  p_value = cox_p,
  hazard_ratio = cox_hazard_ratio
)}

# apply the function to each row of gained_recurrence
start_time <- Sys.time()

results <- gained_recurrence %>%
  mutate(results = pmap(list(stage_simple, mut_motif_id, cancer_type), 
                        ~ perform_cox_analysis(..1, ..2, ..3))) %>% 
  unnest(cols = results)

end_time <- Sys.time()
end_time - start_time


# beepr::beep(2)

# plot -------------------------------------------------------------------------
# plot cox results
p006_cox <- results %>% 
  ggplot(aes(x = hazard_ratio, y = -log10(p_value))) +
  geom_point(shape = 21, fill = "grey", alpha = 0.3, aes(size = n)) +
  geom_point(data = filter(results, p_value < 0.05),
             shape = 21, fill = "orange", alpha = 0.3, aes(size = n)) +
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.2) +
  geom_text_repel(data = filter(results, p_value < 0.05),
                  aes(label = paste(Hugo_Symbol, HGVSp_Short)),
                  box.padding = 0.75,
                  size = 2.5,
                  segment.size = 0.2) +
  theme_bw() +
  theme(text = element_text(color = "black"),
        panel.grid.minor = element_blank()) +
  labs(x = "Hazard Ratio", y = "-log10(p-value)",
       title = paste0("Cox regression, min. recurrence: ", recurrence_threshold)) +
  facet_wrap(~ stage_simple)

p006_cox

# export -----------------------------------------------------------------------
# save plot
ggsave("20240403_gained_motifs/plots/006_cox_results.png", p006_cox, width = 12, height = 6, units = "in")

# save results
save(results, file = "20240403_gained_motifs/processed_data/006_cox_results.RData")

# save modified MC3 data
save(tcga_mc3_mod, file = "20240403_gained_motifs/processed_data/006_tcga_mc3_mod.RData")

# cleanup ----------------------------------------------------------------------
rm(gained_recurrence, all_cancer_stages_mod, all_survival_data, results, 
   tcga_mc3_df_GRCh38, tcga_mc3_mod, end_time, idx, recurrence_threshold,
   start_time, perform_cox_analysis)
