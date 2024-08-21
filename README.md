# The role of recurrent somatic mutations that alter conserved m<sup>6</sup>A motifs in human cancer

## Authors
Oliver Artz<sup>1</sup>, James R. White<sup>2</sup>, Benoit Rousseau<sup>1</sup>, Guillem Argiles<sup>1</sup>, Michael B. Foote<sup>1</sup>, Paul Johannet<sup>1</sup>, Miteshkumar Patel<sup>1</sup>, Somer Abdelfattah<sup>1</sup>, Shrey Patel<sup>1</sup>, Callahan Wilde<sup>1</sup>, David Mieles<sup>1</sup>, and Luis A. Diaz, Jr<sup>1</sup><sup>*</sup>

<sup>1</sup> Division of Solid Tumor Oncology, Department of Medicine, Memorial Sloan Kettering, New York City, NY, USA

<sup>2</sup> Resphera Biosciences, Baltimore, MD, USA

<sup>*</sup> To whom correspondence should be addressed. Email: diazl5@mskcc.org



## Abstract
Will be added after acceptance.

## Usage
This repository contains essential scripts to reproduce the analyses and visualizations of the data presented in the manuscript.

### Data availability
Original data was downloaded from `TCGA` and publications mentioned in `Table_S1` of the manuscript.

### Analysis
Data analysis was performed as described in the manuscript. For additional information, please refer to the scripts below.

The m<sup>6</sup>A data set was compiled using the following scripts:
- `001_make_m6A_data_set.R`
- `001.1_make_motifs_data_set.R`

Subsequently the following scripts were used for further analyses and visualizations: 

- Figure 1
  - `002_plot_n_m6A_sites.R`
  - `021_sites_per_gene.R`
  - `020_homer_results.R`
  - `004_plot_mutation_in_methylated_motifs.R`
  - `005_plot_proportion_of_DRACH.R`

- Figure 2
  - `006_plot_recurrence_mutated_motifs.R`
  - `008_top_mutated_m6A_genes.R`
  - `009_top_mutated_m6A_motifs.R`

- Figure 3
  - `012_expression_DESeq2_gene_level.R`
  - `013_expression_DESeq2_motif_level.R`
  - `014_CTNNB1_gene_model.R`

- Figure 4
  - `015_m6AMut_vs_totalMut.R`
  - `016_fraction_mutation_in_meth_site_cancer_type.R`
  - `017_survival_quantiles_cancer_types.R`
  - `018_KM_quartiles_sig.R`
  - `019_total_mut_comparison.R`

- Figure 5
  - `002_motif_consequence.R`
  - `003_gain_per_cancer_type.R`
  - `005_recurrence_gain.R`
  - `006_survival_analysis.R`
  - `015_DESeq2_results_genes.R`

- Supplemental Figure 1
  - `010_methylated_genes_per_cancer_type_above_threshold_with_loss.R`
  - `011_methylated_motifs_per_cancer_type_above_threshold_with_loss.R`
  - `022_sites_per_cancer_gene.R`
  - `022b_sites_per_cancer_gene_per_kb.R`

- Supplemental Figure 2
  - `023_mutated_sites_per_gene_patient.R`
  - `007_plot_compare_methylation_cancer_genes.R`

- Supplemental Figure 3
  - `007_KM_curves.R`
  - `011_DESeq2_results_motifs.R`
  - `016_recur_sig_genes.R`
  - `017_gain_per_gene_per_patient.R`