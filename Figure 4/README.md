# **Identification of TNBC prognostically significant Active Alternative Promoters (AAPs)**

This repository contains scripts and data used to identify BL1-specific and Basal subtype-specific AAPs in the FUSCC cohort and their subsequent validation in the external dataset 1.

**Key Components**

1. script_01_dataset_split_subtype_stratification.R: Code to stratify the FUSCC TNBC samples into training and validation set.
2. script_02_gene_cutoff_univariate_training.R: Code for calculating maximally selected rank statistics cutoff value by running to obtain 5000 iterations for each gene individually in the training set; performing univariate survival analysis & identification of the nonsignificant genes.
3. script_03_prom_cutoff_univariate_training.R: Code for calculating maximally selected rank statistics cutoff value by running to obtain 5000 iterations for each promoter individually in the training set; performing univariate survival analysis & identification of the significant promoters.
4. script_04_Maxstat_cutoff_boxplots.R: Code for generating the boxplot of all 5000 cutoffs obtained for the validated candidate promoters and their corresponding genes.
5. script_05_gene_univariate_testing.R: Code for running univariate survival analysis on genes in testing set by applying the cutoff calculated from the training set and fetching nonsignificant genes.
6. script_06_prom_univariate_testing.R: Code for running univariate survival analysis on promoters in testing set by applying the cutoff calculated from the training set and fetching significant promoters and subsequent validation of the significant prognostic AAPs.
7. script_07_km_plots_for_candidates.R: Code to generate Kaplan–Meier (KM) plots for *HUWE1* & *FTX* and their corresponding promoters.
8. script_08_forest_plot.R: Code to reproduce the forest plot presented in Figure 4g of the paper.
9. script_09_external_val_only_TNBC_samples.R: Code for validation of *HUWE1* & *FTX* gene expression and promoter activity variations between high risk and low risk TNBC patients of the external dataset 3.
10. script_10_external_val_all_BC_samples.R: Code for validation of *HUWE1* & *FTX* gene expression and promoter activity variations between high risk and low risk breast cancer patients of the external dataset 3.
