# **Identification of TNBC prognostically significant Active Alternative Promoters (AAPs)**

This repository contains scripts and data used to identify BL1-specific and Basal subtype-specific AAPs in the FUSCC cohort and their subsequent validation in the external dataset 1.

**Key Components**

1. *script_01_dataset_split_subtype_stratification.R*: Code to stratify the FUSCC TNBC samples into training and validation set.
2. *script_02_gene_cutoff_univariate_training.R*: Code for calculating maximally selected rank statistics cutoff value to obtain 5,000 iterations for each gene individually in the training set; performing univariate survival analysis & identification of the nonsignificant genes.
3. *script_03_prom_cutoff_univariate_training.R*: Code for calculating maximally selected rank statistics cutoff value to obtain 5,000 iterations for each promoter individually in the training set; performing univariate survival analysis & identification of the significant promoters.
4. *script_04_gene_univariate_testing.R*: Code for running univariate survival analysis on genes in testing set by applying the cutoff calculated from the training set and fetching nonsignificant genes.
5. *script_05_prom_univariate_testing.R*: Code for running univariate survival analysis on promoters in testing set by applying the cutoff calculated from the training set and fetching significant promoters and subsequent validation of the significant prognostic AAPs.
6. *script_07_km_plots_for_candidates.R*: Code to generate Kaplan–Meier (KM) plots for *HUWE1* & *FTX* and their corresponding promoters.
