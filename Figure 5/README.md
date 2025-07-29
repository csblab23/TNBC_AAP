# Construction of Relapse-Associated Risk Assessment Model

This repository contains scripts and data used to construct predictive models for relapse risk assessment in patients.

## Key Components

**1.Panel05/FS/**: Directory containing all model construction scripts organized by individual feature sets (FS = Feature Set). <br>
**2.Input.csv**: Unified input dataset used across all models for training and evaluation. <br>
**3.KM_plots_for_risk_models_training_Set.R**: Code for generating the KM plots on training set as presented in figure 5a-c.<br>
**4.KM_plots_for_risk_models_validation_Set.R**: Code for generating the KM plots on validation set as presented in figure 5d-f.<br>
**5.script_for_metric_calculation.R**: Code for generating the line plots for C-index; IBS and AUROC as presented in figure g-i.<br>

Each feature set is evaluated independently to assess its contribution to relapse prediction. Models may be constructed using:
- Clinical features only
- Promoter features only
- Combined (clinical + promoter) features

##  Model Objectives

- Identify features associated with relapse risk
- Compare performance across different feature sets

##  Workflow Summary

1. Preprocess input data (`Input.csv`)
2. Select relevant features (FS)
3. Train and evaluate models
4. Store outputs in structured folders

---


