# Construction of Relapse-Associated Risk Assessment Model

This repository contains scripts and data used to construct predictive models for relapse risk assessment in patients.

## 📁 Structure and Key Components

- **Panel05/FS/** – Directory containing all model construction scripts organized by individual feature sets (FS = Feature Set).
- **Input.csv** – Unified input dataset used across all models for training and evaluation.

Each feature set is evaluated independently to assess its contribution to relapse prediction. Models may be constructed using:
- Clinical features only
- Promoter features only
- Combined (clinical + promoter) features

## 🛠️ Model Objectives

- Identify features associated with relapse risk
- Compare performance across different feature sets

## 🔄 Workflow Summary

1. Preprocess input data (`Input.csv`)
2. Select relevant features (FS)
3. Train and evaluate models
4. Store outputs in structured folders

---


