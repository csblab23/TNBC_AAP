# script - Feature set - f (6 feature model) ---- how feature selection take place! 

library(tidyverse)  
library(caret)      
library(leaps)      
library(MASS)

#reading input file
final <- as.data.frame(read_csv("./stepwise_regression_14_features_training.csv"))
final <- as.data.frame(final)
genes <- colnames(final)[11:dim(final)[2]]
cox_model <- coxph(as.formula(paste("Surv(RFS_time_Months, RFS_Status) ~", paste(genes, collapse = " + "))), data = final)

step.model <- stepAIC(cox_model, direction = "both", trace = TRUE)
summary(step.model)
# Extract selected genes
selected_genes <- names(coef(step.model))
coef_25 <- as.data.frame(coef(step.model))
coef_25$gene <- rownames(coef_25)
print(selected_genes)
length(selected_genes)

backward_gene_list = selected_genes
length(backward_gene_list) #25
