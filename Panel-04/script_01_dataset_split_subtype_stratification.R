## Code to stratify the FUSCC dataset into Training & Validation (testing) sets based on the mRNA_subtype ##

# Load required packages
library(caret)
library(survival)
library(survminer)
library(proActiv)
library(data.table)
library(dplyr)
library(tidyverse)
library(readxl)
library(limma)
library(edgeR)
library(sva)

set.seed(123)

# Load metafile
metadata <- as.data.frame(read_csv("meta_with_batchinfo.csv"))
metadata <- metadata[, -1]  # Remove extra column
metadata <- metadata %>%
  mutate(Tumor_Normal1 = case_when(
    Tumor_Normal %in% c("TumorTissue", "TumorTissue.2018", "TumorTissue.rep") ~ "Tumor",
    Tumor_Normal %in% c("PairedNormalTissue", "PairedNormalTissue.rep") ~ "AdjNormal"
  ))
rownames(metadata) <- metadata$Run
tumor_metadata <- subset(metadata, Tumor_Normal1 == 'Tumor') # taking only tumor samples

# Take samples with RFS>0
surv_metadata =tumor_metadata[tumor_metadata$RFS_time_Months > 0,] 
surv_metadata=surv_metadata[,c("Project_ID", "Run", "batch", "Tumor_Normal", "Tumor_Normal1", "RFS_Status", "RFS_time_Days", "RFS_time_Months", 'CNA_Subtype', 'Grade', 'mRNA_Subtype')]
filter_data <- surv_metadata %>%
  mutate(CNA_Subtype = replace_na(CNA_Subtype, 'unknown'),
         Grade = replace_na(Grade, 'unknown'))
filter_data$CNA_Subtype = as.factor(filter_data$CNA_Subtype)
filter_data$Grade = as.factor(filter_data$Grade)
filter_data$mRNA_Subtype = as.factor(filter_data$mRNA_Subtype)

# Divide the data in train-test using caret package
inTrain <- createDataPartition(
  y = filter_data$mRNA_Subtype,
  p = .70, # Proportion for the training set
  list = FALSE
)

training <- filter_data[ inTrain,]
testing  <- filter_data[-inTrain,]

# Save files
write.csv(x = training, file = 'training_dataset_sub.csv')
write.csv(x = testing, file = 'testing_dataset_sub.csv')

##############################################################################
##############################################################################