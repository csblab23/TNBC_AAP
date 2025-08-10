## Code to make the forest plot given in Figure 4g ##

# Load required packages
library(survminer)
library(proActiv)
library(data.table)
library(dplyr)
library(tidyverse)
library(readxl)
library(limma)
library(edgeR)
library(sva)
library(ggforestplot)
library(ggplot2)
library(broom)
library(survival)
library(survivalAnalysis)
library(fastDummies)

# Load the metadata
metadata <- as.data.frame(read_csv("meta_with_batchinfo.csv"))
metadata <- metadata[, -1]  # Remove extra column
metadata <- metadata %>%
  mutate(Tumor_Normal1 = case_when(
    Tumor_Normal %in% c("TumorTissue", "TumorTissue.2018", "TumorTissue.rep") ~ "Tumor",
    Tumor_Normal %in% c("PairedNormalTissue", "PairedNormalTissue.rep") ~ "AdjNormal"
  ))
rownames(metadata) <- metadata$Run

# Take only tumor samples and subset required columns
tumor_metadata <- subset(metadata, Tumor_Normal1 == 'Tumor')
tumor_metadata <- tumor_metadata[,c('Project_ID', 'Run', 'batch', 'Tumor_Normal1', 'mRNA_Subtype', 'CNA_Subtype', 'Size_cm', 'RFS_Status', 'RFS_time_Months', 'Ki67')]
tumor_metadata <- tumor_metadata %>%
  mutate(CNA_Subtype = replace_na(CNA_Subtype, 'unknown'))

# Load the result2 file obtained by running proActiv on only tumor samples 
result2 <- readRDS("result2_v2_updated.rds")
result2_promoter = as.data.frame(rowData(result2)) 

# Load absolute promoter activity
Ab_counts <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)
Ab_sub <- Ab_counts[c( '3004',  '41492', '9001', '12191', '19936'), ]  # fetching our promoters of interest
colnames(Ab_sub) <- gsub("_SJ.out", "", colnames(Ab_sub))
Ab_sub <- as.data.frame(t(Ab_sub))

# Match order of samples
tumor_metadata <- tumor_metadata[rownames(Ab_sub),]
summary(tumor_metadata)
tumor_metadata$mRNA_Subtype <- as.factor(tumor_metadata$mRNA_Subtype)
tumor_metadata$CNA_Subtype <- as.factor(tumor_metadata$CNA_Subtype)
summary(tumor_metadata)
tumor_metadata$pr3004_HUWE1 <- Ab_sub$'3004'
tumor_metadata$pr41492_FTX <- Ab_sub$'41492'
tumor_metadata$pr9001_CALD1 <- Ab_sub$'9001'
tumor_metadata$pr12191_TANK <- Ab_sub$'12191'
tumor_metadata$pr19936_TPM4 <- Ab_sub$'19936'
tumor_metadata <- na.omit(tumor_metadata)

tumor_metadata$pr_exp_3004_HUWE1 <- ifelse(tumor_metadata$pr3004_HUWE1 > 3.913305, "1", "0")  # cutoff extracted from the Rdata for promoters 
tumor_metadata$pr_exp_41492_FTX <- ifelse(tumor_metadata$pr41492_FTX > 3.427995, "1", "0") 

tumor_metadata$pr_exp_9001_CALD1 <- ifelse(tumor_metadata$pr9001_CALD1 > 4.112375, "1", "0") 
tumor_metadata$pr_exp_12191_TANK <- ifelse(tumor_metadata$pr12191_TANK > 4.336671, "1", "0") 
tumor_metadata$pr_exp_19936_TPM4 <- ifelse(tumor_metadata$pr19936_TPM4 > 6.645274, "1", "0") 

str(tumor_metadata)
met_sub <- tumor_metadata[,c(1:10, 16:dim(tumor_metadata)[2])]
str(met_sub)
met_sub[,11:dim(met_sub)[2]] <- lapply(met_sub[,11:dim(met_sub)[2]], as.factor)

# Subset columns of interest
df_final <- met_sub[,c("Project_ID","Run","RFS_Status", "RFS_time_Months", "pr_exp_3004_HUWE1", "pr_exp_41492_FTX","pr_exp_9001_CALD1", "pr_exp_12191_TANK", "pr_exp_19936_TPM4", 'CNA_Subtype',  'mRNA_Subtype')]
df_final$mRNA_Subtype <- as.factor(df_final$mRNA_Subtype)

# Remove NA's from the dataframe
cleaned_df_final <- na.omit(df_final)
cleaned_df_final$CNA_Subtype <- as.factor(cleaned_df_final$CNA_Subtype) # convert to factor
str(cleaned_df_final)
cleaned_df_final[,5:dim(cleaned_df_final)[2]] <- lapply(cleaned_df_final[,5:dim(cleaned_df_final)[2]], as.character)
cleaned_df_final[,5:dim(cleaned_df_final)[2]] <- lapply(cleaned_df_final[,5:dim(cleaned_df_final)[2]], as.factor)

# Preparing an input dataframe with all required information
input_df <- data.frame(Sample = cleaned_df_final$Run, time = cleaned_df_final$RFS_time_Months, event = cleaned_df_final$RFS_Status,
                       pr3004_HUWE1 = cleaned_df_final$pr_exp_3004_HUWE1,  
                       pr41492_FTX = cleaned_df_final$pr_exp_41492_FTX,
                       
                       pr9001_CALD1 = cleaned_df_final$pr_exp_9001_CALD1,
                       pr12191_TANK = cleaned_df_final$pr_exp_12191_TANK,
                       pr19936_TPM4 = cleaned_df_final$pr_exp_19936_TPM4,
                       
                       CNA_Subtype = cleaned_df_final$CNA_Subtype,
                       mRNA_Subtype = cleaned_df_final$mRNA_Subtype
)

input_df_2 <- input_df %>%
  mutate(pr41492_FTX = case_when(
    pr41492_FTX == '0' ~ 'Low_pr41492_FTX',
    pr41492_FTX == '1' ~ 'High_pr41492_FTX'
  ))

input_df_2 <- input_df_2 %>%
  mutate(pr3004_HUWE1 = case_when(
    pr3004_HUWE1 == '0' ~ 'Low_pr3004_HUWE1',
    pr3004_HUWE1 == '1' ~ 'High_pr3004_HUWE1'
  ))

input_df_2 <- input_df_2 %>%
  mutate(pr9001_CALD1 = case_when(
    pr9001_CALD1 == '0' ~ 'Low_pr9001_CALD1',
    pr9001_CALD1 == '1' ~ 'High_pr9001_CALD1'
  ))

input_df_2 <- input_df_2 %>%
  mutate(pr12191_TANK = case_when(
    pr12191_TANK == '0' ~ 'Low_pr12191_TANK',
    pr12191_TANK == '1' ~ 'High_pr12191_TANK'
  ))

input_df_2 <- input_df_2 %>%
  mutate(pr19936_TPM4 = case_when(
    pr19936_TPM4 == '0' ~ 'Low_pr19936_TPM4',
    pr19936_TPM4 == '1' ~ 'High_pr19936_TPM4'
  ))

# Converting to factor and defining reference
input_df_2$pr3004_HUWE1 <- as.factor(input_df_2$pr3004_HUWE1)
levels(input_df_2$pr3004_HUWE1)
input_df_2$pr41492_FTX <- as.factor(input_df_2$pr41492_FTX)
levels(input_df_2$pr41492_FTX)
input_df_2$pr41492_FTX <- relevel(input_df_2$pr41492_FTX, ref = 'Low_pr41492_FTX')
input_df_2$mRNA_Subtype <- relevel(input_df_2$mRNA_Subtype, ref = 'IM')
input_df_2$CNA_Subtype <- relevel(input_df_2$CNA_Subtype, ref = 'Low_CIN')

# Cox model specification with factor
cox_model <- coxph(Surv(time = time, event = event) ~ pr3004_HUWE1 + pr41492_FTX + pr9001_CALD1 +pr12191_TANK+ pr19936_TPM4  + mRNA_Subtype + CNA_Subtype, data = input_df_2)
# Summarize the Cox model
summary(cox_model)

# Visualization with ggforest
forest_plot <- ggforest(model = cox_model, 
                        data = input_df_2,  # Add this line
                        main = '',
                        fontsize = 1,
                        
                        refLabel = 'Reference')

p <- forest_plot +
  theme(
    text = element_text(family = "Arial", face = 'bold', size = 18)
  )
p
ggsave("forestplot.png", p, width = 13, height = 12, dpi = 1200)

############################################################################
############################################################################
