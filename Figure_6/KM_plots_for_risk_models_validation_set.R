## Code for Figure 5 d,e,f -->> To construct risk groups and to check the survival of each risk group in Validation set for: Clinical only; Promoter only & All features ##

# Load required packages
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
library(maxstat)
library(fastDummies)

# Read the metafile
metadata <- as.data.frame(read_csv("meta_with_batchinfo.csv"))
metadata <- metadata[, -1]  # Remove extra column
metadata <- metadata %>%
  mutate(Tumor_Normal1 = case_when(
    Tumor_Normal %in% c("TumorTissue", "TumorTissue.2018", "TumorTissue.rep") ~ "Tumor",
    Tumor_Normal %in% c("PairedNormalTissue", "PairedNormalTissue.rep") ~ "AdjNormal"
  ))
rownames(metadata) <- metadata$Run

# Subset tumor samples
tumor_metadata <- subset(metadata, Tumor_Normal1 == 'Tumor')
tumor_metadata <- tumor_metadata[,c('Project_ID', 'Run', 'batch', 'Tumor_Normal1', 'mRNA_Subtype', 'CNA_Subtype', 'Size_cm', 'RFS_Status', 'RFS_time_Months', 'Ki67')]
tumor_metadata <- tumor_metadata %>%
  mutate(CNA_Subtype = replace_na(CNA_Subtype, 'unknown'))  # replace NA with unknown

# Integrate the expression of the promoters 
result2 <- readRDS("result2_v2_updated.rds")
result2_promoter = as.data.frame(rowData(result2)) 
Ab_counts <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)
Ab_sub <- Ab_counts[c( '3004',  '41492'), ]
colnames(Ab_sub) <- gsub("_SJ.out", "", colnames(Ab_sub))
Ab_sub <- as.data.frame(t(Ab_sub))

# Make sure sample order is same 
tumor_metadata <- tumor_metadata[rownames(Ab_sub),]
tumor_metadata$pr_3004_HUWE1 <- Ab_sub$'3004'
tumor_metadata$pr_41492_FTX <- Ab_sub$'41492'

# Converting to binary based on the cut off value
tumor_metadata$pr_exp_3004_HUWE1 <- ifelse(tumor_metadata$pr_3004_HUWE1 > 3.913305, "0", "1")
tumor_metadata$pr_exp_41492_FTX <- ifelse(tumor_metadata$pr_41492_FTX > 3.427995, "1", "0")
met_sub <- tumor_metadata[,c(1:10, 13:14)]
met_sub[,11:12] <- lapply(met_sub[,11:12], as.factor)

# Here, keep only the samples that are present in validation (testing) dataset
val_dataset_sub <- read_csv("testing_dataset_sub.csv")
met_sub <- met_sub[val_dataset_sub$Run,] 

## Find maxstat cut-off for the clinical features
# For tumor size
mod_maxstat <- maxstat.test(Surv(RFS_time_Months, RFS_Status) ~ met_sub$Size_cm,
                            data = met_sub, smethod = "LogRank", pmethod = c("none"))
cutoff_values<- mod_maxstat$estimate
met_sub$Size_cm <- ifelse(met_sub$Size_cm > 3.3, "1", "0")
met_sub$Size_cm <- as.factor(met_sub$Size_cm)

##--

# For Ki67 index
mod_maxstat <- maxstat.test(Surv(RFS_time_Months, RFS_Status) ~ met_sub$Ki67,
                            data = met_sub, smethod = "LogRank", pmethod = c("none"))
cutoff_values<- mod_maxstat$estimate
met_sub$Ki67 <- ifelse(met_sub$Ki67 > 20, "1", "0")
met_sub$Ki67 <- as.factor(met_sub$Ki67)

# Subset required columns
df_final <- met_sub[,c("Project_ID","Run","RFS_Status", "RFS_time_Months", "Size_cm", "Ki67", "pr_exp_3004_HUWE1", "pr_exp_41492_FTX")]

# Take complete cases
df_final <- df_final[complete.cases(df_final), ]
df_final[,5:6] <- lapply(df_final[,5:6], as.character)
df_final[,5:6] <- lapply(df_final[,5:6], as.numeric)

# Find number of patients in each risk group
df_final$risk_group <- rowSums(df_final[,5:6], na.rm = TRUE)

# KM plot 
surv_meta <- df_final[,c(1:4, 9)]
surv_meta$risk_group <- as.factor(surv_meta$risk_group)

## Figure 5d ------- Plot for clinical features only ------------------------------------------------------->>
sfit <- survfit(Surv(surv_meta$RFS_time_Months, surv_meta$RFS_Status)~ surv_meta[,'risk_group'], data=surv_meta)
test <- ggsurvplot(
  sfit,
  pval = TRUE,
  pval.method = TRUE,
  size = 1.5,
  palette = c("purple", "magenta", "forestgreen"),
  title = "Clinical Features Only (Validation)",
  xlab = "Time (Months)",
  ylab = "Relapse Free Survival",
  legend.labs = c("0ev. (n = 15)", "1ev. (n = 76)", "2ev. (n = 9)"),
  legend.title = "Risk Events",
  legend = c(0.75, 0.25),  # Adjust legend position inside the plot area
  risk.table = TRUE,
  risk.table.height = 0.2,
  surv.median.line = "hv",
  risk.table.fontsize = 6,
  fontsize = 16,
  pval.size = 15,
  font.main = 20,
  font.tickslab = 20,
  color.tickslab = 'black',
  size = 1,
  tables.y.text = FALSE,
  font.legend = 20,
  ggtheme = theme_minimal() +
    theme(
      plot.title = element_text(size = 24, face = 'bold', family = 'Arial', hjust = 0.5),
      axis.title.x = element_text(size = 18, face = 'bold', colour = 'black'),
      axis.title.y = element_text(size = 18, face = 'bold', colour = 'black'),
      axis.text.x = element_text(colour = 'black'),
      axis.text.y = element_text(colour = 'black'),
      legend.text = element_text(size = 20, face = 'bold'),
      legend.title = element_text(size = 25, face = 'bold'),
      legend.position = c(0.75, 0.25),  # Position legend inside the plot
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_line(colour = "grey90"),
      panel.grid.minor = element_line(colour = "grey90"),
      panel.border = element_rect(color = "black", size = 1, fill = NA),
      panel.background = element_blank()
    )
)
test
ggsave("24_12_KM_for_clin_features_val.png", test$plot, width = 7, height = 6, dpi = 1200,bg="white")


## Figure 5e ------- Plot for prognostic features (promoters) only ------------------------------------------------------->>
df_final[,7:8] <- lapply(df_final[,7:8], as.character)
df_final[,7:8] <- lapply(df_final[,7:8], as.numeric)

# Find number of patients in each risk group
df_final$risk_group <- rowSums(df_final[,7:8], na.rm = TRUE)

# KM plot 
surv_meta <- df_final[,c(1:4, 9)]
surv_meta$risk_group <- as.factor(surv_meta$risk_group)

sfit <- survfit(Surv(surv_meta$RFS_time_Months, surv_meta$RFS_Status)~ surv_meta[,'risk_group'], data=surv_meta)
test <- ggsurvplot(
  sfit,
  pval = TRUE,
  pval.method = TRUE,
  size = 1.5,
  palette = c("purple", "magenta", "forestgreen"),
  title = "Promoter Features Only (Validation)",
  xlab = "Time (Months)",
  ylab = "Relapse Free Survival",
  legend.labs = c("0ev. (n = 41)", "1ev. (n = 35)", "2ev. (n = 24)"),
  legend.title = "Risk Events",
  legend = c(0.75, 0.15),  # Adjust legend position inside the plot area
  risk.table = TRUE,
  risk.table.height = 0.2,
  surv.median.line = "hv",
  risk.table.fontsize = 6,
  fontsize = 16,
  pval.size = 12,
  font.main = 20,
  font.tickslab = 20,
  color.tickslab = 'black',
  size = 1,
  tables.y.text = FALSE,
  font.legend = 20,
  ggtheme = theme_minimal() +
    theme(
      plot.title = element_text(size = 24, face = 'bold', family = 'Arial', hjust = 0.5),
      axis.title.x = element_text(size = 18, face = 'bold', colour = 'black'),
      axis.title.y = element_text(size = 18, face = 'bold', colour = 'black'),
      axis.text.x = element_text(colour = 'black'),
      axis.text.y = element_text(colour = 'black'),
      legend.text = element_text(size = 20, face = 'bold'),
      legend.title = element_text(size = 25, face = 'bold'),
      legend.position = c(0.70, 0.15),  # Position legend inside the plot
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_line(colour = "grey90"),
      panel.grid.minor = element_line(colour = "grey90"),
      panel.border = element_rect(color = "black", size = 1, fill = NA),
      panel.background = element_blank()
    )
)
test
ggsave("24_12_KM_for_prom_features_val.png", test$plot, width = 7, height = 6, dpi = 1200,bg="white")


## Figure 5f ------- Plot for all features (clinical + promoters) ------------------------------------------------------->>
# Find number of patients in each risk group
df_final$risk_group <- rowSums(df_final[,5:8], na.rm = TRUE)
df_final$risk_group2 <- cut(df_final$risk_group,
                            breaks = c(-Inf, 1, 2, 3, Inf),  # Adjust breaks for your desired intervals
                            labels = c("0-1", "2", "3+", "3+"),  # Apply the labels accordingly
                            include.lowest = TRUE,
                            right = TRUE)  # Right-closed intervals
# KM plot 
df_final$risk_group <- as.factor(df_final$risk_group)
df_final$risk_group2 <- as.factor(df_final$risk_group2)
surv_meta <- df_final[,c(1:4, 10)]

sfit <- survfit(Surv(surv_meta$RFS_time_Months, surv_meta$RFS_Status)~ surv_meta[,'risk_group2'], data=surv_meta)
test <- ggsurvplot(
  sfit,
  pval = TRUE,
  pval.method = TRUE,
  size = 1.5,
  palette = c("purple", "magenta", "forestgreen"),
  title = "All Features (Validation)",
  xlab = "Time (Months)",
  ylab = "Relapse Free Survival",
  legend.labs = c("0-1ev. (n = 44)", "2ev. (n = 32)", "3+ev. (n = 24)"),
  legend.title = "Risk Events",
  legend = c(0.75, 0.25),  # Adjust legend position inside the plot area
  risk.table = TRUE,
  risk.table.height = 0.2,
  surv.median.line = "hv",
  risk.table.fontsize = 6,
  fontsize = 16,
  pval.size = 13,
  font.main = 20,
  font.tickslab = 20,
  color.tickslab = 'black',
  size = 1,
  tables.y.text = FALSE,
  font.legend = 20,
  ggtheme = theme_minimal() +
    theme(
      plot.title = element_text(size = 24, face = 'bold', family = 'Arial', hjust = 0.5),
      axis.title.x = element_text(size = 18, face = 'bold', colour = 'black'),
      axis.title.y = element_text(size = 18, face = 'bold', colour = 'black'),
      axis.text.x = element_text(colour = 'black'),
      axis.text.y = element_text(colour = 'black'),
      legend.text = element_text(size = 20, face = 'bold'),
      legend.title = element_text(size = 25, face = 'bold'),
      legend.position = c(0.75, 0.25),  # Position legend inside the plot
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_line(colour = "grey90"),
      panel.grid.minor = element_line(colour = "grey90"),
      panel.border = element_rect(color = "black", size = 1, fill = NA),
      panel.background = element_blank()
    )
)
test

ggsave("24_12_KM_for_all_features_val.png", test$plot, width = 7, height = 6, dpi = 1200,bg="white")

#######################################################################################################
#######################################################################################################
