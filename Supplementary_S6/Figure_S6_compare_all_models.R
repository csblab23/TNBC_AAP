## Code to make line plots shown in Supplementary figure S6 for comparing the calculated metrics across all the generated models ##

# Load required libraries
library(ggplot2)
library(dplyr)
library(readxl)
library(tidyverse)

set.seed(176092.01)

# Load the dataframe having C-index values for each model:
training_cindex_df <- read_xlsx("all_models_Cindex_training_testing.xlsx")
training_cindex_df <- as.data.frame(training_cindex_df)
training_cindex_df$Feature_Set <- paste0("FS_",training_cindex_df$Feature_Set)

# Round to two digits
training_cindex_df$C_index <- round(training_cindex_df$C_index, 2)

# Define feature order
feature_order <- c("FS_C","FS_A","FS_F","FS_B","FS_E","FS_D")

# Apply order to Feature_Set
training_cindex_df$Feature_Set <- factor(training_cindex_df$Feature_Set, levels = feature_order)
training_cindex_df$Submodel <- factor(training_cindex_df$Submodel,levels = c("All features",
                                                                             "Promoter features only",
                                                                             "Clinical features only"))

# Plot
p <- ggplot(training_cindex_df,
            aes(x = Feature_Set, y = C_index,
                group = interaction(Submodel, Set),
                color = Submodel, linetype = Set)) +
  geom_point(size = 3) +
  geom_line() +
  facet_wrap(~ Method) +
  scale_linetype_manual(values = c("Training" = "dotted", "Testing" = "solid")) +
  labs(
    title = "Comparison of C-index across models",
    x = "Feature Set",
    y = "C-index",
    color = "Submodel",
    linetype = "Data Split"
  ) +
  theme_bw() +
  scale_y_continuous(breaks = pretty(training_cindex_df$C_index, n = 10),limits = c(0.60,0.76)) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1,colour = "black",size=10),
    axis.text.y = element_text(colour = "black",size=12),
    axis.title = element_text(colour = "black",size=15,face = "bold"),
    text = element_text(colour = "black"),
    plot.title = element_text(hjust = 0.5,face = "bold",colour = "black",size = 15),
    strip.text = element_text(colour = "black", face = "bold", size = 12)
  )
p
ggsave("cindex_FINAL_plot.png", plot = p, width = 9, height = 4, dpi = 1000)

#################################################################################

# Load the dataframe having IBS values for each model:
df <- read_xlsx("all_models_IBS_training_testing.xlsx")
df <- as.data.frame(df)
df$Feature_Set <- paste0("FS_",df$Feature_Set)

# Round to two digits
df$IBS <- round(df$IBS, 3)

# Define feature order 
feature_order <- c("FS_C","FS_A","FS_F","FS_B","FS_E","FS_D")

# Apply order to Feature_Set
df$Feature_Set <- factor(df$Feature_Set, levels = feature_order)
df$Submodel <- factor(df$Submodel,levels = c("All features",
                                             "Promoter features only",
                                             "Clinical features only"))
# Plot
p <- ggplot(df,
            aes(x = Feature_Set, y = IBS,
                group = interaction(Submodel, Set),
                color = Submodel, linetype = Set)) +
  geom_point(size = 3) +
  geom_line() +
  facet_wrap(~ Method) +
  scale_linetype_manual(values = c("Training" = "dotted", "Testing" = "solid")) +
  labs(
    title = "Comparison of IBS across models",
    x = "Feature Set",
    y = "IBS",
    color = "Submodel",
    linetype = "Data Split"
  ) +
  theme_bw() +
  # scale_y_continuous(breaks = pretty(df$IBS, n = 10),limits = c(0.60,0.76)) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1,colour = "black",size=10),
    axis.text.y = element_text(colour = "black",size=12),
    axis.title = element_text(colour = "black",size=15,face = "bold"),
    text = element_text(colour = "black"),
    plot.title = element_text(hjust = 0.5,face = "bold",colour = "black",size = 15),
    strip.text = element_text(colour = "black", face = "bold", size = 12)
  )
p
ggsave("ibs_FINAL_plot.png", plot = p, width = 9, height = 4, dpi = 1000)

#################################################################################

# Load the dataframe having AUROC values for each model:
df <- read_xlsx("all_models_AUROC_training_testing.xlsx")
df <- as.data.frame(df)
df$Feature_Set <- paste0("FS_",df$Feature_Set)

# Round to two digits
df$AUROC <- round(df$AUROC, 2)

# Define feature order 
feature_order <- c("FS_C","FS_A","FS_F","FS_B","FS_E","FS_D")

# Apply order to Feature_Set
df$Feature_Set <- factor(df$Feature_Set, levels = feature_order)
df$Submodel <- factor(df$Submodel,levels = c("All features",
                                             "Promoter features only",
                                             "Clinical features only"))

# Plot
p <- ggplot(df,
            aes(x = Feature_Set, y = AUROC,
                group = interaction(Submodel, Set),
                color = Submodel, linetype = Set)) +
  geom_point(size = 3) +
  geom_line() +
  facet_wrap(~ Method) +
  scale_linetype_manual(values = c("Training" = "dotted", "Testing" = "solid")) +
  labs(
    title = "Comparison of AUROC across models",
    x = "Feature Set",
    y = "AUROC",
    color = "Submodel",
    linetype = "Data Split"
  ) +
  theme_bw() +
  scale_y_continuous(breaks = pretty(df$AUROC, n = 10),limits = c(0.58,0.76)) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1,colour = "black",size=10),
    axis.text.y = element_text(colour = "black",size=12),
    axis.title = element_text(colour = "black",size=15,face = "bold"),
    text = element_text(colour = "black"),
    plot.title = element_text(hjust = 0.5,face = "bold",colour = "black",size = 15),
    strip.text = element_text(colour = "black", face = "bold", size = 12)
  )
p
ggsave("auroc_FINAL_plot.png", plot = p, width = 9, height = 4, dpi = 1000)
######################################################################################
######################################################################################
