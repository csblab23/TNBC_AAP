## This script includes the code used for generating C-index; IBS score and AUROC plots shown in Figure 5g,h,i ##

# Load required libraries
library(dplyr)
library(tidyverse)

## 1. For C-index Figure 5g ------------------------------------------>>

# Read the files having scores:
all_features <- read_csv("clin_and_prom.csv")
prom_only <- read_csv("prom_only.csv")
clin_only <- read_csv("clin_only.csv")

# Omit-NA values
all_features <- na.omit(all_features)
prom_only <- na.omit(prom_only)
clin_only <- na.omit(clin_only)

# Calculate the standard deviation of performance on test datasets at different seeds
# Random forest
test_rf_all_sd <- sd(all_features$c_index_rf_test)
test_rf_prom_sd <- sd(prom_only$c_index_rf_test)
test_rf_clin_sd <- sd(clin_only$c_index_rf_test)

# Lasso
test_lasso_all_sd <- sd(all_features$c_index_lasso_test)
test_lasso_prom_sd <- sd(prom_only$c_index_lasso_test)
test_lasso_clin_sd <- sd(clin_only$c_index_lasso_test)

# Elastic net
test_elastic_all_sd <- sd(all_features$c_index_elastic_test)
test_elastic_prom_sd <- sd(prom_only$c_index_elastic_test)
test_elastic_clin_sd <- sd(clin_only$c_index_elastic_test)


# Calculate the mean performance on test dataset at different seeds:
# Random forest
test_rf_all <- mean(all_features$c_index_rf_test)
test_rf_prom <- mean(prom_only$c_index_rf_test)
test_rf_clin <- mean(clin_only$c_index_rf_test)

# Lasso
test_lasso_all <- mean(all_features$c_index_lasso_test)
test_lasso_prom <- mean(prom_only$c_index_lasso_test)
test_lasso_clin <- mean(clin_only$c_index_lasso_test)

# Elastic net
test_elastic_all <- mean(all_features$c_index_elastic_test)
test_elastic_prom <- mean(prom_only$c_index_elastic_test)
test_elastic_clin <- mean(clin_only$c_index_elastic_test)


# Calculate mean of standard deviation on train dataset at different seeds:
# Random forest
sd_tr_rf_all <- mean(all_features$std_c_index_rf)
sd_tr_rf_prom <- mean(prom_only$std_c_index_rf)
sd_tr_rf_clin <- mean(clin_only$std_c_index_rf)

# Lasso
sd_tr_lasso_all <- mean(all_features$std_c_index_lasso)
sd_tr_lasso_prom <- mean(prom_only$std_c_index_lasso)
sd_tr_lasso_clin <- mean(clin_only$std_c_index_lasso)

# Elastic net
sd_tr_elastic_all <- mean(all_features$std_c_index_elastic)
sd_tr_elastic_prom <- mean(prom_only$std_c_index_elastic)
sd_tr_elastic_clin <- mean(clin_only$std_c_index_elastic)

# Calculate mean of mean performance on train dataset at different seeds:
# Random forest
mean_training_c_index_rf_all <- mean(all_features$mean_c_index_rf)
mean_training_c_index_rf_prom <- mean(prom_only$mean_c_index_rf)
mean_training_c_index_rf_clin <- mean(clin_only$mean_c_index_rf)

# Elastic net
mean_training_c_index_elastic_all <- mean(all_features$mean_c_index_elastic)
mean_training_c_index_elastic_prom <- mean(prom_only$mean_c_index_elastic)
mean_training_c_index_elastic_clin <- mean(clin_only$mean_c_index_elastic)

# Lasso
mean_training_c_index_lasso_all <- mean(all_features$mean_c_index_lasso)
mean_training_c_index_lasso_prom <- mean(prom_only$mean_c_index_lasso)
mean_training_c_index_lasso_clin <- mean(clin_only$mean_c_index_lasso)

# Summarising all information in a dataframe
data <- data.frame(
  Method = c("Lasso", "Lasso", "Lasso", "Elastic", "Elastic", "Elastic", "Random_Survival\n_Forest", "Random_Survival\n_Forest", "Random_Survival\n_Forest"),
  Features = c("Clinical Features\n + Promoter Activity", "Promoter\nActivity", "Clinical\nFeatures", 
               "Clinical Features\n + Promoter Activity", "Promoter\nActivity", "Clinical\nFeatures", 
               "Clinical Features\n + Promoter Activity", "Promoter\nActivity", "Clinical\nFeatures"),
  Mean_Training = c(mean_training_c_index_lasso_all, mean_training_c_index_lasso_prom, mean_training_c_index_lasso_clin,
                    mean_training_c_index_elastic_all, mean_training_c_index_elastic_prom,mean_training_c_index_elastic_clin,
                    mean_training_c_index_rf_all, mean_training_c_index_rf_prom, mean_training_c_index_rf_clin),
  SD_Training = c(sd_tr_lasso_all, sd_tr_lasso_prom, sd_tr_lasso_clin, 
                  sd_tr_elastic_all, sd_tr_elastic_prom,sd_tr_elastic_clin,
                  sd_tr_rf_all, sd_tr_rf_prom, sd_tr_rf_clin),
  Validation = c(test_lasso_all, test_lasso_prom, test_lasso_clin,
                 test_elastic_all, test_elastic_prom, test_elastic_clin,
                 test_rf_all, test_rf_prom, test_rf_clin),
  Validation_sd = c(test_lasso_all_sd, test_lasso_prom_sd, test_lasso_clin_sd,
                    test_elastic_all_sd, test_elastic_prom_sd, test_elastic_clin_sd,
                    test_rf_all_sd, test_rf_prom_sd, test_rf_clin_sd)
)

data$Features <- factor(data$Features, levels = c("Clinical\nFeatures","Promoter\nActivity", "Clinical Features\n + Promoter Activity" ))
data$Method <- as.factor(data$Method)

# Reshape the data to long format
data_long <- data %>%
  pivot_longer(cols = c("Mean_Training", "Validation"),
               names_to = "Metric",
               values_to = "Mean") %>%
  pivot_longer(cols = c("SD_Training", "Validation_sd"),
               names_to = "Metric_sd",
               values_to = "SD") %>%
  filter((Metric == "Mean_Training" & Metric_sd == "SD_Training") | 
           (Metric == "Validation" & Metric_sd == "Validation_sd")) %>%
  select(-Metric_sd)

# Plot:
p = ggplot(data_long, aes(x = Features, y = Mean, color = Metric, group = Metric)) +
  geom_line(size = 1.2) +  # Line for each metric (training and validation)
  geom_point(size = 3) +  # Points at each data point
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.1, color = "black", size = 0.2) +  # Error bars
  facet_grid(Method ~ ., scales = "free_y", switch = "y") + 
  labs(x = "Feature Set",
       y = "C-index (Mean ± SD)") +
  theme_grey(base_size = 14) +
  theme(
    axis.title = element_text(size = 16, face = "bold",colour = "black"),
    axis.text = element_text(size = 15,colour = "black"),
    axis.text.x = element_text(angle = 40, hjust = 1,face = "bold"), 
    strip.text = element_text(size = 14, face = "bold"),
    plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_blank()
  ) +
  scale_color_manual(values =c("Mean_Training" = "#298c8c", "Validation" = "#ea801c"),
                     labels = c("Training Mean", "Validation Mean")) +
  scale_x_discrete(expand = c(0.05,  0.05)) 
p
ggsave("c_index.png", plot = p, width = 4, height = 8, dpi = 1000)

######################################################################################
######################################################################################


## 2. For IBS Figure 5h:------------------------------------>>

# Read the files having scores:
all_features <- read_csv("clin_and_prom.csv")
prom_only <- read_csv("prom_only.csv")
clin_only <- read_csv("clin_only.csv")

# Omit-NA values
all_features <- na.omit(all_features)
prom_only <- na.omit(prom_only)
clin_only <- na.omit(clin_only)

# Calculate standard deviation of performance on test datasets at different seeds
# Random forest
test_rf_all_sd <- sd(all_features$ibs_test_rf)
test_rf_prom_sd <- sd(prom_only$ibs_test_rf)
test_rf_clin_sd <- sd(clin_only$ibs_test_rf)

# Lasso
test_lasso_all_sd <- sd(all_features$ibs_test_lasso)
test_lasso_prom_sd <- sd(prom_only$ibs_test_lasso)
test_lasso_clin_sd <- sd(clin_only$ibs_test_lasso)

# Elastic net
test_elastic_all_sd <- sd(all_features$ibs_test_elastic)
test_elastic_prom_sd <- sd(prom_only$ibs_test_elastic)
test_elastic_clin_sd <- sd(clin_only$ibs_test_elastic)

# Calculate mean performance on test dataset at different seeds
# Random forest
test_rf_all <- mean(all_features$ibs_test_rf)
test_rf_prom <- mean(prom_only$ibs_test_rf)
test_rf_clin <- mean(clin_only$ibs_test_rf)

# Lasso
test_lasso_all <- mean(all_features$ibs_test_lasso)
test_lasso_prom <- mean(prom_only$ibs_test_lasso)
test_lasso_clin <- mean(clin_only$ibs_test_lasso)

# Elastic net
test_elastic_all <- mean(all_features$ibs_test_elastic)
test_elastic_prom <- mean(prom_only$ibs_test_elastic)
test_elastic_clin <- mean(clin_only$ibs_test_elastic)

# Mean of standard deviation on train dataset at different seeds
# Random forest
sd_tr_rf_all <- mean(all_features$std_ibs_rf)
sd_tr_rf_prom <- mean(prom_only$std_ibs_rf)
sd_tr_rf_clin <- mean(clin_only$std_ibs_rf)

# Lasso
sd_tr_lasso_all <- mean(all_features$std_ibs_lasso)
sd_tr_lasso_prom <- mean(prom_only$std_ibs_lasso)
sd_tr_lasso_clin <- mean(clin_only$std_ibs_lasso)

# Elastic net
sd_tr_elastic_all <- mean(all_features$std_ibs_elastic)
sd_tr_elastic_prom <- mean(prom_only$std_ibs_elastic)
sd_tr_elastic_clin <- mean(clin_only$std_ibs_elastic)

# Calculate mean of mean performance on train dataset at different seeds
# Random forest
mean_training_ibs_rf_all <- mean(all_features$mean_ibs_rf)
mean_training_ibs_rf_prom <- mean(prom_only$mean_ibs_rf)
mean_training_ibs_rf_clin <- mean(clin_only$mean_ibs_rf)

# Elastic net
mean_training_ibs_elastic_all <- mean(all_features$mean_ibs_elastic)
mean_training_ibs_elastic_prom <- mean(prom_only$mean_ibs_elastic)
mean_training_ibs_elastic_clin <- mean(clin_only$mean_ibs_elastic)

# Lasso
mean_training_ibs_lasso_all <- mean(all_features$mean_ibs_lasso)
mean_training_ibs_lasso_prom <- mean(prom_only$mean_ibs_lasso)
mean_training_ibs_lasso_clin <- mean(clin_only$mean_ibs_lasso)

# Summarising all information in a dataframe
data <- data.frame(
  Method = c("Lasso", "Lasso", "Lasso", "Elastic", "Elastic", "Elastic", "Random_Survival\n_Forest", "Random_Survival\n_Forest", "Random_Survival\n_Forest"),
  Features = c("Clinical Features\n + Promoter Activity", "Promoter\nActivity", "Clinical\nFeatures", 
               "Clinical Features\n + Promoter Activity", "Promoter\nActivity", "Clinical\nFeatures", 
               "Clinical Features\n + Promoter Activity", "Promoter\nActivity", "Clinical\nFeatures"),
  Mean_Training = c(mean_training_ibs_lasso_all, mean_training_ibs_lasso_prom, mean_training_ibs_lasso_clin,
                    mean_training_ibs_elastic_all, mean_training_ibs_elastic_prom,mean_training_ibs_elastic_clin,
                    mean_training_ibs_rf_all, mean_training_ibs_rf_prom, mean_training_ibs_rf_clin),
  SD_Training = c(sd_tr_lasso_all, sd_tr_lasso_prom, sd_tr_lasso_clin, 
                  sd_tr_elastic_all, sd_tr_elastic_prom,sd_tr_elastic_clin,
                  sd_tr_rf_all, sd_tr_rf_prom, sd_tr_rf_clin),
  Validation = c(test_lasso_all, test_lasso_prom, test_lasso_clin,
                 test_elastic_all, test_elastic_prom, test_elastic_clin,
                 test_rf_all, test_rf_prom, test_rf_clin),
  
  Validation_sd = c(test_lasso_all_sd, test_lasso_prom_sd, test_lasso_clin_sd,
                    test_elastic_all_sd, test_elastic_prom_sd, test_elastic_clin_sd,
                    test_rf_all_sd, test_rf_prom_sd, test_rf_clin_sd)
  
)

data$Features <- factor(data$Features, levels = c("Clinical\nFeatures","Promoter\nActivity", "Clinical Features\n + Promoter Activity" ))
data$Method <- as.factor(data$Method)

# Reshape the data to long format
data_long <- data %>%
  pivot_longer(cols = c("Mean_Training", "Validation"),
               names_to = "Metric",
               values_to = "Mean") %>%
  pivot_longer(cols = c("SD_Training", "Validation_sd"),
               names_to = "Metric_sd",
               values_to = "SD") %>%
  filter((Metric == "Mean_Training" & Metric_sd == "SD_Training") | 
           (Metric == "Validation" & Metric_sd == "Validation_sd")) %>%
  select(-Metric_sd)

# Plot:
p = ggplot(data_long, aes(x = Features, y = Mean, color = Metric, group = Metric)) +
  geom_line(size = 1.2) +  
  geom_point(size = 3) +  
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.1, color = "black", size = 0.2) +  # for Error bars
  facet_grid(Method ~ ., scales = "free_y", switch = "y") + 
  labs(x = "Feature Set",
       y = "IBS (Mean ± SD)") +
  theme_grey(base_size = 14) +
  theme(
    axis.title = element_text(size = 16, face = "bold",colour = "black"),
    axis.text = element_text(size = 15,colour = "black"),
    axis.text.x = element_text(angle = 40, hjust = 1,face = "bold"),  
    strip.text = element_text(size = 14, face = "bold"),
    plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm")
  ) +
  scale_color_manual(values =c("Mean_Training" = "#298c8c", "Validation" = "#ea801c"),
                     labels = c("Training Mean", "Validation Mean")) +
  scale_x_discrete(expand = c(0.05,  0.05)) 
p
ggsave("IBS_score.png", plot = p, width = 4, height = 8, dpi = 1000)

#####################################################################################
#####################################################################################


## 3. For AUROC Figure 5i------------------------------------->>

# Read the files having scores:
all_features <- read_csv("clin_and_prom.csv")
prom_only <- read_csv("prom_only.csv")
clin_only <- read_csv("clin_only.csv")

# Omit-NA values
all_features <- na.omit(all_features)
prom_only <- na.omit(prom_only)
clin_only <- na.omit(clin_only)

# Calculate standard deviation of performance on test datasets at different seeds
# Random forest
test_rf_all_sd <- sd(all_features$auc_test_rf)
test_rf_prom_sd <- sd(prom_only$auc_test_rf)
test_rf_clin_sd <- sd(clin_only$auc_test_rf)

# Lasso
test_lasso_all_sd <- sd(all_features$auc_test_lasso)
test_lasso_prom_sd <- sd(prom_only$auc_test_lasso)
test_lasso_clin_sd <- sd(clin_only$auc_test_lasso)

# Elastic net
test_elastic_all_sd <- sd(all_features$auc_test_elastic)
test_elastic_prom_sd <- sd(prom_only$auc_test_elastic)
test_elastic_clin_sd <- sd(clin_only$auc_test_elastic)

# Calculate mean performance on test dataset at different seeds
# Random forest
test_rf_all <- mean(all_features$auc_test_rf)
test_rf_prom <- mean(prom_only$auc_test_rf)
test_rf_clin <- mean(clin_only$auc_test_rf)

# Lasso
test_lasso_all <- mean(all_features$auc_test_lasso)
test_lasso_prom <- mean(prom_only$auc_test_lasso)
test_lasso_clin <- mean(clin_only$auc_test_lasso)

# Elastic net
test_elastic_all <- mean(all_features$auc_test_elastic)
test_elastic_prom <- mean(prom_only$auc_test_elastic)
test_elastic_clin <- mean(clin_only$auc_test_elastic)


# Calculating mean of standard deviation on train dataset at different seeds
# Random forest
sd_tr_rf_all <- mean(all_features$std_auc_rf)
sd_tr_rf_prom <- mean(prom_only$std_auc_rf)
sd_tr_rf_clin <- mean(clin_only$std_auc_rf)

# Lasso
sd_tr_lasso_all <- mean(all_features$std_auc_lasso)
sd_tr_lasso_prom <- mean(prom_only$std_auc_lasso)
sd_tr_lasso_clin <- mean(clin_only$std_auc_lasso)

# Elastic net
sd_tr_elastic_all <- mean(all_features$std_auc_elastic)
sd_tr_elastic_prom <- mean(prom_only$std_auc_elastic)
sd_tr_elastic_clin <- mean(clin_only$std_auc_elastic)


# Calculating mean of mean performance on train dataset at different seeds
# Random forest
mean_training_auroc_rf_all <- mean(all_features$mean_auc_rf)
mean_training_auroc_rf_prom <- mean(prom_only$mean_auc_rf)
mean_training_auroc_rf_clin <- mean(clin_only$mean_auc_rf)

# Elastic net
mean_training_auroc_elastic_all <- mean(all_features$mean_auc_elastic)
mean_training_auroc_elastic_prom <- mean(prom_only$mean_auc_elastic)
mean_training_auroc_elastic_clin <- mean(clin_only$mean_auc_elastic)

# Lasso
mean_training_auroc_lasso_all <- mean(all_features$mean_auc_lasso)
mean_training_auroc_lasso_prom <- mean(prom_only$mean_auc_lasso)
mean_training_auroc_lasso_clin <- mean(clin_only$mean_auc_lasso)

# Summarising all information in a dataframe
data <- data.frame(
  Method = c("Lasso", "Lasso", "Lasso", "Elastic", "Elastic", "Elastic", "Random_Survival\n_Forest", "Random_Survival\n_Forest", "Random_Survival\n_Forest"),
  Features = c("Clinical Features\n + Promoter Activity", "Promoter\nActivity", "Clinical\nFeatures", 
               "Clinical Features\n + Promoter Activity", "Promoter\nActivity", "Clinical\nFeatures", 
               "Clinical Features\n + Promoter Activity", "Promoter\nActivity", "Clinical\nFeatures"),
  Mean_Training = c(mean_training_auroc_lasso_all, mean_training_auroc_lasso_prom, mean_training_auroc_lasso_clin,
                    mean_training_auroc_elastic_all, mean_training_auroc_elastic_prom,mean_training_auroc_elastic_clin,
                    mean_training_auroc_rf_all, mean_training_auroc_rf_prom, mean_training_auroc_rf_clin),
  SD_Training = c(sd_tr_lasso_all, sd_tr_lasso_prom, sd_tr_lasso_clin, 
                  sd_tr_elastic_all, sd_tr_elastic_prom,sd_tr_elastic_clin,
                  sd_tr_rf_all, sd_tr_rf_prom, sd_tr_rf_clin),
  Validation = c(test_lasso_all, test_lasso_prom, test_lasso_clin,
                 test_elastic_all, test_elastic_prom, test_elastic_clin,
                 test_rf_all, test_rf_prom, test_rf_clin),
  Validation_sd = c(test_lasso_all_sd, test_lasso_prom_sd, test_lasso_clin_sd,
                    test_elastic_all_sd, test_elastic_prom_sd, test_elastic_clin_sd,
                    test_rf_all_sd, test_rf_prom_sd, test_rf_clin_sd)
  
)

data$Features <- factor(data$Features, levels = c("Clinical\nFeatures","Promoter\nActivity", "Clinical Features\n + Promoter Activity" ))
data$Method <- as.factor(data$Method)

# Reshape the data to long format
data_long <- data %>%
  pivot_longer(cols = c("Mean_Training", "Validation"),
               names_to = "Metric",
               values_to = "Mean") %>%
  pivot_longer(cols = c("SD_Training", "Validation_sd"),
               names_to = "Metric_sd",
               values_to = "SD") %>%
  filter((Metric == "Mean_Training" & Metric_sd == "SD_Training") | 
           (Metric == "Validation" & Metric_sd == "Validation_sd")) %>%
  select(-Metric_sd)

# Plot:
p = ggplot(data_long, aes(x = Features, y = Mean, color = Metric, group = Metric)) +
  geom_line(size = 1.2) +  # Line for each metric (training and validation)
  geom_point(size = 3) +  # Points at each data point
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.1, color = "black", size = 0.2) +  # Error bars
  facet_grid(Method ~ ., scales = "free_y", switch = "y") + 
  labs(x = "Feature Set",
       y = "Time-dependent Area under the ROC (Mean ± SD)") +
  theme_grey(base_size = 14) +
  theme(
    axis.title = element_text(size = 16, face = "bold",colour = "black"),
    axis.text = element_text(size = 15,colour = "black"),
    axis.text.x = element_text(angle = 40, hjust = 1,face = "bold"), 
    plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"),
    strip.text = element_text(size = 14, face = "bold")
  ) +
  scale_color_manual(values =c("Mean_Training" = "#298c8c", "Validation" = "#ea801c"),
                     labels = c("Training Mean", "Validation Mean")) +
  scale_x_discrete(expand = c(0.05,  0.05))  
p
ggsave("auroc_score.png", plot = p, width = 4, height = 8, dpi = 1000)

######################################################################################
######################################################################################

## 1. For NRI Figure 5j ------------------------------------------>>

## FOR TRAINING SET RESULTS:-------------------------------------------------##
# Read the files having scores:
clin_vs_all_train <- as.data.frame(read_csv("NRI_train/NRI_CLIN_vs_allfeature_train.csv"))
# Omit-NA values
clin_vs_all_train <- na.omit(clin_vs_all_train)

## Subset for different methods
# Random forest
random_forest_train <- clin_vs_all_train[clin_vs_all_train$model=="random_forest",]

# Elastic net
elastic_net_train <- clin_vs_all_train[clin_vs_all_train$model=="ELASTIC_NET",]

# Lasso
lasso_train <- clin_vs_all_train[clin_vs_all_train$model=="LASSO",]

# Combine all training df
df_train <- rbind(random_forest_train, elastic_net_train, lasso_train)
df_train$Set <- "Training"

## FOR TESTING SET RESULTS:-------------------------------------------------##
# Read the files having scores:
clin_vs_all_test <- as.data.frame(read_csv("NRI_test/NRI_CLIN_vs_all_features.csv"))
# Omit-NA values
clin_vs_all_test <- na.omit(clin_vs_all_test)

## Subset for different methods
# Random forest
random_forest_test <- clin_vs_all_test[clin_vs_all_test$model=="random_forest",]

# Elastic net
elastic_net_test <- clin_vs_all_test[clin_vs_all_test$model=="ELASTIC_NET",]

# Lasso
lasso_test <- clin_vs_all_test[clin_vs_all_test$model=="LASSO",]

# Combine testing df
df_test <- rbind(random_forest_test, elastic_net_test, lasso_test)
df_test$Set <- "Testing"


## Combine both training and testing dataframes
df_combined <- rbind(df_train, df_test)

# Group & calculate the standard deviation
nri_summary <- df_combined %>%
  group_by(Set, model, percentile) %>%
  summarise(
    value = mean(nri_total, na.rm = TRUE),
    .groups = "drop"
  )

# Format labels
nri_summary$model <- gsub("ELASTIC_NET","Elastic net",nri_summary$model)
nri_summary$model <- gsub("LASSO","Lasso",nri_summary$model)
nri_summary$model <- gsub("random_forest","Random forest",nri_summary$model)

# Add 'th' to each percentile value
nri_summary$percentile <- paste0(nri_summary$percentile, "th")

# Set the order of percentiles for plotting
nri_summary$percentile <- factor(nri_summary$percentile, levels = c("25th", "50th", "75th"))

# Round to two digits
nri_summary$value <- round(nri_summary$value, 2)

# Plot
p <- ggplot(nri_summary, aes(
  x = percentile,
  y = value,
  group = interaction(model, Set),
  color = model,
  linetype = Set
)) +
  geom_point(size = 3) +
  geom_line(size = 1.2) +  
  labs(
    y = "NRI",
    x = "Percentile",
    color = "Classifier"
  ) +
  scale_y_continuous(
    limits = c(0.3, 0.7),
    breaks = seq(0.3, 0.7, by = 0.02),  # More frequent breaks
    expand = c(0, 0)
  ) +
  scale_color_manual(values = c(
    "Elastic net" = "#1f77b4",
    "Lasso" = "#ff7f0e",
    "Random forest" = "#2ca02c"
  )) +
  theme_bw() +
  # theme_grey(base_size = 14) +
  theme(
    axis.text.x = element_text(colour = "black", face = "bold", size = 12),
    axis.text.y = element_text(colour = "black", size = 14),
    plot.margin = margin(t = 0.2, r = 0.005, b = 0.05, l = 0.55, unit = "cm"),
    axis.title = element_text(colour = "black", face = "bold", size = 16),
    text = element_text(colour = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
  )
p
ggsave("C:/Users/Simran/Desktop/Promoter activity project/PAPER/NPJ_Revision1/NRI_FINAL_updated.png", plot = p, width = 10, height = 6, dpi = 1000)
###################################################################################################################################################
###################################################################################################################################################
