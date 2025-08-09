## This code is to make the box plots presented in Supplementary S7a for the Feature set C,A,F comparison ##

# Load required libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(EnvStats)

# For C-index:--------------------------------->>

# Read the files having scores:
# Define function to load, clean, and label a dataset
load_model_data <- function(path, model_name) {
  read_csv(path) %>%
    select(c_index_elastic_test, c_index_rf_test, c_index_lasso_test) %>%  # Only keep relevant columns
    na.omit() %>%
    mutate(Model = model_name)
}

# Load each model's data
fsc <- load_model_data("5_feature_model/clin_and_prom.csv", "FS_C")
fsa <- load_model_data("4_feature_model/clin_and_prom.csv", "FS_A")
fsf <- load_model_data("6_feature_updated/clin_and_prom.csv", "FS_F")

# Combine all into one dataframe
combined <- bind_rows(fsc, fsa, fsf)

# Reshape from wide to long format
long_data <- combined %>%
  pivot_longer(cols = starts_with("c_index"),
               names_to = "Method",
               values_to = "CIndex") %>%
  mutate(
    Method = case_when(
      Method == "c_index_elastic_test" ~ "Elastic net",
      Method == "c_index_rf_test" ~ "Random forest",
      Method == "c_index_lasso_test" ~ "Lasso"
    )
  )

# Define comparisons between models 
model_comparisons <- list(
  c("FS_C", "FS_A"),
  c("FS_C", "FS_F"),
  c("FS_A", "FS_F")
)
# Converting to factor to define order
long_data$Model <- factor(long_data$Model,levels = c("FS_C","FS_A","FS_F"))

# Boxplot
p <- ggplot(long_data, aes(x = Model, y = CIndex, fill = Model)) +
  geom_boxplot() +
  stat_compare_means(
    comparisons = model_comparisons,
    method = "wilcox.test",
    label = "p.format"
  ) +
  facet_wrap(~ Method) +
  labs(title = "C-Index Comparison",
       x = "Feature Set",
       y = "C-Index") +
  theme_bw() +
  scale_fill_brewer(palette = "Pastel1") +
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.title = element_text(face = "bold",size = 15),
        plot.title = element_text(colour = "black",size = 15,hjust = 0.5,face = "bold"),
        text = element_text(colour = "black"),
        strip.text = element_text(colour = "black",face = "bold",size=14),
        legend.position = "none")

p
ggsave("c_index_compare.png", plot = p, width = 8, height = 5, dpi = 1000)

##########################################################################

# For IBS:--------------------------------->>

# Define function to load, clean, and label a dataset
load_model_data <- function(path, model_name) {
  read_csv(path) %>%
    select(ibs_test_elastic, ibs_test_rf, ibs_test_lasso) %>%  # Only keep relevant columns
    na.omit() %>%
    mutate(Model = model_name)
}

# Load each model's data
fsc <- load_model_data("5_feature_model/clin_and_prom.csv", "FS_C")
fsa <- load_model_data("4_feature_model/clin_and_prom.csv", "FS_A")
fsf <- load_model_data("6_feature_updated/clin_and_prom.csv", "FS_F")

# Combine all into one dataframe
combined <- bind_rows(fsc, fsa, fsf)

# Reshape from wide to long format
long_data <- combined %>%
  pivot_longer(cols = starts_with("ibs"),
               names_to = "Method",
               values_to = "IBS") %>%
  mutate(
    Method = case_when(
      Method == "ibs_test_elastic" ~ "Elastic net",
      Method == "ibs_test_rf" ~ "Random forest",
      Method == "ibs_test_lasso" ~ "Lasso"
    )
  )

# Define comparisons between models
model_comparisons <- list(
  c("FS_C", "FS_A"),
  c("FS_C", "FS_F"),
  c("FS_A", "FS_F")
)
# Converting to factor to define order
long_data$Model <- factor(long_data$Model,levels = c("FS_C","FS_A","FS_F"))

# Boxplot
p <- ggplot(long_data, aes(x = Model, y = IBS, fill = Model)) +
  geom_boxplot() +
  stat_compare_means(
    comparisons = model_comparisons,
    method = "wilcox.test",
    label = "p.format"
  ) +
  facet_wrap(~ Method) +
  labs(title = "IBS Comparison",
       x = "Feature Set",
       y = "IBS") +
  theme_bw() +
  scale_fill_brewer(palette = "Pastel1") +
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.title = element_text(face = "bold",size = 15),
        plot.title = element_text(colour = "black",size = 15,hjust = 0.5,face = "bold"),
        text = element_text(colour = "black"),
        strip.text = element_text(colour = "black",face = "bold",size=14),
        legend.position = "none")

p
ggsave("ibs_compare.png", plot = p, width = 8, height = 5, dpi = 1000)

##########################################################################

# For AUROC:--------------------------------->>

# Define function to load, clean, and label a dataset
load_model_data <- function(path, model_name) {
  read_csv(path) %>%
    select(auc_test_elastic, auc_test_rf, auc_test_lasso) %>%  # Only keep relevant columns
    na.omit() %>%
    mutate(Model = model_name)
}

# Load each model's data
fsc <- load_model_data("5_feature_model/clin_and_prom.csv", "FS_C")
fsa <- load_model_data("4_feature_model/clin_and_prom.csv", "FS_A")
fsf <- load_model_data("6_feature_updated/clin_and_prom.csv", "FS_F")

# Combine all into one dataframe
combined <- bind_rows(fsc, fsa, fsf)

# Reshape from wide to long format
long_data <- combined %>%
  pivot_longer(cols = starts_with("auc"),
               names_to = "Method",
               values_to = "AUROC") %>%
  mutate(
    Method = case_when(
      Method == "auc_test_elastic" ~ "Elastic net",
      Method == "auc_test_rf" ~ "Random forest",
      Method == "auc_test_lasso" ~ "Lasso"
    )
  )

# Define comparisons between models
model_comparisons <- list(
  c("FS_C", "FS_A"),
  c("FS_C", "FS_F"),
  c("FS_A", "FS_F")
)
# Converting to factor to define order
long_data$Model <- factor(long_data$Model,levels = c("FS_C","FS_A","FS_F"))

# Boxplot
p <- ggplot(long_data, aes(x = Model, y = AUROC, fill = Model)) +
  geom_boxplot() +
  stat_compare_means(
    comparisons = model_comparisons,
    method = "wilcox.test",
    label = "p.format"
  ) +
  facet_wrap(~ Method) +
  labs(title = "AUROC Comparison",
       x = "Feature Set",
       y = "AUROC") +
  theme_bw() +
  scale_fill_brewer(palette = "Pastel1") +
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.title = element_text(face = "bold",size = 15),
        plot.title = element_text(colour = "black",size = 15,hjust = 0.5,face = "bold"),
        text = element_text(colour = "black"),
        strip.text = element_text(colour = "black",face = "bold",size=14),
        legend.position = "none")

p
ggsave("auroc_compare.png", plot = p, width = 8, height = 5, dpi = 1000)

#########################################################################################
#########################################################################################
