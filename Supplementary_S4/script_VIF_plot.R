#### NOTE: This script contains the code used for calculating Mean Variance Inflation Factor (VIF) per feature & generate the Supplementary figure S4a ####

# Load required libraries
library(dplyr)
library(tidyverse)
library(ggplot2)

# Load the calculated VIF values dataframe
vif <- as.data.frame(read_csv('100_states_vif.csv'))

# Calculate the average VIF per feature
avg_vif <- aggregate(VIF ~ Feature, data = vif, FUN = mean)
sd_vif <- aggregate(VIF ~ Feature, data = vif, FUN = sd)

# Merge them together
vif_summary <- merge(avg_vif, sd_vif, by = "Feature")
colnames(vif_summary) <- c("Feature", "Mean_VIF", "SD_VIF")
vif_summary <- vif_summary[vif_summary$Feature != 'const',]  # Remove rows corresponding to constant

# Sort Feature factor levels by ascending Mean_VIF
vif_summary$Feature <- factor(vif_summary$Feature, levels = vif_summary$Feature[order(vif_summary$Mean_VIF)])

# Plot
p <- ggplot(vif_summary, aes(x = Feature, y = Mean_VIF, group = 1)) +
  geom_line(color = "#2C3E50", size = 1) +
  geom_point(size = 3, color = "#1F77B4") +
  geom_errorbar(aes(ymin = Mean_VIF - SD_VIF, ymax = Mean_VIF + SD_VIF), 
                width = 0.2, color = "black", alpha = 0.8) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Mean Variance Inflation Factor (VIF) per Feature",
    x = "Feature",
    y = "Mean VIF"
  )
p
ggsave("VIF_plot.png", plot = p, width = 8, height = 5, dpi = 1000)
##################################################################################
##################################################################################
