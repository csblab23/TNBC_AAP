## This script contains the code to generate Supplementary figure S2 ##

# Load required packages
library(ggplot2)
library(dplyr)
library(tidyverse)

# NOTE: BED files were used to identify shared (via bedtools intersect) and unique (via bedtools subtract) regions between pr1077- and pr1079-associated transcript
# Load the dataframe based on the above analysis
df <- data.frame(
  Region = rep(c("5'UTR", "CDS", "3'UTR"), each = 3),
  Category = rep(c("Shared", "Unique_normal", "Unique_cancer"), times = 3),
  Value = c(0, 286, 245,   # 5'UTR
            1548, 174, 141, # CDS
            373, 1152, 0)   # 3'UTR
)

# Define total lengths of each functional region for both transcripts
total_normal <- c("5'UTR" = 286, "CDS" = 1722, "3'UTR" = 1525)
total_cancer <- c("5'UTR" = 245, "CDS" = 1689, "3'UTR" = 373)

# Function to get appropriate total length
get_total <- function(region, category) {
  if (category %in% c("Shared", "Unique_normal")) {
    return(total_normal[region])
  } else if (category == "Unique_cancer") {
    return(total_cancer[region])
  } else {
    return(NA)
  }
}

# Apply to get total and percentage
df$Total <- mapply(get_total, df$Region, df$Category)
df$Percentage <- round(100 * df$Value / df$Total, 1)

# Define order
df$Region <- factor(df$Region, levels = c("5'UTR", "CDS","3'UTR"))

# Plot
p <- ggplot(df, aes(x = Category, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  facet_wrap(~ Region, nrow = 1) +
  scale_fill_manual(values = c("Shared" = "gray80", 
                               "Unique_normal" = "#2ca02c", 
                               "Unique_cancer" = "#ff7f0e")) +
  labs(y = "% change in functional sequence") +
  theme(legend.position = "none", axis.text.x = element_text(angle=30)) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour = "black",size = 12.5),
    axis.title.y = element_text(colour = "black",size = 13,face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_text(colour = "black",size = 13,face = "bold"),
    legend.text = element_text(colour = "black", size = 12),
    strip.text = element_text(face = "bold", colour = "black", size = 13)
  )
p
ggsave("stackedbar_sequence.png", plot = p, width = 6, height = 3.4, dpi = 1000)

################################################################################
################################################################################
