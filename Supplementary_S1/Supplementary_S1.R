############## This script includes all the codes used for SUPPLEMENTARY S1 ################

# Load required libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(webr)
library(dplyr)

## For Figure S1a:--------------->>

# Read the result2 file for all 448 FUSCC samples obtained after proActiv
result2 <- readRDS("result2_tum_vs_adjnormal.rds")
result2_promoter=rowData(result2)
result2_promoter=as.data.frame(result2_promoter)

# Not considering internal promoters for downstream analysis
result_prom <- result2_promoter[result2_promoter$internalPromoter == 'FALSE',]

# Check for single-promoter and multi-promoter genes
data_plot <- as.data.frame(table(unname(table(result_prom$geneId))))
colnames(data_plot) <- c('No_of_promoters', "Frequency")
data_plot$No_of_promoters <- as.numeric(data_plot$No_of_promoters) 
data_plot_modified2 <- data_plot %>%
  mutate(No_of_promoters = ifelse(No_of_promoters >= 4, 4, No_of_promoters))

# Obtain the total number of genes
total_genes <- sum(data_plot_modified2$Frequency)
data_plot_modified2$percent_of_genes <- (data_plot_modified2$Frequency / total_genes) * 100

# Compile multi-promoter genes
data <- data_plot_modified2[1:4,]
sum(data_plot_modified2$percent_of_genes[4:15]) # 1.754201
data$percent_of_genes[4] <- sum(data_plot_modified2$percent_of_genes[4:15])
data$Frequency[4] <- sum(data_plot_modified2$Frequency[4:15]) 
data$No_of_promoters <- as.character(data$No_of_promoters)
data$No_of_promoters <- gsub("4",">=4",data$No_of_promoters)
data$No_of_promoters <- factor(data$No_of_promoters, levels = c('1', '2', '3', '>=4'))
data$final_freq <- round(data$percent_of_genes, 1)

# Calculate the difference to make the sum 100
diff <- 100 - sum(data$final_freq)
# Add the difference to the smallest positive value
data$final_freq[which.min(data$final_freq)] <- data$final_freq[which.min(data$final_freq)] + diff

# Define a color palette for the levels
color_palette <- c("1" = "peachpuff2", 
                   "2" = "#A0522D", 
                   "3" = "#CD853F", 
                   ">=4" = "#E3B76D")

# Plot:
g<-ggplot(data, aes(x = No_of_promoters, y = final_freq,fill=No_of_promoters)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "",
       x = "Number of Promoters",
       y = "% of Genes") +
  scale_x_discrete(labels = c('1', '2', '3', '>=4')) +
  scale_y_continuous(limits = c(0, 100),expand = c(0, 0)) +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 32, face = "bold", color = 'black'),
    axis.title.y = element_text(size = 32, face = "bold", color = 'black'),
    axis.text = element_text(size = 30, color = 'black'),
    axis.text.x = element_text(hjust = 0.5, size = 30),
    legend.title = element_blank(),
    legend.text = element_blank(),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "gray80"),
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', color = "gray90"),
    axis.line = element_line(size = 0.6, color = "black")
  ) +
  geom_text(aes(label = sprintf("%.2f%%", final_freq )),
            vjust = -0.5, size = 10, color = "black",fontface="bold")
g
ggsave("suppl_fig1_barplot.png", g, width = 10, height = 10, dpi = 1000)

#########################################################################
#########################################################################

## For Supplementary S1b & c:------------------>>

# Making position-wise plots separately for major and minor promoters:
# First taking only major promoters
major <- pdata2_summary[pdata2_summary$Tumor.class=='Major',]

# Keep position 'first position' and 'downstream' (only these two categories):
major$promoterPosition2 <- ifelse(major$promoterPosition==1,'First_position','Downstream')
major <- major %>%
  group_by(promoterPosition2) %>%
  mutate(prop = sum(count) / sum(major$count) * 100) %>%
  ungroup()

# Create a dataframe
df_major <- data.frame(promoterPosition2=c('First_position','Downstream'),
                       prop = c(56.9,43.1),
                       Tumor.class = 'Major')
# Convert to factor
df_major$promoterPosition2 <- factor(df_major$promoterPosition2,levels = c('First_position','Downstream'))

# Plot
gp <- ggplot(df_major) +
  geom_col(aes(x = promoterPosition2, y = prop, fill = Tumor.class)) +
  geom_text(aes(x = promoterPosition2, y = prop, label = paste0(prop, "%")), 
            position = position_stack(vjust = 0.5), size = 15, fontface="bold" ,family = "Arial", color = "black") +
  ylab('Proportion of promoter types') + 
  labs(x = expression(bold(Promoter ~ Position ~ "(5'" %->% " 3')"))) +
  scale_fill_manual(values = c("Major" = "cadetblue")) +
  scale_y_continuous(limits = c(0, 80), expand = c(0, 0)) +  # Set limits and remove expansion
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 28, family = "Arial", colour = 'black'),
    axis.title.x = element_text(size = 30, face = 'bold', colour = "black"),
    axis.title = element_text(size = 30, face = "bold", family = "Arial", color = 'black'),
    legend.text = element_text(size = 30,colour = "black"),
    axis.text.x = element_text(size = 24, hjust = 0.5, colour = "black", family = "Arial", margin = margin(t = 5)),
    axis.title.y = element_text(size = 30, face = 'bold', colour = "black"),
    legend.title = element_text(size = 30, face = 'bold', colour = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  guides(fill = guide_legend(title = "Promoter\nCategory"))
gp
ggsave("suppl_fig1_major_position.png", gp, width = 9, height = 10, dpi = 1000)

###########-----------

# Now taking minor promoters
minor <- pdata2_summary[pdata2_summary$Tumor.class=='Minor',]

# Keep position 'first position' and 'downstream' (only these two categories):
minor$promoterPosition2 <- ifelse(minor$promoterPosition==1,'First_position','Downstream')
minor <- minor %>%
  group_by(promoterPosition2) %>%
  mutate(prop = sum(count) / sum(minor$count) * 100) %>%
  ungroup()

# Create a dataframe:
df_minor <- data.frame(promoterPosition2=c('First_position','Downstream'),
                       prop = c(36.0,64.0),
                       Tumor.class = 'Minor')
# Convert to factor
df_minor$promoterPosition2 <- factor(df_minor$promoterPosition2,levels = c('First_position','Downstream'))

# Plot:
gp <- ggplot(df_minor) +
  geom_col(aes(x = promoterPosition2, y = prop, fill = Tumor.class)) +
  geom_text(aes(x = promoterPosition2, y = prop, label = paste0(prop, "%")), 
            position = position_stack(vjust = 0.5), size = 15, fontface="bold" ,family = "Arial", color = "black") +
  ylab('Proportion of promoter types') + 
  labs(x = expression(bold(Promoter ~ Position ~ "(5'" %->% " 3')"))) +
  scale_fill_manual(values = c("Minor" = "#CCC591")) +
  scale_y_continuous(limits = c(0, 80), expand = c(0, 0)) +  
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 28, family = "Arial", colour = 'black'),
    axis.title.x = element_text(size = 30, face = 'bold', colour = "black"),
    axis.title = element_text(size = 30, face = "bold", family = "Arial", color = 'black'),
    legend.text = element_text(size = 30,colour = "black"),
    axis.text.x = element_text(size = 24, hjust = 0.5, colour = "black", family = "Arial", margin = margin(t = 5)),
    axis.title.y = element_text(size = 30, face = 'bold', colour = "black"),
    legend.title = element_text(size = 30, face = 'bold', colour = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  guides(fill = guide_legend(title = "Promoter\nCategory"))

gp
ggsave("suppl_fig1_minor_position.png", gp, width = 9, height = 10, dpi = 1000)

##################################################################################
##################################################################################
