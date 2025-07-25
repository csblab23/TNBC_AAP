## Code to make a box plot for the 5000 iteration cutoff values for obtained 5 AAPs and their genes ##

# Load libraries
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyverse)

## 1. For promoters
# Load the dataframe having all cutoff values
prom_cutoff <- as.data.frame(read_csv("promoters_cutoff_median_working_with.csv"))

# Keep just the selected 5 promoters
keep <- c('pr3004_ENSG00000086758.17', 'pr9001_ENSG00000122786.20', 
          'pr12191_ENSG00000136560.14', 'pr19936_ENSG00000167460.17',
          'pr41492_ENSG00000230590.13')

prom_cutoff_sel <- prom_cutoff[prom_cutoff$promoter %in% keep,]

# Adding gene names corresponding to ensembl ids
genenames_list <- c('pr3004_HUWE1', 'pr9001_CALD1', 'pr12191_TANK', 'pr19936_TPM4', 'pr41492_FTX' )
prom_cutoff_sel$pr_genes <- genenames_list 
prom_cutoff_sel <- prom_cutoff_sel %>% relocate(pr_genes)

# Split the column 3 having comma-separated 5000 iteration cutoff values
df_expanded <- prom_cutoff_sel %>%
  rowwise() %>%
  mutate(cutoff_vec = strsplit(all_cutoffs, ",")) %>%
  unnest_longer(cutoff_vec) %>%
  mutate(cutoff_vec = as.numeric(cutoff_vec)) %>%
  ungroup()
df_expanded <- as.data.frame(df_expanded)
df_expanded$pr_genes <- factor(df_expanded$pr_genes,levels = c("pr3004_HUWE1","pr41492_FTX",
                                                               "pr9001_CALD1","pr12191_TANK","pr19936_TPM4"))

# Make boxplot
g <- ggplot(df_expanded, aes(x = pr_genes, y = cutoff_vec)) +
  geom_boxplot(outlier.size = 0.5,size = 0.6) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(size = 22, colour = "black",angle = 45, hjust = 1),
    axis.text.y = element_text(size = 22, colour = "black"),
    plot.title = element_text(size = 25, hjust = 0.5,face = "bold"),  
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    legend.text = element_blank(),
    legend.title = element_blank()
  ) +  labs(title = "Cutoff Distributions per AAP",
       x = "AAP_Gene",
       y = "5000 Iterations Cutoff Values")
g
ggsave("Promoter_5000iterations_boxplot.png", plot = g, width = 8, height = 7, units = "in",
       bg = "white")

#############################################################################################
#############################################################################################

## 2. For the gene:
# Load the dataframe having all cutoff values
gene_cutoff <- as.data.frame(read_csv("gene_cutoff_median_working_with.csv"))
colnames(gene_cutoff)[1]<-'gene'

# Keep just the selected 5 promoters
keep <- c('ENSG00000086758.17', 'ENSG00000122786.20', 'ENSG00000136560.14', 'ENSG00000167460.17', 'ENSG00000230590.13')
gene_cutoff_sel <- gene_cutoff[gene_cutoff$gene %in% keep,]

# Adding gene names corresponding to ensembl ids
genenames_list <- c('HUWE1', 'CALD1', 'TANK', 'TPM4', 'FTX' )
gene_cutoff_sel$genenames <- genenames_list 
gene_cutoff_sel <- gene_cutoff_sel %>% relocate(genenames)

# Split the column 3 having comma-separated 5000 iteration cutoff values
df_expanded <- gene_cutoff_sel %>%
  rowwise() %>%
  mutate(cutoff_vec = strsplit(all_cutoffs, ",")) %>%
  unnest_longer(cutoff_vec) %>%
  mutate(cutoff_vec = as.numeric(cutoff_vec)) %>%
  ungroup()
df_expanded <- as.data.frame(df_expanded)
df_expanded$genenames <- factor(df_expanded$genenames,levels = c("HUWE1","FTX",
                                                               "CALD1","TANK","TPM4"))
# Make boxplot
g <- ggplot(df_expanded, aes(x = genenames, y = cutoff_vec)) +
  geom_boxplot(outlier.size = 0.5,size = 0.6) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(size = 22, colour = "black",angle = 45, hjust = 1),
    axis.text.y = element_text(size = 22, colour = "black"),
    plot.title = element_text(size = 25, hjust = 0.5,face = "bold"),  
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
    legend.text = element_blank(),
    legend.title = element_blank()
  ) +  labs(title = "Cutoff Distributions per Gene",
            x = "Genes",
            y = "5000 Iterations Cutoff Values")
g
ggsave("Genes_5000iterations_boxplot.png", plot = g, width = 8, height = 7, units = "in",
       bg = "white")

#############################################################################################
#############################################################################################