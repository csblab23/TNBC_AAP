## Code to make the KM PLOT of Supplementary S3c ##

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

# Read input files 
result2 <- readRDS('result2_v2_updated.rds')
metadata <- as.data.frame(read_csv("meta_with_batchinfo.csv"))
count <- as.data.frame(read_csv("Featurecounts_combined_fuscc.csv"))

# Data preprocessing
metadata[1:5,1:5]
metadata <- metadata[, -1]  # Remove extra column
metadata <- metadata %>%
  mutate(Tumor_Normal1 = case_when(
    Tumor_Normal %in% c("TumorTissue", "TumorTissue.2018", "TumorTissue.rep") ~ "Tumor",
    Tumor_Normal %in% c("PairedNormalTissue", "PairedNormalTissue.rep") ~ "AdjNormal"
  ))
rownames(metadata) <- metadata$Run

# Take tumor samples
tumor_metadata <- subset(metadata, Tumor_Normal1 == 'Tumor')

# Not taking internal promoters' genes for downstream analysis
result2_promoter = rowData(result2) 
filtered_df <- as.data.frame(result2_promoter)
filtered_df <- filtered_df %>% 
  filter(internalPromoter == 'FALSE')

# Taking genes with >=2 promoters 
filtered_df <- filtered_df %>%
  group_by(geneId) %>%
  filter(n() >= 2) %>%
  mutate(ID = paste('pr', promoterId, '_', geneId, sep=''))
filtered_df <- filtered_df[, c('promoterId','geneId', 'ID')]   # subset required columns

# Load the promoter counts
promoter_counts <- as.data.frame(result2@assays@data@listData$promoterCounts)
keep = intersect(filtered_df$promoterId, rownames(promoter_counts))
promoter_counts <- promoter_counts[keep,]
colnames(promoter_counts) <- gsub("_SJ.out", "", colnames(promoter_counts))

# Merge the Dataframes
p_counts <- merge( filtered_df, promoter_counts, by.x = "promoterId", by.y = 'row.names')

gene_lists <- list()

for (sample_col in colnames(p_counts)[4:ncol(p_counts)]) {
  cat("Processing sample:", sample_col, "\n")
  filtered_genes <- c()  # Initialize an empty vector to store filtered genes for this sample
  # Loop through each unique gene in the dataset
  for (gene in unique(p_counts$geneId)) {
    
    # Extract expression counts for the current gene in the current sample
    gene_counts <- p_counts[p_counts$geneId == gene, sample_col]
    
    # Count how many promoters of that gene have expression values greater than 10
    num_ids_gt_threshold <- sum(gene_counts > 10)
    
    # If the gene has at least 2 promoters with expression > 10, retain it
    if (num_ids_gt_threshold >= 2) {
      filtered_genes <- c(filtered_genes, gene)
    }
  }
  # Store the filtered gene list for the current sample in a named list
  gene_lists[[sample_col]] <- filtered_genes
}

## For taking genes that are present in atleast 60% of samples ##---------------------------LEVEL-03
n_samples <- length(gene_lists)
ensg_freq <- table(unlist(gene_lists))
threshold <- 0.6 * n_samples
ensg_60percent <- names(ensg_freq[ensg_freq >= threshold])
filtered_df <- filtered_df[filtered_df$geneId %in% ensg_60percent,]

# Load the relative promoter activity
Rl_counts <- as.data.frame(result2@assays@data@listData$relativePromoterActivity) ##----------LEVEL-04
Rl_counts <- Rl_counts[complete.cases(Rl_counts),]
keep = intersect(filtered_df$promoterId, rownames(Rl_counts))
Rl_counts <- Rl_counts[keep,]
colnames(Rl_counts) <- gsub("_SJ.out", "", colnames(Rl_counts))

## For taking genes that have promoters with relative promoter activity > 0.25 in 10% of samples
tobetaken_rel <- rowSums(Rl_counts > 0.25) >= dim(Rl_counts)[2]*0.1
Rl_counts=Rl_counts[tobetaken_rel,]
filtered_df <- filtered_df[filtered_df$promoterId %in% rownames(Rl_counts),]
ensg_freq_2 <- table(filtered_df$geneId)
ensg_final <- names(ensg_freq_2[ensg_freq_2 >= 2])
filtered_df <- filtered_df[filtered_df$geneId %in% ensg_final,]

# Now we will pre-process feature count matrix here only, why here?? -> because we have to take genes that are present in the filtered step till here
count <- read_csv("Featurecounts_combined_fuscc.csv") %>%
  as.data.frame() %>%
  column_to_rownames("Geneid") %>% 
  subset(., select = -1) %>% 
  filter(rownames(.) %in% filtered_df$geneId)
cpm <- cpm(count)
is.exprs <- rowSums(cpm>1) >= dim(tumor_metadata)[1]*0.1
count <- count[is.exprs, ]
x <- DGEList(counts=count )
x <- calcNormFactors(x,method = "TMM")
v <- voom(x, plot=F)
vMat <- v$E

# Batch adjustment
vMat_2 <- sva::ComBat(vMat, batch=metadata$batch)
vMat_2 <- vMat_2[, tumor_metadata$Run] # removing normal samples
p_info <- filtered_df[filtered_df$geneId %in% row.names(vMat_2),]

# Taking Absolute Promoter Activity
Ab_counts <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)
Ab_counts <- Ab_counts[rownames(Ab_counts) %in% p_info$promoterId,]
colnames(Ab_counts) <- gsub("_SJ.out", "", colnames(Ab_counts))
Ab_counts <- merge( p_info, Ab_counts, by.x= 'promoterId', by.y = 'row.names')
rownames(Ab_counts) <- Ab_counts$ID
Ab_counts <- Ab_counts[, c(-1, -2, -3)]

# TAKING SAMPLES THAT HAVE RFS_time_months > 0 
surv_metadata =tumor_metadata[tumor_metadata$RFS_time_Months > 0,]
surv_metadata=surv_metadata[,c("Project_ID", "Run", "batch", "Tumor_Normal", "Tumor_Normal1", "RFS_Status", "RFS_time_Days", "RFS_time_Months")]
Ab_counts <- Ab_counts[, surv_metadata$Run]
vMat_2 <- vMat_2[, surv_metadata$Run]
tobetaken <- rowSums(Ab_counts > 0) >= dim(tumor_metadata)[1]*0.1
Ab_counts=Ab_counts[tobetaken,]

#######

# Function for performing matrix processing
perform_analysis_for_matrix <- function(metadata, matrix) {
  mm <- merge(metadata, matrix, by = 0)
  rownames(mm) <- mm$Row.names
  mm <- mm[, -1]
}

prom_matrix <- perform_analysis_for_matrix(surv_metadata, t(Ab_counts))

########----------------->>

# Loading the cut-off values & stratifying samples into high and low groups
load('promoters_cut_off_sub_int_false.RData')

mm_check_binary <- prom_matrix[, 9:dim(prom_matrix)[2]]
for (i in 1:ncol(mm_check_binary)) {
  mm_check_binary[, i] <- ifelse(prom_matrix[, i + 8] > final_cutoff_values[[i]], "high", "low")
}
mm_meta <- prom_matrix[, 1:8]
mm_finaldf <- merge(mm_meta, mm_check_binary, by = 0)
mm_finaldf <- tibble::column_to_rownames(mm_finaldf, "Row.names")
for (i in 9:dim(mm_finaldf)[2]) {
  if (is.character(mm_finaldf[, i])) {
    mm_finaldf[, i] <- factor(mm_finaldf[, i])
  }
}

## Selecting the gene of interest -- 'FTX'-> it has 3 promoters out of which 2 are non-significant and 1 is significant:

# For pr41491_FTX
grep('ENSG00000230590.13', colnames(mm_finaldf))
colnames(mm_finaldf)[123]
sfit <- survfit(Surv(mm_finaldf$RFS_time_Days, mm_finaldf$RFS_Status)~ mm_finaldf[,123], data=mm_finaldf)
test <- ggsurvplot(
  sfit,
  pval = TRUE,
  pval.method = TRUE,
  size = 1.5,
  palette = c("darkorchid3", "#298c8c"),
  title = "pr41491_FTX",
  xlab = "Time (Days)",
  ylab = "Relapse Free Survival",
  legend.labs = c("High (n = 220)", "Low (n = 138)"),
  legend.title = "Risk Groups",
  legend = c(0.75, 0.25),  # Adjust legend position inside the plot area
  risk.table = TRUE,
  risk.table.height = 0.2,
  surv.median.line = "hv",
  risk.table.fontsize = 6,
  fontsize = 16,
  pval.size = 11,
  font.main = 30,
  font.tickslab = 20,
  color.tickslab = 'black',
  size = 1,
  tables.y.text = FALSE,
  font.legend = 22,
  ggtheme = theme_minimal() +
    theme(
      plot.title = element_text(size = 18, face = 'bold', family = 'Arial', hjust = 0.5),
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
ggsave("20_12_pr_41491_FTX.png", test$plot, width = 7, height = 6, dpi = 1200,bg = "white")

############################################################################################
############################################################################################
