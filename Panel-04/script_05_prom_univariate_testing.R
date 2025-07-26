## Code to use the cutoff value calculated from the training dataset on the testing dataset by performing univariate survival analysis ##

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
metadata <- metadata[, -1]  # Remove extra column
metadata <- metadata %>%
  mutate(Tumor_Normal1 = case_when(
    Tumor_Normal %in% c("TumorTissue", "TumorTissue.2018", "TumorTissue.rep") ~ "Tumor",
    Tumor_Normal %in% c("PairedNormalTissue", "PairedNormalTissue.rep") ~ "AdjNormal"
  ))
rownames(metadata) <- metadata$Run
# Taking tumor samples
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
# Subset required columns
filtered_df <- filtered_df[, c('promoterId','geneId', 'ID')]

# Load promoter counts
promoter_counts <- as.data.frame(result2@assays@data@listData$promoterCounts)
keep = intersect(filtered_df$promoterId, rownames(promoter_counts))
promoter_counts <- promoter_counts[keep,]
colnames(promoter_counts) <- gsub("_SJ.out", "", colnames(promoter_counts))
# Merge the Dataframes
p_counts <- merge( filtered_df, promoter_counts, by.x = "promoterId", by.y = 'row.names')

gene_lists <- list()

# Loop through each sample column in the dataset
for (sample_col in colnames(p_counts)[4:ncol(p_counts)]) {
  cat("Processing sample:", sample_col, "\n")
  # Initialize an empty vector to store filtered genes for this sample
  filtered_genes <- c()
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

# Load relative promoter activity
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

# Now we will pre-process feature count matrix here only, why here?? -->> because we have to take genes that are present in the filtered step till here
count <- read_csv("Featurecounts_combined_fuscc.csv") %>%
  as.data.frame() %>%
  column_to_rownames("Geneid") %>% 
  subset(., select = -1) %>% 
  filter(rownames(.) %in% filtered_df$geneId)
cpm <- cpm(count) # Calculate cpm
is.exprs <- rowSums(cpm>1) >= dim(tumor_metadata)[1]*0.1
count <- count[is.exprs, ]
x <- DGEList(counts=count )
x <- calcNormFactors(x,method = "TMM")
v <- voom(x, plot=F)
vMat <- v$E

# Batch adjustment
vMat_2 <- sva::ComBat(vMat, batch=metadata$batch)
vMat_2 <- vMat_2[, tumor_metadata$Run] #removing normal sample
p_info <- filtered_df[filtered_df$geneId %in% row.names(vMat_2),]

# Reading Absolute promoter activity
Ab_counts <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)
Ab_counts <- Ab_counts[rownames(Ab_counts) %in% p_info$promoterId,]
colnames(Ab_counts) <- gsub("_SJ.out", "", colnames(Ab_counts))
Ab_counts <- merge( p_info, Ab_counts, by.x= 'promoterId', by.y = 'row.names')
rownames(Ab_counts) <- Ab_counts$ID
Ab_counts <- Ab_counts[, c(-1, -2, -3)]

# TAKING SAMPLES THAT HAVE RFS_time_months > 0 --- survival data we are working with
surv_metadata =tumor_metadata[tumor_metadata$RFS_time_Months > 0,]
surv_metadata=surv_metadata[,c("Project_ID", "Run", "batch", "Tumor_Normal", "Tumor_Normal1", "RFS_Status", "RFS_time_Days", "RFS_time_Months")]

# Making sure that the sample order is same in Absolute promoter activity matrix and gene expression matrix 
Ab_counts <- Ab_counts[, surv_metadata$Run]
vMat_2 <- vMat_2[, surv_metadata$Run]
tobetaken <- rowSums(Ab_counts > 0) >= dim(tumor_metadata)[1]*0.1
Ab_counts=Ab_counts[tobetaken,]

# Step-01 -> Read the testing dataset
testing_dataset <- read_csv('testing_dataset_sub.csv')
rownames(testing_dataset) <- NULL
testing_dataset <- testing_dataset %>% tibble::column_to_rownames('Run')
testing_dataset <- testing_dataset[,-1]
testing_dataset$Run <- rownames(testing_dataset)

# Step-02 -> In feature count matrix take only those samples that are present in testing dataset
keep = intersect(colnames(vMat_2), testing_dataset$Run )
vMAT_up <- vMat_2[, keep]

# Step-03 -> Similarly, take only samples present in testing datasets for metafile too
keep = intersect(surv_metadata$Run, testing_dataset$Run )
surv_metadata <- surv_metadata[keep,]
dim(surv_metadata) 

# Loading the cut-off values for promoters
load('promoters_cut_off_sub_int_false.RData')

# Perfoming univariate survival analysis on testing dataset samples --using the threshold values calculated using training dataset
perform_analysis_for_matrix <- function(metadata, matrix) {
  mm <- merge(metadata, matrix, by = 0)
  rownames(mm) <- mm$Row.names
  mm <- mm[, -1]
}
prom_matrix <- perform_analysis_for_matrix(surv_metadata, t(Ab_counts_up))
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

# Univariate survival analysis
gene_names <- colnames(mm_finaldf)[9:dim(mm_finaldf)[2]]
gene_p_values1 <- numeric(length = length(gene_names))
for (i in 1:length(gene_names)) {
  gene_expression1 <- mm_finaldf[, i + 8]
  fit <- coxph(Surv(mm_finaldf$RFS_time_Days, mm_finaldf$RFS_Status) ~ gene_expression1,
               data = mm_finaldf)
  gene_p_values1[i] <- summary(fit)$coefficients['gene_expression1low', "Pr(>|z|)"]
  cat("P-value for", "gene", i, ":", gene_p_values1[i], "\n")
}

prom_fdr_p_values = p.adjust(gene_p_values1, method = 'fdr')
result_prom <- data.frame(Gene_names = colnames(mm_finaldf)[c(9:dim(mm_finaldf)[2])], Adj_pvalue = prom_fdr_p_values, PValue = gene_p_values1)
result_prom=result_prom %>% separate(Gene_names,into =c("Promoter", "Gene_names2"),sep = "_", remove = FALSE)
write.csv(result_prom, 'all_prom.csv')

# Fetch significant promoters and save result
prom_sig_result <- filter(result_prom, PValue < 0.05)
write.csv(prom_sig_result, "prom_sig_pvalue.csv")


##---------------------------------------- FINAL PROGNOSTIC AAP CANDIDATE FINDING -------------------------------------------------##
# Identifying the candidate genes ---- that are in testing dataset 
list1 <- read_csv('prom_sig_pvalue.csv')
com_p_sig_g_nonsig_pv = (intersect(list1$Gene_names2, gene_nonsig_result$Gene_names))
write.csv(com_p_sig_g_nonsig_pv, 'com_p_sig_g_nonsig_pv.csv')

# To check intersection between training and testing dataset
list2_TEST <- read_csv('com_p_sig_g_nonsig_pv.csv')
list3_TRAIN <- read_csv('com_p_sig_g_nonsig_pv_v2.csv')
intersect(list2_TEST$x, list3_TRAIN$x)

#####################################################################################################################################
#####################################################################################################################################
