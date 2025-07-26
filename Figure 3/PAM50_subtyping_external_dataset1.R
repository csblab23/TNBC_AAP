## Code to run PAM50 subtyping on the External Dataset 1 ##

# Load required libraries
library(dplyr)
library(tidyverse)
library(data.table)
library(genefu)
library(xtable)
library(rmeta)
library(Biobase)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(caret)
library(gprofiler2)
library(genekitr)

set.seed(123)

# NOTE: Since for this dataset we did not have subtyping information so we run PAM50 using the NON-TNBCs of the same dataset

# Read the metadata
meta <- as.data.frame(fread('../metadata_42tnbc_21normal.txt'))
keep <- 'ER+ Breast Cancer Primary Tumor'
meta_nontnbc <- meta[meta$tissue %in% keep,]

# Save the NON-TNBC metadata file
write.csv(meta_nontnbc,'NonTNBC_42_metadata.csv')

# Filter out the 3 samples which failed the QC analysis:
remove <- c('SRR1313117','SRR1313118','SRR1313119')
meta_nontnbc_sel <- meta_nontnbc[!meta_nontnbc$Run %in% remove,]
rownames(meta_nontnbc_sel) <- NULL
write.csv(meta_nontnbc_sel,'NonTNBC_39_metadata.csv')

meta_nontnbc_sel2 <- meta_nontnbc_sel[,c('Run','tissue')]  # subset required columns

meta_nontnbc_sel2 <- meta_nontnbc_sel2 %>%
  mutate(Tumor = case_when(
    tissue == "ER+ Breast Cancer Primary Tumor" ~ "NON_TNBC",
  ))

# Read the metadata for TNBCs and merge it with above df for NON-TNBCs:
meta_tnbc <- as.data.frame(read_csv("meta_sel_61samples.csv"))
meta_tnbc <- meta_tnbc[,-1]
meta_tnbc <- meta_tnbc[meta_tnbc$tumor_normal=="TNBC",]

# Merge with Non tnbc metadata:
colnames(meta_tnbc)[3] <- 'Tumor'
all_tumor_metadata <- rbind(meta_tnbc,meta_nontnbc_sel2)
write.csv(all_tumor_metadata,"All_tumor_metadata.csv")  # this is the combined metadata of TNBC + NON-TNBC

####----

## Now perform PAM50 subtyping:--------------->>
# we have 40 TNBCs and 39 NON-TNBCs
# Read the metadata for all tumors:
metadata <- as.data.frame(read_csv("All_tumor_metadata.csv"))
metadata <- metadata[,-1]

# Read the FPKM matrix:
# For NONTNBC
exp_mat_nontnbc=as.data.frame(read_csv("DF_39NonTNBC_FPKM.csv"))
exp_mat_nontnbc<-exp_mat_nontnbc[,-1]
# for TNBC
exp_mat_tnbc <- as.data.frame(read_csv('DF_FPKM_61samples.csv'))
exp_mat_tnbc<-exp_mat_tnbc[,-1]

# Merge to combine 
merged_mat <- merge(exp_mat_tnbc,exp_mat_nontnbc,by='gene_id')
colnames(merged_mat) <- gsub(".genes.results","",colnames(merged_mat))
merged_mat<- tibble::column_to_rownames(merged_mat,var = 'gene_id')

# Keep the same samples as in the metadata:
merged_mat=merged_mat[,colnames(merged_mat) %in% metadata$Run]
merged_mat <- merged_mat[,metadata$Run]

## Convert the ensembl ids to gene names:
all_gene_ids_2 <- rownames(merged_mat)
all_gene_ids_2 <- gsub("\\..*", "", all_gene_ids_2)
genes <- gprofiler2::gconvert(all_gene_ids_2,organism = 'hsapiens', target = 'ENSG')
gene_names=genes[,c(2,5)]  # fetching just the id and gene name columns
rownames(merged_mat) = gsub("\\..*", "", rownames(merged_mat))
rownames(gene_names) = gene_names$input

# Merge to get the final matrix having gene names:
final_fpkm_matrix = merge(gene_names, merged_mat, by = 0)
final_fpkm_matrix <- final_fpkm_matrix[,-1]

# Check for duplicate gene names:
table(duplicated(final_fpkm_matrix$name))
# Remove duplicates by variance:
e_g=final_fpkm_matrix[,1:2]
e_g$Variance <- apply(final_fpkm_matrix[,3:dim(final_fpkm_matrix)[2]], 1, var)
result <- e_g %>%
  group_by(name) %>%
  arrange(desc(Variance)) %>%
  dplyr::slice(1) %>%
  ungroup()

# Keep the same ensembl ids in the matrix:
final_fpkm_matrix=final_fpkm_matrix[final_fpkm_matrix$input %in% result$input,]
rownames(final_fpkm_matrix)=final_fpkm_matrix$name
final_fpkm_matrix=final_fpkm_matrix[,3:dim(final_fpkm_matrix)[2]]

# Calculate log2(fpkm+1):
numeric_cols <- sapply(final_fpkm_matrix, is.numeric)
final_fpkm_transformed <- log2(final_fpkm_matrix + 1)

# Now separate out TNBC and NON-TNBCs:-----------------
meta_tnbc <- metadata[metadata$Tumor=="TNBC",]
rownames(meta_tnbc)<-meta_tnbc$Run
##
meta_nontnbc <- metadata[metadata$Tumor=="NON_TNBC",]
rownames(meta_nontnbc) <- meta_nontnbc$Run

# Preparing separate matrices for each:
final_tnbc_mat <- final_fpkm_transformed[,colnames(final_fpkm_transformed) %in% rownames(meta_tnbc)]
final_tnbc_mat <- final_tnbc_mat[,rownames(meta_tnbc)] 

##

final_nontnbc_mat <- final_fpkm_transformed[,colnames(final_fpkm_transformed) %in% rownames(meta_nontnbc)]
final_nontnbc_mat <- final_nontnbc_mat[,rownames(meta_nontnbc)] 


## NOW split out samples here----------------------->>
# Three rounds of subtyping will be done:

# Function to rename genes to match PAM50
rename_genes <- function(df) {
  replacements <- c("NUF2" = "CDCA1", "NDC80" = "KNTC2", "ORC6" = "ORC6L")
  for (old in names(replacements)) {
    new <- replacements[old]
    colnames(df) <- gsub(old, new, colnames(df))
    rownames(df) <- gsub(old, new, rownames(df))
    if ("probe" %in% colnames(df)) {
      df$probe <- gsub(old, new, df$probe)
    }
  }
  return(df)
}

# Function to prepare annotation
prepare_annotations <- function(annotation_data) {
  colnames(annotation_data)[1] <- "ensembl"
  id <- annotation_data$ensembl
  genes <- transId(id, transTo = c("ens", "ent"))
  merged <- merge(annotation_data, genes, by = "ensembl")
  clean <- merged[complete.cases(merged$entrezid), ]
  
  # Keep highest variance per entrez and gene name
  clean <- clean %>%
    group_by(entrezid) %>% arrange(desc(Variance)) %>% slice(1) %>% ungroup() %>%
    group_by(name) %>% arrange(desc(Variance)) %>% slice(1) %>% ungroup()
  
  colnames(clean)[2] <- "probe"
  colnames(clean)[5] <- "EntrezGene.ID"
  return(clean)
}

# Function to run PAM50 subtyping
run_pam50 <- function(expression_matrix, annotation, round_id) {
  data(pam50)
  pam50_genes <- rownames(pam50$centroids)
  
  # Filter matrix and annotation
  expression_matrix <- expression_matrix[, colnames(expression_matrix) %in% annotation$probe]
  annotation <- annotation[annotation$probe %in% pam50_genes, ]
  
  # Rename genes
  annotation <- rename_genes(annotation)
  colnames(expression_matrix) <- gsub(".*\\.", "", colnames(expression_matrix))  # cleanup if needed
  
  # Final subset
  counts <- expression_matrix[, colnames(expression_matrix) %in% annotation$probe]
  rownames(annotation) <- annotation$probe
  
  # Rescale expression
  rescale_matrix <- function(mat, q = 0.05) {
    apply(mat, 2, function(x) rescale(x, q = q))
  }
  counts_scaled <- rescale_matrix(counts)
  
  # PAM50 prediction
  preds <- molecular.subtyping(sbt.model = "pam50",
                               data = counts_scaled,
                               annot = as.matrix(annotation),
                               do.mapping = TRUE)
  
  preds <- as.data.frame(preds)
  preds <- tibble::rownames_to_column(preds, var = "Sample")
  write.csv(preds, paste0("pam50_prediction_genefu_set", round_id, ".csv"), row.names = FALSE)
  return(preds)
}

# Setup: Split TNBC sets
tnbc_sets <- list(
  final_tnbc_mat[, 1:13],
  final_tnbc_mat[, 14:26],
  final_tnbc_mat[, 27:ncol(final_tnbc_mat)]
)

## Main Loop:--------------------------->>
all_preds <- list()
for (i in seq_along(tnbc_sets)) {
  cat("Running Round", i, "\n")
  
  mat <- merge(final_nontnbc_mat, tnbc_sets[[i]], by = 0)
  mat <- column_to_rownames(mat, "Row.names")
  mat <- t(mat)
  
  # Annotation prep (only once, or reload from `result` if different each time)
  annot_clean <- prepare_annotations(result)
  
  # Run PAM50
  preds <- run_pam50(mat, annot_clean, i)
  all_preds[[i]] <- preds
}

# Merge PAM50 results
metadata <- metadata %>% column_to_rownames("Run")
all_tnbc_preds <- lapply(seq_along(all_preds), function(i) {
  df <- all_preds[[i]][, 1:2]  # Sample + subtype
  rownames(df) <- df$Sample
  merged <- merge(df, metadata, by = 0)
  tnbc_only <- merged[merged$Tumor == "TNBC", c("Sample", "subtype", "Tumor")]
  return(tnbc_only)
})

# Save the results
final_tnbc_pam50 <- bind_rows(all_tnbc_preds)
write.csv(final_tnbc_pam50, "TNBC_PAM50results.csv", row.names = FALSE)

#####################################################################################
#####################################################################################
