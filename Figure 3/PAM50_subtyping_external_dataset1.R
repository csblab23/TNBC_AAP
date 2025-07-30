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

# Read the metadata for all tumors:
metadata <- as.data.frame(read_csv("All_tumor_metadata.csv"))
metadata <- metadata[,-1]

# Read the FPKM matrix:
# for NONTNBC
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
# 1. 39 nontnbc + 13 tnbc
# 2. 39 nontnbc + 13 tnbc
# 3. 39 nontnbc + 14 tnbc

tnbc_set1 <- final_tnbc_mat[,1:13]
tnbc_set2 <- final_tnbc_mat[,14:26]
tnbc_set3 <- final_tnbc_mat[,27:ncol(final_tnbc_mat)] 

#######

# Preparing df for round 1:
mat_set1 <- merge(final_nontnbc_mat,tnbc_set1,by=0)
mat_set1 <- tibble::column_to_rownames(mat_set1,var = 'Row.names')

# Preparing gene id data so as to convert ens ids to entrez ids:
annota_data=result
colnames(annota_data)[1] ="ensembl"
id <- annota_data$ensembl
genes <- transId(id, transTo = c("ens", "ent"))

# Merge
all_annot=merge(annota_data,genes,by="ensembl")

# Removing those rows that have NA in entrez id
all_annot_clean <- all_annot[complete.cases(all_annot$entrezid), ]

# Again keeping only those entrez ids that have high variance 
all_annot_clean_2 <- all_annot_clean %>%
  group_by(entrezid) %>%
  arrange(desc(Variance)) %>%
  dplyr::slice(1) %>%
  ungroup()

all_annot_clean_final <- all_annot_clean_2 %>%
  group_by(name) %>%
  arrange(desc(Variance)) %>%
  dplyr::slice(1) %>%
  ungroup()

mat_set1=t(mat_set1) # transpose counts matrix
colnames(all_annot_clean_final)[2]="probe"
colnames(all_annot_clean_final)[5]="EntrezGene.ID"

data(pam50)
data(pam50.scale)
data(pam50.robust)
pam50_genes <- as.data.frame(pam50[["centroids"]])
mat_set1=mat_set1[,colnames(mat_set1) %in% all_annot_clean_final$probe]

all_annot_clean_final=as.data.frame(all_annot_clean_final)
rownames(all_annot_clean_final) <- all_annot_clean_final$probe
counts <- mat_set1

# Renaming the genes:
colnames(counts) <- ifelse(grepl("NUF2", colnames(counts)), "CDCA1", colnames(counts))
colnames(counts) <- ifelse(grepl("NDC80", colnames(counts)), "KNTC2", colnames(counts))
colnames(counts) <- ifelse(grepl("ORC6", colnames(counts)), "ORC6L", colnames(counts))

rownames(all_annot_clean_final) <- ifelse(grepl("NUF2", rownames(all_annot_clean_final)), "CDCA1", rownames(all_annot_clean_final))
rownames(all_annot_clean_final) <- ifelse(grepl("NDC80", rownames(all_annot_clean_final)), "KNTC2", rownames(all_annot_clean_final))
rownames(all_annot_clean_final) <- ifelse(grepl("ORC6", rownames(all_annot_clean_final)), "ORC6L", rownames(all_annot_clean_final))

all_annot_clean_final$probe <- ifelse(grepl("NUF2", all_annot_clean_final$probe), "CDCA1", all_annot_clean_final$probe)
all_annot_clean_final$probe <- ifelse(grepl("NDC80", all_annot_clean_final$probe), "KNTC2", all_annot_clean_final$probe)
all_annot_clean_final$probe <- ifelse(grepl("ORC6", all_annot_clean_final$probe), "ORC6L", all_annot_clean_final$probe)

all_annot_clean_final <- all_annot_clean_final[all_annot_clean_final$probe %in% rownames(pam50_genes),]
all_annot_clean_final=as.matrix(all_annot_clean_final)
counts_final <- counts[,colnames(counts) %in% rownames(pam50_genes)]

# Subtyping
rescale_matrix <- function(mat, q = 0.05) {
  apply(mat, 2, function(x) rescale(x, q = q))
}
rescaled_matrix <- rescale_matrix(counts_final)
PAM50Preds <- molecular.subtyping(sbt.model = "pam50",data = rescaled_matrix,
                                  annot = all_annot_clean_final,do.mapping = TRUE)

# Get sample counts pertaining to each subtype
PAM50Preds <- as.data.frame(PAM50Preds)
PAM50Preds <- tibble::rownames_to_column(PAM50Preds,var="Sample")
write.csv(PAM50Preds,"pam50_prediction_genefu_set1.csv")

#### Round 2:
mat_set2 <- merge(final_nontnbc_mat,tnbc_set2,by=0)
mat_set2 <- tibble::column_to_rownames(mat_set2,var = 'Row.names')

# Preparing gene id data so as to convert ens ids to entrez ids:
annota_data=result
colnames(annota_data)[1] ="ensembl"
id <- annota_data$ensembl
genes <- transId(id, transTo = c("ens", "ent"))

# Merge
all_annot=merge(annota_data,genes,by="ensembl")

# Removing those rows that has na in entrez id
all_annot_clean <- all_annot[complete.cases(all_annot$entrezid), ]

# Again keeping only those entrez ids that has high variance 
all_annot_clean_2 <- all_annot_clean %>%
  group_by(entrezid) %>%
  arrange(desc(Variance)) %>%
  dplyr::slice(1) %>%
  ungroup()

all_annot_clean_final <- all_annot_clean_2 %>%
  group_by(name) %>%
  arrange(desc(Variance)) %>%
  dplyr::slice(1) %>%
  ungroup()

mat_set2=t(mat_set2) # transpose counts matrix
colnames(all_annot_clean_final)[2]="probe"
colnames(all_annot_clean_final)[5]="EntrezGene.ID"

data(pam50)
data(pam50.scale)
data(pam50.robust)
pam50_genes <- as.data.frame(pam50[["centroids"]])
mat_set2=mat_set2[,colnames(mat_set2) %in% all_annot_clean_final$probe]

all_annot_clean_final=as.data.frame(all_annot_clean_final)
rownames(all_annot_clean_final) <- all_annot_clean_final$probe
counts <- mat_set2

# Renaming the genes:
colnames(counts) <- ifelse(grepl("NUF2", colnames(counts)), "CDCA1", colnames(counts))
colnames(counts) <- ifelse(grepl("NDC80", colnames(counts)), "KNTC2", colnames(counts))
colnames(counts) <- ifelse(grepl("ORC6", colnames(counts)), "ORC6L", colnames(counts))

rownames(all_annot_clean_final) <- ifelse(grepl("NUF2", rownames(all_annot_clean_final)), "CDCA1", rownames(all_annot_clean_final))
rownames(all_annot_clean_final) <- ifelse(grepl("NDC80", rownames(all_annot_clean_final)), "KNTC2", rownames(all_annot_clean_final))
rownames(all_annot_clean_final) <- ifelse(grepl("ORC6", rownames(all_annot_clean_final)), "ORC6L", rownames(all_annot_clean_final))

all_annot_clean_final$probe <- ifelse(grepl("NUF2", all_annot_clean_final$probe), "CDCA1", all_annot_clean_final$probe)
all_annot_clean_final$probe <- ifelse(grepl("NDC80", all_annot_clean_final$probe), "KNTC2", all_annot_clean_final$probe)
all_annot_clean_final$probe <- ifelse(grepl("ORC6", all_annot_clean_final$probe), "ORC6L", all_annot_clean_final$probe)

all_annot_clean_final <- all_annot_clean_final[all_annot_clean_final$probe %in% rownames(pam50_genes),]
all_annot_clean_final=as.matrix(all_annot_clean_final)
counts_final <- counts[,colnames(counts) %in% rownames(pam50_genes)]

# Subtyping
rescale_matrix <- function(mat, q = 0.05) {
  apply(mat, 2, function(x) rescale(x, q = q))
}
rescaled_matrix <- rescale_matrix(counts_final)
PAM50Preds <- molecular.subtyping(sbt.model = "pam50",data = rescaled_matrix,
                                  annot = all_annot_clean_final,do.mapping = TRUE)


# Get sample counts pertaining to each subtype
PAM50Preds <- as.data.frame(PAM50Preds)
PAM50Preds <- tibble::rownames_to_column(PAM50Preds,var="Sample")
write.csv(PAM50Preds,"pam50_prediction_genefu_set2.csv")

#### Round 3:
mat_set3 <- merge(final_nontnbc_mat,tnbc_set3,by=0)
mat_set3 <- tibble::column_to_rownames(mat_set3,var = 'Row.names')

# Preparing gene id data so as to convert ens ids to entrez ids:
annota_data=result
colnames(annota_data)[1] ="ensembl"
id <- annota_data$ensembl
genes <- transId(id, transTo = c("ens", "ent"))

# Merge
all_annot=merge(annota_data,genes,by="ensembl")

# Removing those rows that has na in entrez id
all_annot_clean <- all_annot[complete.cases(all_annot$entrezid), ]

# Again keeping only those entrez ids that has high variance 
all_annot_clean_2 <- all_annot_clean %>%
  group_by(entrezid) %>%
  arrange(desc(Variance)) %>%
  dplyr::slice(1) %>%
  ungroup()

all_annot_clean_final <- all_annot_clean_2 %>%
  group_by(name) %>%
  arrange(desc(Variance)) %>%
  dplyr::slice(1) %>%
  ungroup()

mat_set3=t(mat_set3) # transpose counts matrix
colnames(all_annot_clean_final)[2]="probe"
colnames(all_annot_clean_final)[5]="EntrezGene.ID"

data(pam50)
data(pam50.scale)
data(pam50.robust)
pam50_genes <- as.data.frame(pam50[["centroids"]])
mat_set3=mat_set3[,colnames(mat_set3) %in% all_annot_clean_final$probe]

all_annot_clean_final=as.data.frame(all_annot_clean_final)
rownames(all_annot_clean_final) <- all_annot_clean_final$probe
counts <- mat_set3

# Renaming the genes:
colnames(counts) <- ifelse(grepl("NUF2", colnames(counts)), "CDCA1", colnames(counts))
colnames(counts) <- ifelse(grepl("NDC80", colnames(counts)), "KNTC2", colnames(counts))
colnames(counts) <- ifelse(grepl("ORC6", colnames(counts)), "ORC6L", colnames(counts))

rownames(all_annot_clean_final) <- ifelse(grepl("NUF2", rownames(all_annot_clean_final)), "CDCA1", rownames(all_annot_clean_final))
rownames(all_annot_clean_final) <- ifelse(grepl("NDC80", rownames(all_annot_clean_final)), "KNTC2", rownames(all_annot_clean_final))
rownames(all_annot_clean_final) <- ifelse(grepl("ORC6", rownames(all_annot_clean_final)), "ORC6L", rownames(all_annot_clean_final))

all_annot_clean_final$probe <- ifelse(grepl("NUF2", all_annot_clean_final$probe), "CDCA1", all_annot_clean_final$probe)
all_annot_clean_final$probe <- ifelse(grepl("NDC80", all_annot_clean_final$probe), "KNTC2", all_annot_clean_final$probe)
all_annot_clean_final$probe <- ifelse(grepl("ORC6", all_annot_clean_final$probe), "ORC6L", all_annot_clean_final$probe)

all_annot_clean_final <- all_annot_clean_final[all_annot_clean_final$probe %in% rownames(pam50_genes),]
all_annot_clean_final=as.matrix(all_annot_clean_final)
counts_final <- counts[,colnames(counts) %in% rownames(pam50_genes)]

# Subtyping
rescale_matrix <- function(mat, q = 0.05) {
  apply(mat, 2, function(x) rescale(x, q = q))
}
rescaled_matrix <- rescale_matrix(counts_final)
PAM50Preds <- molecular.subtyping(sbt.model = "pam50",data = rescaled_matrix,
                                  annot = all_annot_clean_final,do.mapping = TRUE)

# Get sample counts pertaining to each subtype
PAM50Preds <- as.data.frame(PAM50Preds)
PAM50Preds <- tibble::rownames_to_column(PAM50Preds,var="Sample")
write.csv(PAM50Preds,"pam50_prediction_genefu_set3.csv")

##--

## Merging the output results now and making one output of pam50 subtyping result:
output1 <- as.data.frame(read_csv('pam50_prediction_genefu_set1.csv'))
output1 <- output1[,-1]
output1 <- output1[,1:2]
rownames(output1) <- output1$Sample
rownames(metadata)<-metadata$Run
output1_merged <- merge(output1,metadata,by=0)


# Separating tnbcs from here:
output1_tnbc <- output1_merged[output1_merged$Tumor=="TNBC",]
output1_tnbc <- output1_tnbc[,c(2,3,6)] # subsetting required columns

##---
output2 <- as.data.frame(read_csv('pam50_prediction_genefu_set2.csv'))
output2 <- output2[,-1]
output2 <- output2[,1:2]
rownames(output2) <- output2$Sample
rownames(metadata)<-metadata$Run
output2_merged <- merge(output2,metadata,by=0)


# Separating tnbcs from here:
output2_tnbc <- output2_merged[output2_merged$Tumor=="TNBC",]
output2_tnbc <- output2_tnbc[,c(2,3,6)] # subsetting required columns

##--
output3 <- as.data.frame(read_csv('pam50_prediction_genefu_set3.csv'))
output3 <- output3[,-1]
output3 <- output3[,1:2]
rownames(output3) <- output3$Sample

rownames(metadata)<-metadata$Run
output3_merged <- merge(output3,metadata,by=0)

# Separating tnbcs from here:
output3_tnbc <- output3_merged[output3_merged$Tumor=="TNBC",]
output3_tnbc <- output3_tnbc[,c(2,3,6)] # subsetting required columns

###########

# Combine all the three TNBC outputs & Save the results:
TNBC_pam50_results <- rbind(output1_tnbc,output2_tnbc,output3_tnbc)
rownames(TNBC_pam50_results) <- NULL
write.csv(TNBC_pam50_results,'TNBC_PAM50results.csv')

#####################################################################################
#####################################################################################