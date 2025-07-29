## Code for Finding BL1-specific AAPs in FUSCC and validating in external dataset 1 ##

## STEPS INCLUDED:
# 1. Find BL1-specific AAPs in FUSCC dataset
# 2. Find BL1-specific AAPs in external dataset 1 (referred to as the validation dataset)
# 3. Intersect to find common AAPs

####---------------------------------------------------------------------------------####

# Load required libraries
library(readxl)
library(qtl2)
library(dplyr)
library(tidyverse)
library(data.table)

# Find BL1-specific AAPs in FUSCC dataset:---->>

# DEG between BL1 vs REST 
# Load the metafile  
metadata <- as.data.frame(read_csv("fuscc_merged_meta.csv"))
metadata=tibble::column_to_rownames(metadata,var = '...1')
subset_metadata <- metadata

# Read the feature-counts matrix 
count <- as.data.frame(read_csv("Featurecounts_combined_fuscc.csv"))
count=count[,-1] # remove extra s.no. column
count <- tibble::column_to_rownames(count, var = "Geneid")

# Keep same samples in count matrix as in metadata:
count <- count[,colnames(count) %in% subset_metadata$Run]
count <- count[,subset_metadata$Run]  # make order of samples same 

######

# Now calculate cpm  
cpm <- cpm(count)
is.exprs <- rowSums(cpm>1) >= 20 # Apply filtering criteria
counts2 <- count[is.exprs, ]

# Prepare DGEList and voom normalize the matrix
x <- DGEList(counts=counts2)
x <- calcNormFactors(x, method='TMM')
v <- voom(x, plot=T)
vMat <- v$E

# Adjust batch by sva combat here
subset_metadata$batch=as.factor(subset_metadata$batch)
exprMat <- sva::ComBat(vMat, batch=subset_metadata$batch)

####----

# Run the Wilcoxon rank-sum test for each gene
count_norm <- as.data.frame(exprMat)

# Now make the conditions column 
conditions <- subset_metadata[,c("Run","bl1_vs_rest")]
count_norm <- count_norm[,conditions$Run]
conditions <- t(conditions)
conditions <- conditions[2,]
conditions<-factor(conditions)

# Make sure order of samples is same:
pvalues <- sapply(1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
  p=wilcox.test(gene~conditions, data)$p.value
  return(p)
})

fdr=p.adjust(pvalues,method = "fdr")

# Calculate the fold-change for each gene
conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))] # bl1
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))] # rest
foldChanges=log2(rowMeans(dataCon1)/rowMeans(dataCon2))

# Output results
outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)
write.csv(outRst,"DEG_wilcox_BL1_vs_REST.csv")

# Fetch upregulated genes
upregulated_genes <- outRst[outRst$log2foldChange > 1 & outRst$pValues < 0.05,]
write.csv(upregulated_genes,"upreg_deg_BL1_vs_REST.csv")

# Fetch downregulated genes
downregulated_genes <- outRst[outRst$log2foldChange < -1 & outRst$pValues < 0.05,]
write.csv(downregulated_genes,"downreg_deg_BL1_vs_REST.csv")

##--

# Read the result2 file obtained by running proactiv
result2 <- readRDS("result2_updated_subype_fuscc.rds")

# Absolute promoter activity:-------------------------------->>
abs_pa <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)  
colnames(abs_pa) <- gsub("_SJ.out","",colnames(abs_pa))

# Filtering
result2_promoter = rowData(result2) 
filtered_df <- as.data.frame(result2_promoter)
filtered_df <- filtered_df[, c('promoterId','geneId','internalPromoter')]

# Filter 1: keep internal promoters = FALSE ------------------------------
filtered_df <- filtered_df[filtered_df$internalPromoter == "FALSE",]

# Filter 2: keep the genes which have >=2 prom ---------------------------
# for this we will fetch the rowData
filtered_df <- filtered_df %>%
  group_by(geneId) %>%
  filter(n() >= 2) %>%
  mutate(ID = paste('pr', promoterId, '_', geneId, sep='')) 

# Keep the same promoter ids in the absolute counts matrix
abs_pa_filtered <- abs_pa[rownames(abs_pa) %in% filtered_df$promoterId,]

# Run the Wilcoxon rank-sum test for each gene
count_norm <- as.data.frame(abs_pa_filtered)

# Now make the conditions column
conditions <- subset_metadata[,c("Run","bl1_vs_rest")]
count_norm <- count_norm[,conditions$Run]
conditions <- t(conditions)
conditions <- conditions[2,]
conditions<-factor(conditions)

# Make sure to have the order of samples same:
count_norm <- count_norm[,subset_metadata$Run]

pvalues <- sapply(1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
  p=wilcox.test(gene~conditions, data)$p.value
  return(p)
})

fdr=p.adjust(pvalues,method = "fdr")

# Calculate the fold-change for each gene
conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))] # BL1
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))] # rest
foldChanges=log2(rowMeans(dataCon1)/rowMeans(dataCon2))

# Output results 
outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)
write.csv(outRst,"DRP_wilcox_AbsPA.csv")

# Fetch upregulated DRPs
upregulated_genes <- outRst[outRst$log2foldChange > 1.2 & outRst$pValues < 0.05,]
write.csv(upregulated_genes,"upreg_drp_absPA.csv")

# Fetch downregulated DRPs
downregulated_genes <- outRst[outRst$log2foldChange < -1.2 & outRst$pValues < 0.05,]
write.csv(downregulated_genes,"down_drp_absPA.csv")

##--

# Relative promoter activity:-------------------------------->>
rel_pa <- as.data.frame(result2@assays@data@listData$relativePromoterActivity)  
colnames(rel_pa) <- gsub("_SJ.out","",colnames(rel_pa))

# Filtering
result2_promoter = rowData(result2) 
filtered_df <- as.data.frame(result2_promoter)
filtered_df <- filtered_df[, c('promoterId','geneId','internalPromoter')]

# Filter 1: keep internal promoters = FALSE ------------------------------
filtered_df <- filtered_df[filtered_df$internalPromoter == "FALSE",]

# Filter 2: keep the genes which have >=2 prom ---------------------------
# for this we will fetch the rowData
filtered_df <- filtered_df %>%
  group_by(geneId) %>%
  filter(n() >= 2) %>%
  mutate(ID = paste('pr', promoterId, '_', geneId, sep='')) 

# Keep the same promoter ids in the absolute counts matrix
rel_pa_filtered <- rel_pa[rownames(rel_pa) %in% filtered_df$promoterId,]

# There are many NA values so just take complete cases:
rel_pa_filtered2 <- rel_pa_filtered[complete.cases(rel_pa_filtered),]

##--

# Run the Wilcoxon rank-sum test for each gene
count_norm <- as.data.frame(rel_pa_filtered2)

# Now make the conditions column 
conditions <- subset_metadata[,c("Run","bl1_vs_rest")]
count_norm <- count_norm[,conditions$Run]
conditions <- t(conditions)
conditions <- conditions[2,]
conditions<-factor(conditions)

# Make sure to have the order of samples same:
count_norm <- count_norm[,subset_metadata$Run]

pvalues <- sapply(1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
  p=wilcox.test(gene~conditions, data)$p.value
  return(p)
})

fdr=p.adjust(pvalues,method = "fdr")

# Calculate the fold-change for each gene
conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))] # bl1
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))] # rest
foldChanges=log2(rowMeans(dataCon1)/rowMeans(dataCon2))

# Output results 
outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)
write.csv(outRst,"DRP_wilcox_RelPA.csv")

# Fetch upregulated DRPs
upregulated_genes <- outRst[outRst$log2foldChange > 1.2 & outRst$pValues < 0.05,]
write.csv(upregulated_genes,"upreg_drp_RelPA.csv")

# Fetch downregulated DRPs
downregulated_genes <- outRst[outRst$log2foldChange < -1.2 & outRst$pValues < 0.05,]
write.csv(downregulated_genes,"down_drp_RelPA.csv")

#################################---------------------------------------

## FETCH THE FUSCC BL1-specific AAPS:

genes_out <- as.data.frame(read_csv('DEG_wilcox_BL1_vs_REST.csv'))
abspa_out <- as.data.frame(read_csv('DRP_wilcox_AbsPA.csv'))
relpa_out <- as.data.frame(read_csv('DRP_wilcox_RelPA.csv'))

# Fetch nonsignificant genes
colnames(genes_out)[1] <- 'geneID'
genes_out_nonsig <- genes_out[genes_out$pValues>0.05,]  

# Keep nonsignificant genes only
filtered_df2 <- filtered_df[filtered_df$geneId%in%genes_out_nonsig$geneID,]

# Keep DRPs significant by absolute promoter activity
colnames(abspa_out)[1] <- 'promId'
abspa_out2 <- abspa_out[abspa_out$pValues<0.05,]

# Keep DRPs significant by relative promoter activity
colnames(relpa_out)[1] <- 'promId'
relpa_out2 <- relpa_out[relpa_out$pValues<0.05,]

# Find DRPs significant by both absolute and relative promoter activity
sig_p <- intersect(abspa_out2$promId,relpa_out2$promId)
filtered_df3 <- filtered_df2[filtered_df2$promoterId%in%sig_p,]
write.csv(filtered_df3,'final_BL1_AAPs_fuscc.csv')

######################################################################################
######################################################################################

# Find BL1-specific AAPs in Validation dataset

# Load the metafile  
metadata <- as.data.frame(read_csv("val_merged_meta.csv"))
metadata <- tibble::column_to_rownames(metadata,var = "...1")

# Find DEG between BL1 vs REST ---->>
# Read the feature-counts matrix 
count <- as.data.frame(read_csv("Raw_counts_61samples.csv"))
count=count[,-1] # remove extra s.no. column
count <- tibble::column_to_rownames(count, var = "Geneid")

# Keep just same samples in count matrix as in meta:
count <- count[,colnames(count) %in% metadata$Run]
count <- count[,metadata$Run]  # make order of samples same 

######

# Now calculate cpm  
cpm <- cpm(count)
is.exprs <- rowSums(cpm>1) >= 20 # Apply filtering criteria
counts2 <- count[is.exprs, ]

# Prepare DGEList and voom normalize the matrix
x <- DGEList(counts=counts2)
x <- calcNormFactors(x, method='TMM')
v <- voom(x, plot=T)
vMat <- v$E

####----

# Run the Wilcoxon rank-sum test for each gene
count_norm <- as.data.frame(vMat)

# Now make the conditions column 
rownames(metadata) <- metadata$Run
conditions <- metadata[,c("Run","bl1_vs_rest")]
count_norm <- count_norm[,conditions$Run]
conditions <- t(conditions)
conditions <- conditions[2,]
conditions<-factor(conditions)

# Make sure order of samples is same:
pvalues <- sapply(1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
  p=wilcox.test(gene~conditions, data)$p.value
  return(p)
})

fdr=p.adjust(pvalues,method = "fdr")

# Calculate the fold-change for each gene
conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]  # bl1
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]  # rest
foldChanges=log2(rowMeans(dataCon1)/rowMeans(dataCon2))

# Output results based on FDR threshold
outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)
write.csv(outRst,"DEG_wilcox_bl1_vs_REST.csv")

# Fetch upregulated genes
upregulated_genes <- outRst[outRst$log2foldChange > 1 & outRst$pValues < 0.05,]
write.csv(upregulated_genes,"upreg_deg_bl1_vs_REST.csv")

# Fetch downregulated genes
downregulated_genes <- outRst[outRst$log2foldChange < -1 & outRst$pValues < 0.05,]
write.csv(downregulated_genes,"downreg_deg_bl1_vs_REST.csv")

####----

# 2. Find DRPs between BL1 vs REST

# Read the result2 file obtained by running proactiv on tnbc type4 (UNS removed)
result2 <- readRDS("result2_updated_subype_val.rds")

# Absolute promoter activity:-------------------------------->>
abs_pa <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)  
colnames(abs_pa) <- gsub("_SJ.out","",colnames(abs_pa))

# Filtering
result2_promoter = rowData(result2) 
filtered_df <- as.data.frame(result2_promoter)
filtered_df <- filtered_df[, c('promoterId','geneId','internalPromoter')]

# Filter 1: keep internal promoters = FALSE ------------------------------
filtered_df <- filtered_df[filtered_df$internalPromoter == "FALSE",]

# Filter 2: keep the genes which have >=2 prom ---------------------------
# for this we will fetch the rowData
filtered_df <- filtered_df %>%
  group_by(geneId) %>%
  filter(n() >= 2) %>%
  mutate(ID = paste('pr', promoterId, '_', geneId, sep='')) 

# Now keep the same promoter ids in the absolute counts matrix
abs_pa_filtered <- abs_pa[rownames(abs_pa) %in% filtered_df$promoterId,]

####----

# Run the Wilcoxon rank-sum test for each gene
count_norm <- as.data.frame(abs_pa_filtered)

# Define conditions
conditions <- metadata[,c("Run","bl1_vs_rest")]
count_norm <- count_norm[,conditions$Run]
conditions <- t(conditions)
conditions <- conditions[2,]
conditions<-factor(conditions)

# Make sure to have the order of samples same:
count_norm <- count_norm[,metadata$Run]

pvalues <- sapply(1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
  p=wilcox.test(gene~conditions, data)$p.value
  return(p)
})

fdr=p.adjust(pvalues,method = "fdr")

# Calculate the fold-change for each gene
conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]  # bl1
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]  # rest
foldChanges=log2(rowMeans(dataCon1)/rowMeans(dataCon2))

# Output results 
outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)
write.csv(outRst,"DRP_wilcox_AbsPA.csv")

# Fetch upregulated DRPs
upregulated_genes <- outRst[outRst$log2foldChange > 1.2 & outRst$pValues < 0.05,]
write.csv(upregulated_genes,"upreg_drp_absPA.csv")

# Fetch downregulated DRPs
downregulated_genes <- outRst[outRst$log2foldChange < -1.2 & outRst$pValues < 0.05,]
write.csv(downregulated_genes,"down_drp_absPA.csv")

##--

# Relative promoter activity:-------------------------------->>
rel_pa <- as.data.frame(result2@assays@data@listData$relativePromoterActivity)  
colnames(rel_pa) <- gsub("_SJ.out","",colnames(rel_pa))

# Filtering
result2_promoter = rowData(result2) 
filtered_df <- as.data.frame(result2_promoter)
filtered_df <- filtered_df[, c('promoterId','geneId','internalPromoter')]

# Filter 1: keep internal promoters = FALSE ------------------------------
filtered_df <- filtered_df[filtered_df$internalPromoter == "FALSE",]

# Filter 2: keep the genes which have >=2 prom ---------------------------
# for this we will fetch the rowData
filtered_df <- filtered_df %>%
  group_by(geneId) %>%
  filter(n() >= 2) %>%
  mutate(ID = paste('pr', promoterId, '_', geneId, sep='')) 

# Now keep the same promoter ids in the absolute counts matrix
rel_pa_filtered <- rel_pa[rownames(rel_pa) %in% filtered_df$promoterId,]

# There are many NA values so just take complete cases:
rel_pa_filtered2 <- rel_pa_filtered[complete.cases(rel_pa_filtered),]

##--

# Run the Wilcoxon rank-sum test for each gene
count_norm <- as.data.frame(rel_pa_filtered2)

# Define conditions  
conditions <- metadata[,c("Run","bl1_vs_rest")]
count_norm <- count_norm[,conditions$Run]
conditions <- t(conditions)
conditions <- conditions[2,]
conditions<-factor(conditions)

# Make sure to have the order of samples same:
count_norm <- count_norm[,metadata$Run]

pvalues <- sapply(1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
  p=wilcox.test(gene~conditions, data)$p.value
  return(p)
})

fdr=p.adjust(pvalues,method = "fdr")

# Calculate the fold-change for each gene
conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]  # bl1
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]  # rest
foldChanges=log2(rowMeans(dataCon1)/rowMeans(dataCon2))

# Output results 
outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)
write.csv(outRst,"DRP_wilcox_RelPA.csv")

# Fetch upregulated DRPs
upregulated_genes <- outRst[outRst$log2foldChange > 1.2 & outRst$pValues < 0.05,]
write.csv(upregulated_genes,"upreg_drp_RelPA.csv")

# Fetch downregulated DRPs
downregulated_genes <- outRst[outRst$log2foldChange < -1.2 & outRst$pValues < 0.05,]
write.csv(downregulated_genes,"down_drp_RelPA.csv")

##--

#### Applying filters to fetch AAPs:
genes_out <- as.data.frame(read_csv('DEG_wilcox_bl1_vs_REST.csv'))
abspa_out <- as.data.frame(read_csv('DRP_wilcox_AbsPA.csv'))
relpa_out <- as.data.frame(read_csv('DRP_wilcox_RelPA.csv'))

# Fetch nonsignificant genes:
colnames(genes_out)[1] <- 'geneID'
genes_out_nonsig <- genes_out[genes_out$pValues>0.05,]

# Filter and keep only nonsignificant genes:
filtered_df2 <- filtered_df[filtered_df$geneId%in%genes_out_nonsig$geneID,]

# Fetch promoters significant by absolute promoter activity
colnames(abspa_out)[1] <- 'promId'
abspa_out2 <- abspa_out[abspa_out$pValues<0.05,]

# Fetch promoters significant by relative promoter activity
colnames(relpa_out)[1] <- 'promId'
relpa_out2 <- relpa_out[relpa_out$pValues<0.05,]

# Identify promoters significant by both absolute and relative promoter activity
sig_p <- intersect(abspa_out2$promId,relpa_out2$promId)
filtered_df3 <- filtered_df2[filtered_df2$promoterId%in%sig_p,]
write.csv(filtered_df3,'final_BL1_AAPs_val.csv')

#################################################################################
#################################################################################

# Intersect to find common AAPs:

# Load FUSCC BL1-specific AAPs
fuscc_aap <- read_csv('final_BL1_AAPs_fuscc.csv')
fuscc_aap <- fuscc_aap[,-1]

# Load external dataset 1 BL1-specific AAPs
val_aap <- read_csv('final_BL1_AAPs_val.csv')
val_aap <- val_aap[,-1]

# Find common 
intersect(fuscc_aap$promoterId,val_aap$promoterId)
# 1002  1003  9366 10530 10533 10653 24831 25747 30507

#################################################################################
#################################################################################
