## Code to find BASAL-specific AAPs in FUSCC and their validation in external dataset 1 ##

# This script includes the following:
# 1. Find BASAL AAPs in FUSCC
# 2. Find BASAL AAPs in external dataset 1 (regarded as the validation dataset in the code)
# 3. Fetch common AAPs between FUSCC and validation dataset

##-------------------------------------------------------------------------------------------##

# Load required libraries
library(readxl)
library(qtl2)
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
library(edgeR)
library(genekitr)

# Find DEGs and DRPs between BASAL vs NON_BASAL (OTHER) for FUSCC:
# Load the FUSCC metafile  
metadata <- as.data.frame(read_csv("./meta_with_batchinfo.csv"))
metadata=metadata[,-1]  # removing the extra s.no. column
metadata$Tumor_Normal1=metadata$Tumor_Normal
metadata =metadata %>% 
  mutate(Tumor_Normal1 = case_when(Tumor_Normal=="TumorTissue"  ~ "Tumor",
                                   Tumor_Normal=="TumorTissue.2018"  ~ "Tumor",
                                   Tumor_Normal=="TumorTissue.rep"  ~ "Tumor",
                                   Tumor_Normal=="PairedNormalTissue"  ~ "AdjNormal",
                                   Tumor_Normal=="PairedNormalTissue.rep"  ~ "AdjNormal"))
rownames(metadata) <- metadata$Run

# Take all tumor samples only:
metadata <- metadata[metadata$Tumor_Normal1=='Tumor',]
subset_metadata <- metadata[,c('Run','batch','Intrinsic_Subtype')] # just keep the required columns

#####

# Read the feature-counts matrix 
count <- as.data.frame(read_csv("Featurecounts_combined_fuscc.csv"))
count=count[,-1] # remove extra s.no. column
count <- tibble::column_to_rownames(count, var = "Geneid")

# Keep same samples in count matrix as in meta:
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
conditions <- subset_metadata[,c("Run","Intrinsic_Subtype")]
count_norm <- count_norm[,conditions$Run]
conditions <- t(conditions)
conditions <- conditions[2,]
conditions<-factor(conditions)

# make sure order of samples is same:
pvalues <- sapply(1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
  p=wilcox.test(gene~conditions, data)$p.value
  return(p)
})

fdr=p.adjust(pvalues,method = "fdr")

# Calculate the fold-change for each gene
conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))] # basal
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))] # other
foldChanges=log2(rowMeans(dataCon1)/rowMeans(dataCon2))

# Output results based on FDR threshold
outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)

# filter genes here 
write.csv(outRst,"DEG_wilcox_BASAL_vs_OTHER.csv")

# Fetch upregulated DEGs
upregulated_genes <- outRst[outRst$log2foldChange > 1 & outRst$pValues < 0.05,]
write.csv(upregulated_genes,"upreg_deg_BASAL_vs_OTHER.csv")

# Fetch downregulated DEGs
downregulated_genes <- outRst[outRst$log2foldChange < -1 & outRst$pValues < 0.05,]
write.csv(downregulated_genes,"downreg_deg_BASAL_vs_OTHER.csv")

####-----

# For DRP between BASAL vs NON-BASAL (OTHER) using Absolute Promoter Activity and Relative Promoter Activity:

# Read the result2 file obtained by running proactiv
result2 <- readRDS("result2_intrinsic_subtype.rds")

# Part 1: Use Absolute promoter activity
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

####----

# Run the Wilcoxon rank-sum test for each gene
count_norm <- as.data.frame(abs_pa_filtered)

# Now make the conditions column
conditions <- subset_metadata[,c("Run","Intrinsic_Subtype")]
count_norm <- count_norm[,conditions$Run]
conditions <- t(conditions)
conditions <- conditions[2,]
conditions<-factor(conditions)

# make sure to have the order of samples same:
count_norm <- count_norm[,subset_metadata$Run]

pvalues <- sapply(1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
  p=wilcox.test(gene~conditions, data)$p.value
  return(p)
})

fdr=p.adjust(pvalues,method = "fdr")

# Calculate the fold-change for each gene
conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))] # basal
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))] # other
foldChanges=log2(rowMeans(dataCon1)/rowMeans(dataCon2))

# Output results 
outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)

# filter genes here 
write.csv(outRst,"DRP_wilcox_AbsPA.csv")

# Fetch upregulated DRPs
upregulated_genes <- outRst[outRst$log2foldChange > 1.2 & outRst$pValues < 0.05,]
write.csv(upregulated_genes,"upreg_drp_absPA.csv")

# Fetch downregulated DRPs
downregulated_genes <- outRst[outRst$log2foldChange < -1.2 & outRst$pValues < 0.05,]
write.csv(downregulated_genes,"down_drp_absPA.csv")

######------

# Part 2: Use Relative promoter activity
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
conditions <- subset_metadata[,c("Run","Intrinsic_Subtype")]
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
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))] # basal
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))] # other
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

####----

## Identify the Basal-specific AAPS for FUSCC dataset:------------------------->>

# Read the result2 file obtained by running proactiv 
result2 <- readRDS("result2_intrinsic_subtype.rds")

# Remove internal promoters & keep multi-promoter genes only  
result2_promoter = rowData(result2) 
filtered_df <- as.data.frame(result2_promoter)
filtered_df <- filtered_df[, c('promoterId','geneId','internalPromoter')]

# Filter 1: keep internal promoters = FALSE ------------------------------
filtered_df <- filtered_df[filtered_df$internalPromoter == "FALSE",]

# Filter 2: keep the genes which have >=2 prom ---------------------------
filtered_df <- filtered_df %>%
  group_by(geneId) %>%
  filter(n() >= 2) %>%
  mutate(ID = paste('pr', promoterId, '_', geneId, sep='')) 

filtered_df2 <- filtered_df[,c('promoterId', 'geneId')] # subset columns

# Read the output file for DEG:
gene_output <- as.data.frame(read_csv('DEG_wilcox_BASAL_vs_OTHER.csv'))
colnames(gene_output)[1] <- 'Geneid'

# Keep only nonsignificant genes:
gene_output_nonsig <- gene_output[gene_output$pValues > 0.05,]
nonsig_g_filtered_df2 <- filtered_df2[filtered_df2$geneId%in%gene_output_nonsig$Geneid,]

####

# Read the Absolute promoter activity wilcox result:
abspa <- as.data.frame(read_csv('DRP_wilcox_AbsPA.csv'))
colnames(abspa)[1] <- 'promid'
signif_abspa <- abspa[abspa$pValues<0.05,]  # fetch significant DRPs

# Read the Relative promoter activity wilcox result:
relpa <- as.data.frame(read_csv('DRP_wilcox_RelPA.csv'))
colnames(relpa)[1] <- 'promid'
signif_relpa <- relpa[relpa$pValues<0.05,]   # fetch significant DRPs

# Intersect to get the DRPs significant by both absolute and relative promoter activity
signif_prom_all <- intersect(signif_abspa$promid,signif_relpa$promid)

# Keep them in the above dataframe:
nonsig_g_sig_p_filtered_df2 <- nonsig_g_filtered_df2[nonsig_g_filtered_df2$promoterId%in%signif_prom_all,]
write.csv(nonsig_g_sig_p_filtered_df2,"nonsig_g_sig_p_fuscc.csv")

#################################################################################

## Similarly, Script to find BASAL-sepcific AAPs in validation/external dataset 1:

# 1. find DEG between BASAL vs NON_BASAL (OTHER) ------->>
# Load the metafile  
metadata <- as.data.frame(read_csv("TNBC_mergedmeta.csv"))
metadata=tibble::column_to_rownames(metadata,var = '...1')

# Read the subtypes information file:
subtype_meta <- read_csv("subtyping_4_result.csv")
subtype_meta_sel <- subtype_meta[,c('Original_name','tnbctype_4')]  # just keep required columns
colnames(subtype_meta_sel)[1]<-'Run'

# Merge both files
combined <- merge(metadata,subtype_meta_sel,by='Run')

# Read the feature-counts matrix 
count <- as.data.frame(read_csv("Raw_counts_61samples.csv"))
count=count[,-1] # remove extra s.no. column
count <- tibble::column_to_rownames(count, var = "Geneid")

# keep just same samples in count matrix as in meta:
count <- count[,colnames(count) %in% metadata$Run]
count <- count[,metadata$Run]  # make order of samples same 

######

# Now calculate cpm  
cpm <- cpm(count)
is.exprs <- rowSums(cpm>1) >= 5 # Apply filtering criteria
counts2 <- count[is.exprs, ]

# Prepare DGEList and voom normalize the matrix
x <- DGEList(counts=counts2)
x <- calcNormFactors(x, method='TMM')
v <- voom(x, plot=T)
vMat <- v$E

##--

# Run the Wilcoxon rank-sum test for each gene
count_norm <- as.data.frame(vMat)

# Define conditions  
conditions <- metadata[,c("Run","intrinsic_subtype")]
count_norm <- count_norm[,conditions$Run]
conditions <- t(conditions)
conditions <- conditions[2,]
conditions<-factor(conditions)

pvalues <- sapply(1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
  p=wilcox.test(gene~conditions, data)$p.value
  return(p)
})

fdr=p.adjust(pvalues,method = "fdr")

# Calculate the fold-change for each gene
conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))] # basal
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))] # other
foldChanges=log2(rowMeans(dataCon1)/rowMeans(dataCon2))

# Output results based on FDR threshold
outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)
write.csv(outRst,"DEG_wilcox_BASAL_vs_OTHER.csv")

# Fetch upregulated genes
upregulated_genes <- outRst[outRst$log2foldChange > 1 & outRst$pValues < 0.05,]
write.csv(upregulated_genes,"upreg_deg_BASAL_vs_OTHER.csv")

# Fetch downregulated genes
downregulated_genes <- outRst[outRst$log2foldChange < -1 & outRst$pValues < 0.05,]
write.csv(downregulated_genes,"downreg_deg_BASAL_vs_OTHER.csv")

####-----

# 2. find DRP between BASAL vs NON_BASAL (OTHER) using Absolute and Relative Promoter Activity:

# Read the result2 file obtained by running proactiv
result2 <- readRDS("result2_35basal_vs_5others.rds")

# part 1: we will first use Absolute promoter activity
abs_pa <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)  
colnames(abs_pa) <- gsub("_SJ.out","",colnames(abs_pa))

# we have to do filtering now so let us fetch the rowData from result2:
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

####----

# Run the Wilcoxon rank-sum test for each gene
count_norm <- as.data.frame(abs_pa_filtered)

# Now make the conditions column
conditions <- metadata[,c("Run","intrinsic_subtype")]
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
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))] # basal
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))] # other
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

####----

# part 2: we will now use Relative promoter activity
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

####----

# Run the Wilcoxon rank-sum test for each gene
count_norm <- as.data.frame(rel_pa_filtered2)

# Define conditions
conditions <- metadata[,c("Run","intrinsic_subtype")]
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
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))] # basal
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))] # other
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

##############################################################################################################

#-------------------------------------------- FETCH THE AAPS -----------------------------------------------##

# Read the result2 file obtained by running proactiv on basal vs others
result2 <- readRDS("result2_35basal_vs_5others.rds")

# Remove internal promoters & keep multi-promoter genes only 
result2_promoter = rowData(result2) 
filtered_df <- as.data.frame(result2_promoter)
filtered_df <- filtered_df[, c('promoterId','geneId','internalPromoter')]

# Filter 1: keep internal promoters = FALSE ------------------------------
filtered_df <- filtered_df[filtered_df$internalPromoter == "FALSE",]

# Filter 2: keep the genes which have >=2 prom ---------------------------
filtered_df <- filtered_df %>%
  group_by(geneId) %>%
  filter(n() >= 2) %>%
  mutate(ID = paste('pr', promoterId, '_', geneId, sep='')) 

filtered_df2 <- filtered_df[,c('promoterId', 'geneId')] # subset required columns

# Read the genes wilcox result file
genes <- as.data.frame(read_csv("DEG_wilcox_BASAL_vs_OTHER.csv"))
colnames(genes)[1] <- 'geneid'

# Fetch the nonsignificant genes
nonsig_genes <- genes[genes$pValues>0.05,]

# Keep them in the above dataframe:
filtered_df2_nonsig_g <- filtered_df2[filtered_df2$geneId%in%nonsig_genes$geneid,]

####

# Read the absolute promoter activity wilcox result file
abspa <- as.data.frame(read_csv('DRP_wilcox_AbsPA.csv'))
colnames(abspa)[1] <- 'promid'
signif_abspa <- abspa[abspa$pValues<0.05,]  # fetch significant DRPs

# Read the relative promoter activity wilcox result file
relpa <- as.data.frame(read_csv('DRP_wilcox_RelPA.csv'))
colnames(relpa)[1] <- 'promid'
signif_relpa <- relpa[relpa$pValues<0.05,]  # fetch significant DRPs

# Intersect to get the DRPs significant by both absolute and relative promoter activity
signif_prom_all <- intersect(signif_abspa$promid,signif_relpa$promid)

# Keep them in the above gene dataframe:
filtered_df2_nonsig_g_sig_p <- filtered_df2_nonsig_g[filtered_df2_nonsig_g$promoterId%in%signif_prom_all,]

################################################################################################

## VALIDATE BASAL-specific AAPs:---->>

intersect(val$promoterId,fuscc$promoterId)
intersect(val$geneId,fuscc$geneId)
com <- fuscc[fuscc$promoterId%in%val$promoterId,]

################################################################################################
################################################################################################
