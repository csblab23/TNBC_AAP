############## This script includes all the codes used for obtaining Differentially Expressed Genes (DEGs) & Differentially Regulated Promoters (DRPs) in TNBC ################

# Load required libraries
library(edgeR)
library(dplyr)
library(tidyverse)

# Find differentially expressed genes (DEG) between TNBC vs Normal in FUSCC dataset 
# Load the metafile  
metadata <- as.data.frame(read_csv("meta_with_batchinfo.csv"))
metadata <- metadata[,-1]
metadata$Tumor_Normal1=metadata$Tumor_Normal
metadata =metadata %>% 
  mutate(Tumor_Normal1 = case_when(Tumor_Normal=="TumorTissue"  ~ "Tumor",
                                   Tumor_Normal=="TumorTissue.2018"  ~ "Tumor",
                                   Tumor_Normal=="TumorTissue.rep"  ~ "Tumor",
                                   Tumor_Normal=="PairedNormalTissue"  ~ "AdjNormal",
                                   Tumor_Normal=="PairedNormalTissue.rep"  ~ "AdjNormal"))
rownames(metadata) <- metadata$Run  # making the sample names as rownames
subset_meta <- metadata[,c('Run', 'batch','Tumor_Normal1')]   # just keep the required columns

# Read the feature-counts matrix 
count <- as.data.frame(read_csv("Featurecounts_combined_fuscc.csv"))
count=count[,-1] # remove extra s.no. column
count <- tibble::column_to_rownames(count, var = "Geneid")

# Now calculate cpm  
cpm <- cpm(count)
is.exprs <- rowSums(cpm>1) >= 20 # Apply filtering criteria
counts2 <- count[is.exprs, ]

# Make samples order same in both the files 
com_exp_mat=counts2[,subset_meta$Run]

# Prepare DGEList and voom normalize the matrix
x <- DGEList(counts=com_exp_mat)
x <- calcNormFactors(x, method='TMM')
v <- voom(x, plot=T)
vMat <- v$E

# Adjust batch by sva combat here
subset_meta$batch=as.factor(subset_meta$batch)
exprMat <- sva::ComBat(vMat, batch=subset_meta$batch)

####-----

# Run the Wilcoxon rank-sum test for each gene
count_norm <- as.data.frame(exprMat)

# Define conditions 
conditions <- subset_meta[,c("Run","Tumor_Normal1")]
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
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]  # AdjNormal
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]  # Tumor
foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))

# Output results based on FDR threshold
outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)
write.csv(outRst,"DEG_wilcoxon_tum_vs_nor.csv")

# Fetch upregulated genes
upregulated_genes <- outRst[outRst$log2foldChange > 1 & outRst$pValues < 0.05,]
write.csv(upregulated_genes,"upreg_2443deg_bypvalue.csv")

# Fetch downregulated genes
downregulated_genes <- outRst[outRst$log2foldChange < -1 & outRst$pValues < 0.05,]
write.csv(downregulated_genes,"down_2039deg_bypvalue.csv")

#############################################################################

## Find differentially regulated promoters (DRP) by absolute promoter activity 
# Filter applied on Absolute promoter activity: Removing promoters with median activity < 0.25 

# Load the metafile  
metadata <- as.data.frame(read_csv("meta_with_batchinfo.csv"))
metadata <- metadata[,-1]
metadata$Tumor_Normal1=metadata$Tumor_Normal
metadata =metadata %>% 
  mutate(Tumor_Normal1 = case_when(Tumor_Normal=="TumorTissue"  ~ "Tumor",
                                   Tumor_Normal=="TumorTissue.2018"  ~ "Tumor",
                                   Tumor_Normal=="TumorTissue.rep"  ~ "Tumor",
                                   Tumor_Normal=="PairedNormalTissue"  ~ "AdjNormal",
                                   Tumor_Normal=="PairedNormalTissue.rep"  ~ "AdjNormal"))
rownames(metadata) <- metadata$Run  # making the sample names as rownames
subset_meta <- metadata[,c('Run', 'batch','Tumor_Normal1')]   # let us just keep the required columns

# Read the result2 file:
result2 <- readRDS("result2_tum_vs_adjnormal.rds")

# Fetch Absolute promoter activity
abs_pa <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)  
colnames(abs_pa) <- gsub("_SJ.out","",colnames(abs_pa))

# Filtering
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

# Now keep the same promoter ids in the absolute counts matrix
abs_pa_filtered <- abs_pa[rownames(abs_pa) %in% filtered_df$promoterId,]

# Checking the median activity of the promoters to apply filter in next step
row_summaries <- t(apply(abs_pa_filtered, 1, summary))
row_summaries <- as.data.frame(row_summaries)

# Filter promoters with median activity < 0.25
abs_pa_filtered_2 <- abs_pa_filtered[apply(abs_pa_filtered, 1, median) > 0.25, ]

####----

# Run the Wilcoxon rank-sum test for each gene
count_norm <- as.data.frame(abs_pa_filtered_2)

# Define conditions 
conditions <- subset_meta[,c("Run","Tumor_Normal1")]
count_norm <- count_norm[,conditions$Run]
conditions <- t(conditions)
conditions <- conditions[2,]
conditions<-factor(conditions)

# Make sure to have the order of samples same:
count_norm <- count_norm[,subset_meta$Run]

pvalues <- sapply(1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
  p=wilcox.test(gene~conditions, data)$p.value
  return(p)
})

fdr=p.adjust(pvalues,method = "fdr")

# Calculate the fold-change for each gene
conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))] # adjnormal
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))] # tumor
foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))

# Output results 
outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)
write.csv(outRst,"DRP_wilcox_AbsPA_filtered.csv")

# Fetch upregulated DRPs
upregulated_genes <- outRst[outRst$log2foldChange > 1.2 & outRst$pValues < 0.05,]
write.csv(upregulated_genes,"upreg_drp_absPA_bypvalue.csv")

# Fetch downregulated DRPs
downregulated_genes <- outRst[outRst$log2foldChange < -1.2 & outRst$pValues < 0.05,]
write.csv(downregulated_genes,"down_drp_absPA_bypvalue.csv")

#####################################################################
#####################################################################
