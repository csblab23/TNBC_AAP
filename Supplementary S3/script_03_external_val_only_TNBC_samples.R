## Code to make Supplementary S3d,e,f,g ## 
# This includes PROGNOSTIC RESULTS VALIDATION in GSE240671 cohort --->> only TNBC samples taken and proActiv ran using condition = NULL

# Load required libraries
library(survival)
library(survminer)
library(proActiv)
library(dplyr)
library(tidyverse)
library(readxl)
library(data.table)
library(ggplot2)
library(reshape2)
library(EnvStats)
library(ggpubr)
library(edgeR)
library(limma)

# Load the metafile
metadata <- as.data.frame(fread("./META_info_all_samples.txt"))

# Taking only the pre-treatment samples
metadata <- metadata[metadata$timing_biopsies == 'pre-surgery',]

# Keeping only the TNBC samples
metadata <- metadata[metadata$molecular_category=='TNBC',]
rownames(metadata) <- metadata$Run

## Prepare the input data
# list the STAR junction files as input
files <- list.files(path = './prog_ext_vali_GSE240671/all_files/', full.names = TRUE)  # 27 files
files2=as.data.frame(files)
files2$files=gsub("./prog_ext_vali_GSE240671/all_files//","",files2$files)
files2$files=gsub("_SJ.out.tab","",files2$files)
files2$original = files

## Prepare promoter annotation
gtf.file <-('gencode.v44.basic.annotation.gtf')
promoterAnnotation.gencode.v44.subset <- preparePromoterAnnotation(file = gtf.file,
                                                                   species = 'Homo_sapiens')

# Before running ProActiv make sure to match the order of samples in metafile and files2 
metadata = metadata[files2$files,]
allfilenames = files2$original  # ASSIGNING ALL THE ORIGINAL NAMES TO A NEW VARIABLE

# Define the conditions  --->> No conditions here as all are TNBC samples

promoterAnnotation <- promoterAnnotation.gencode.v44.subset

## Run ProActiv
result <- proActiv(files = allfilenames, promoterAnnotation = promoterAnnotation, condition = NULL)
rowData(result)
result_promoter <- rowData(result)
result_promoter = as.data.frame(result_promoter)

# For cleaner downstream analysis remove single-exon transcripts / promoters by eliminating promoter counts that are NA
result2 <- result[complete.cases(assays(result)$promoterCounts),]
result2_promoter = rowData(result2)
result2_promoter = as.data.frame(result2_promoter)
saveRDS(object = result2, file = 'result_GSE240671_TNBC_27_SAMPLE.rds')

####----

## Visualise the expression of promoters in high and low risk group patients ----->>
# Here:
# High risk group patients = Residual disease --- moderate or extensive
# Low risk group patients = residual disease --- minimal or no

# Step-01 -->> defining high and low risk group patients
metadata$risk_group <- ifelse(metadata$rcb_category %in% c('moderate residual disease', 'extensive residual disease') , "High", "Low")

# Step-02 -->> take expression of pr41492 in those patients --- and then make the boxplot of it
result2 <- readRDS('result_GSE240671_TNBC_27_SAMPLE.rds')
ab_counts <- as.data.frame(result2@assays@data@listData[["absolutePromoterActivity"]])
colnames(ab_counts) <- gsub("_SJ.out", "", colnames(ab_counts))
ab_counts <- ab_counts[c( '3003','3004',  '41492', '41490', '41491'), ]
rownames(ab_counts) <- paste0("pr", rownames(ab_counts))
ab_counts <- as.data.frame(t(ab_counts))


# Merge
df <- merge(metadata, ab_counts, by = 0)
rownames(df) <- df$Row.names
df <- df[,c(50:55)]

# Step-03 -->> now make boxplot of expression all the promoters in high and low category
df$risk_group <- as.factor(df$risk_group)

# Boxplot for FTX promoters------------------------------------------->>
subsetdf <- df[,c(1,3,4)]

# Reshape the data from wide to long format using melt from reshape2
data_long <- melt(subsetdf, id.vars = "risk_group", 
                  variable.name = "gene", 
                  value.name = "expression")

# Plot
p <- ggboxplot(data_long, x = "risk_group", y = "expression", 
               color = "risk_group", fill = "risk_group" ,size = 1.5,
               facet.by = "gene",notch = FALSE,add="jitter") +
  scale_color_manual(values = c( "darkorchid3","#298c8c")) +  # Custom color palette
  scale_fill_manual(values = c("thistle2","#b3d7d7")) +
  labs(x = "Risk Group", 
       y = "Absolute Promoter activity") + 
  theme(
    axis.text = element_text(size = 18, color = "black"),  # Adjust axis text size
    axis.title = element_text(size = 18,face = "bold"),  # Adjust axis title size
    strip.text = element_text(size = 16,face = "bold"),  # Facet label size
    plot.title = element_text(hjust = 0.5, size = 16),  # Centered title
    legend.position = "none"  # Remove legend
  ) + scale_y_continuous(limits = c(-0.6, 7))

p <- p + stat_compare_means(aes(group = risk_group),
                            size=6,label.y = 6.5) + stat_n_text(size = 7,
                                                                color = "black",
                                                                fontface = "bold") 
p
ggsave("tnbc_FTXplot.png", p, width = 3.2, height = 5.3, dpi = 1000)

####----

# Boxplot for HUWE1 promoters------------------------------------------------------>>
subsetdf <- df[,1:2]

# Reshape the data from wide to long format using melt from reshape2
data_long <- melt(subsetdf, id.vars = "risk_group", 
                  variable.name = "gene", 
                  value.name = "expression")

# Plot
p <- ggboxplot(data_long, x = "risk_group", y = "expression", 
               color = "risk_group", fill = "risk_group" ,size = 1.5,
               facet.by = "gene",notch = FALSE,add="jitter") +
  scale_color_manual(values = c( "darkorchid3","#298c8c")) +  # Custom color palette
  scale_fill_manual(values = c("thistle2","#b3d7d7")) +
  labs(x = "Risk Group", 
       y = "Absolute Promoter activity") + 
  theme(
    axis.text = element_text(size = 18, color = "black"),  # Adjust axis text size
    axis.title = element_text(size = 18,face = "bold"),  # Adjust axis title size
    strip.text = element_text(size = 16,face = "bold"),  # Facet label size
    plot.title = element_text(hjust = 0.5, size = 16),  # Centered title
    legend.position = "none"  # Remove legend
  ) + scale_y_continuous(limits = c(-0.6, 7))

p <- p + stat_compare_means(aes(group = risk_group),
                            size=6,label.y = 6.5) + stat_n_text(size = 7,
                                                                color = "black",
                                                                fontface = "bold") 
p
ggsave("tnbc_HUWE1plot.png", p, width = 3.2, height = 5.3, dpi = 1000)

##################################################################################################

## FOR GENE EXPRESSION: 

# Checking what is the expression of gene in the high and low category ->
gene_exp <- as.data.frame(fread('output.txt'))
colnames(gene_exp) <- gene_exp[1,]
gene_exp <- gene_exp[-1,]
colnames(gene_exp) <- gsub("GSE240671/rna_seq_pipeline/star/", "",colnames(gene_exp))
colnames(gene_exp) <- gsub("_Aligned.sortedByCoord.out.bam", "",colnames(gene_exp))
rownames(gene_exp) <- gene_exp$Geneid
gene_exp <- gene_exp[,-1]
gene_exp[] <- lapply(gene_exp, function(x) as.numeric(as.character(x)))

x <- DGEList(counts=gene_exp )
x <- calcNormFactors(x,method = "TMM")
v <- voom(x, plot=F)
vMat <- v$E

# Batch adjustment
metadata$Sequencing_batch <- as.factor(metadata$Sequencing_batch)
vMat_2 <- sva::ComBat(vMat, batch=metadata$Sequencing_batch)
vMat_2 <- as.data.frame(vMat_2)

vMat_2 <- vMat_2[c('ENSG00000086758.17','ENSG00000230590.13'), ]  # take genes of interest
vMat_2 <- as.data.frame(t(vMat_2))
df <- merge(metadata, vMat_2, by = 0)
rownames(df) <- df$Row.names
df <- df[,c(50:52)]
df$risk_group <- as.factor(df$risk_group)

# Reshape the data from wide to long format using melt from reshape2
data_long <- melt(df, id.vars = "risk_group", 
                  variable.name = "gene", 
                  value.name = "expression")

p <- ggboxplot(data_long, x = "risk_group", y = "expression", 
               color = "risk_group", size = 1.5,
               facet.by = "gene",notch = FALSE,add="jitter") +
  labs(x = "Risk Group", 
       y = "Normalized gene expression") + 
  scale_color_manual(values = c( "#b67eff","aquamarine3")) + 
  theme(
    text = element_text(size = 25),
    axis.title = element_text(size = 18,face = "bold"),
    axis.text = element_text(size = 18, color = "black"),
    strip.text = element_blank(),
    legend.position = "none"
  ) + scale_y_continuous(limits = c(6.9, 12),expand = c(0,0))

p <- p + stat_compare_means(aes(group = risk_group),
                            size=6,label.y = 11.4) + stat_n_text(size = 7,
                                                                 color = "black",
                                                                 fontface = "bold") 
p
ggsave("tnbc_geneboxplot.png", p, width = 5.5, height = 5, dpi = 1000)

###################################################################################
###################################################################################
