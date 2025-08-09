#### NOTE: This script contains the code used for running the ProActiv pipeline on BL1 vs REST of the TNBC subtypes in FUSCC & EXTERNAL DATASET 1 ####

# Input files: Metafile, SJ.out files & gtf file 

# Load required libraries:
library(proActiv)
library(dplyr)
library(tidyverse)
library(readxl)

##--------------------For FUSCC dataset------------------->>

# BL1 vs REST (TNBC subtypes)

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
rownames(metadata) <- metadata$Run

# Take all tumor (TNBC) samples only:
metadata <- metadata[metadata$Tumor_Normal1=='Tumor',]
metadata <- metadata[,c('Run','batch','Tumor_Normal1')]

# Read the subtyping result file:
subtype_meta <- as.data.frame(read_csv('TNBC_TYPING_360_FINAL.csv'))
subtype_meta_sel <- subtype_meta[,c('Original_name','tnbctype_4')]
colnames(subtype_meta_sel)[1] <-'Run'
# Remove UNS samples:
remove <- 'UNS'
subtype_meta_sel <- subtype_meta_sel[!subtype_meta_sel$tnbctype_4 %in% remove,]

# Merge the two files:
merged_meta <- merge(metadata, subtype_meta_sel, by='Run')

# Make the bl1 vs rest column here:
merged_meta <- merged_meta %>% mutate(bl1_vs_rest = ifelse(tnbctype_4 == 'BL1', 'BL1','REST'))
write.csv(merged_meta,'fuscc_merged_meta.csv')
rownames(merged_meta) <-merged_meta$Run

# list the STAR junction files as input
files <- list.files(path = './all_SJ_files', full.names = TRUE,pattern = "_SJ.out.tab")  # 448 files
files2=as.data.frame(files)
files2$files=gsub("./all_SJ_files/","",files2$files)
files2$files=gsub("_SJ.out.tab","",files2$files)
files2$original = files

# Check if all samples are intersecting
files2_new = files2[files2$files %in% merged_meta$Run,]

## Prepare promoter annotation
gtf.file <-('gencode.v44.basic.annotation.gtf')
promoterAnnotation.gencode.v44.subset <- preparePromoterAnnotation(file = gtf.file,
                                                                   species = 'Homo_sapiens')

# Before running ProActiv make sure to match the order of samples in metafile and files2 
merged_meta = merged_meta[files2_new$files,]
allfilenames = files2_new$original  # ASSIGNING ALL THE ORIGINAL NAMES TO A NEW VARIABLE

# Define the conditions
conditions <- merged_meta$bl1_vs_rest
promoterAnnotation <- promoterAnnotation.gencode.v44.subset

## Run ProActiv
result <- proActiv(files = allfilenames, promoterAnnotation = promoterAnnotation, condition = conditions)
result_promoter <- rowData(result)
result_promoter = as.data.frame(result_promoter)

# For cleaner downstream analysis remove single-exon transcripts / promoters by eliminating promoter counts that are NA
result2 <- result[complete.cases(assays(result)$promoterCounts),]
result2_promoter = rowData(result2)
result2_promoter = as.data.frame(result2_promoter)
saveRDS(result2,"result2_updated_subype.rds")

############################################################################

##--------------------For EXTERNAL DATASET 1------------------->>

# 2. BL1 vs REST (TNBC subtypes)

# Load the metafile  
metadata <- as.data.frame(read_csv("meta_sel_61samples.csv"))
metadata <- tibble::column_to_rownames(metadata,var = "...1")

# Keeping only the TNBC samples:
metadata <- metadata[metadata$tumor_normal=="TNBC",]

# Read the subtypes information file:
subtype_meta <- read_csv("subtyping_4_result.csv")
subtype_meta_sel <- subtype_meta[,c('Original_name','tnbctype_4')]  # just keep required columns
colnames(subtype_meta_sel)[1]<-'Run'
# Remove UNS samples:
remove <- 'UNS'
subtype_meta_sel <- subtype_meta_sel[!subtype_meta_sel$tnbctype_4 %in% remove,]

# Merge the two files:
merged_meta <- merge(metadata, subtype_meta_sel, by='Run')
rownames(merged_meta)<- merged_meta$Run

# Make the BL1 vs REST column here:
merged_meta <- merged_meta %>% mutate(bl1_vs_rest = ifelse(tnbctype_4 == 'BL1', 'BL1','REST'))
write.csv(merged_meta,'val_merged_meta.csv')
rownames(merged_meta) <-merged_meta$Run

# list the STAR junction files as input
files <- list.files(path = './SJ_files_new', full.names = TRUE,pattern = "_SJ.out.tab")  # 61 files
files2=as.data.frame(files)
files2$files=gsub("./SJ_files_new/","",files2$files)
files2$files=gsub("_SJ.out.tab","",files2$files)
files2$original = files

# Check if all samples are intersecting
files2_new = files2[files2$files %in% merged_meta$Run,]

## Prepare promoter annotation
gtf.file <-('gencode.v44.basic.annotation.gtf')
promoterAnnotation.gencode.v44.subset <- preparePromoterAnnotation(file = gtf.file,
                                                                   species = 'Homo_sapiens')

# Before running ProActiv make sure to match the order of samples in metafile and files2 
merged_meta = merged_meta[files2_new$files,]
allfilenames = files2_new$original  # ASSIGNING ALL THE ORIGINAL NAMES TO A NEW VARIABLE

# Define the conditions
conditions <- merged_meta$bl1_vs_rest
promoterAnnotation <- promoterAnnotation.gencode.v44.subset

## Run ProActiv
result <- proActiv(files = allfilenames, promoterAnnotation = promoterAnnotation, condition = conditions)
result_promoter <- rowData(result)
result_promoter = as.data.frame(result_promoter)

# For cleaner downstream analysis remove single-exon transcripts / promoters by eliminating promoter counts that are NA
result2 <- result[complete.cases(assays(result)$promoterCounts),]
result2_promoter = rowData(result2)
result2_promoter = as.data.frame(result2_promoter)
saveRDS(result2,"result2_updated_subype_val.rds")

##########################################################################################################################
##########################################################################################################################
