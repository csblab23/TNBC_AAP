#### NOTE: This script contains the code used for running the ProActiv pipeline on TNBC and ADJACENT NORMALS in FUSCC & EXTERNAL DATASET 1 ####

# Input files: Metafile, SJ.out files & gtf file 

# Load required libraries:
library(proActiv)
library(dplyr)
library(tidyverse)
library(readxl)

##--------------------For FUSCC dataset------------------->>

# TNBC vs ADJACENT NORMAL

# Load the metafile:
metadata <- as.data.frame(read_csv("meta_with_batchinfo.csv"))
metadata=metadata[,-1]  # removing the extra s.no. column
metadata$Tumor_Normal1=metadata$Tumor_Normal

# cleaning the tissue type column
metadata =metadata %>% 
  mutate(Tumor_Normal1 = case_when(Tumor_Normal=="TumorTissue"  ~ "Tumor",
                                   Tumor_Normal=="TumorTissue.2018"  ~ "Tumor",
                                   Tumor_Normal=="TumorTissue.rep"  ~ "Tumor",
                                   Tumor_Normal=="PairedNormalTissue"  ~ "AdjNormal",
                                   Tumor_Normal=="PairedNormalTissue.rep"  ~ "AdjNormal"))
rownames(metadata) <- metadata$Run

# Prepare the input data for running proActiv
# list the STAR junction files as input
files <- list.files(path = './all_SJ_files/', full.names = TRUE)  # 448 files
files2=as.data.frame(files)
files2$files=gsub("./all_SJ_files//","",files2$files)
files2$files=gsub("_SJ.out.tab","",files2$files)
files2$original = files

# Prepare promoter annotation:
gtf.file <-('gencode.v44.basic.annotation.gtf')
promoterAnnotation.gencode.v44.subset <- preparePromoterAnnotation(file = gtf.file,
                                                                   species = 'Homo_sapiens')

# Match the order of samples in metafile and SJ.out.tab files 
metadata = metadata[files2$files,]
allfilenames = files2$original  # ASSIGNING ALL THE ORIGINAL NAMES TO A NEW VARIABLE

# Define the conditions of interest
conditions <- metadata$Tumor_Normal1
promoterAnnotation <- promoterAnnotation.gencode.v44.subset

## Run ProActiv:
result <- proActiv(files = allfilenames, promoterAnnotation = promoterAnnotation, condition = conditions)
result_promoter <- rowData(result)
result_promoter = as.data.frame(result_promoter)

# For cleaner downstream analysis remove single-exon transcripts / promoters by eliminating promoter counts that are NA
result2 <- result[complete.cases(assays(result)$promoterCounts),]
result2_promoter = rowData(result2)
result2_promoter = as.data.frame(result2_promoter)
saveRDS(result2,'result2_tum_vs_adjnormal.rds')

###########################################################

##--------------------For EXTERNAL DATASET 1------------------->>

# TNBC vs ADJACENT NORMAL

# Read the metadata:
meta <- as.data.frame(read_csv("metadata_subset.csv"))
meta <- tibble::column_to_rownames(meta, '...1')
# we had removed 2 TNBC samples due to poor quality reads and rsem issue due to files getting truncated
remove <- c('SRR1313139', 'SRR1313140')
meta_sel <- meta[!meta$Run %in% remove,]
write.csv(meta_sel,"meta_sel_61samples.csv")

####----

# list the STAR junction files as input
files <- list.files(path = './SJ_files', full.names = TRUE)  # 63 files
files2=as.data.frame(files)
files2$files=gsub("./SJ_files/","",files2$files)
files2$files=gsub("_SJ.out.tab","",files2$files)
files2$original = files

# Check if all samples are intersecting
files2_new = files2[files2$files %in% meta_sel$Run,]

## Prepare promoter annotation
gtf.file <-('gencode.v44.basic.annotation.gtf')
promoterAnnotation.gencode.v44.subset <- preparePromoterAnnotation(file = gtf.file,
                                                                   species = 'Homo_sapiens')

# Before running ProActiv make sure to match the order of samples in metafile and files2 
meta_sel = meta_sel[files2_new$files,]
allfilenames = files2_new$original  # ASSIGNING ALL THE ORIGINAL NAMES TO A NEW VARIABLE

# Define the conditions
conditions <- meta_sel$tumor_normal
promoterAnnotation <- promoterAnnotation.gencode.v44.subset

## Run ProActiv
result <- proActiv(files = allfilenames, promoterAnnotation = promoterAnnotation, condition = conditions)
result_promoter <- rowData(result)
result_promoter = as.data.frame(result_promoter)

# For cleaner downstream analysis remove single-exon transcripts / promoters by eliminating promoter counts that are NA
result2 <- result[complete.cases(assays(result)$promoterCounts),]
result2_promoter = rowData(result2)
result2_promoter = as.data.frame(result2_promoter)
saveRDS(result2,"result2_40tnbc_vs_21normal.rds")

##########################################################################################################################
##########################################################################################################################
