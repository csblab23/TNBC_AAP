#### NOTE: This script contains the code used for running the ProActiv pipeline on BASAL and NON BASAL samples in FUSCC & EXTERNAL DATASET 1 ####

# Input files: Metafile, SJ.out files & gtf file 

# Load required libraries:
library(proActiv)
library(dplyr)
library(tidyverse)
library(readxl)

##--------------------For FUSCC dataset------------------->>

# BASAL vs REST (using PAM50 Intrinsic subtype information from FUSCC metafile)

# Load the metafile  
metadata <- as.data.frame(read_csv("meta_with_batchinfo.csv"))
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

## Prepare the input data
# list the STAR junction files as input
files <- list.files(path = './all_SJ_files/', full.names = TRUE)  # 448 files
files2=as.data.frame(files)
files2$files=gsub("./all_SJ_files//","",files2$files)
files2$files=gsub("_SJ.out.tab","",files2$files)
files2$original = files

## Prepare promoter annotation
gtf.file <-('gencode.v44.basic.annotation.gtf')
promoterAnnotation.gencode.v44.subset <- preparePromoterAnnotation(file = gtf.file,
                                                                   species = 'Homo_sapiens')

# Before running ProActiv make sure to match the order of samples in metafile and files2 
files2 <- files2[files2$files%in%metadata$Run,]
metadata = metadata[files2$files,]
allfilenames = files2$original  # ASSIGNING ALL THE ORIGINAL NAMES TO A NEW VARIABLE

# Define the conditions
conditions <- metadata$Intrinsic_Subtype
promoterAnnotation <- promoterAnnotation.gencode.v44.subset

## Run ProActiv
result <- proActiv(files = allfilenames, promoterAnnotation = promoterAnnotation, condition = conditions)
result_promoter <- rowData(result)
result_promoter = as.data.frame(result_promoter)

# For cleaner downstream analysis remove single-exon transcripts / promoters by eliminating promoter counts that are NA
result2 <- result[complete.cases(assays(result)$promoterCounts),]
result2_promoter = rowData(result2)
result2_promoter = as.data.frame(result2_promoter)
saveRDS(result2,'result2_intrinsic_subtype.rds')

############################################################################

##--------------------For EXTERNAL DATASET 1------------------->>

# BASAL vs REST

# Read the metadata:
meta <- as.data.frame(read_csv("All_tumor_metadata.csv"))
meta <- meta[,-1]
rownames(meta) <- meta$Run

# just keep the TNBC samples:
meta_tnbc <- meta[meta$Tumor=='TNBC',]

# Add the PAM50subtype information:
pam50 <- as.data.frame(read_csv('TNBC_PAM50results.csv'))
pam50 <- pam50[,-1]
colnames(pam50)[1] <- 'Run'
pam50 <- pam50[,1:2]

## merge the two dfs:
merged_meta <- merge(meta_tnbc, pam50, by='Run')
rownames(merged_meta)<-merged_meta$Run
merged_meta <- merged_meta %>%
  mutate(intrinsic_subtype = ifelse(subtype == 'Basal', 'Basal', 'Other'))
write.csv(merged_meta,"TNBC_mergedmeta.csv")

# list the STAR junction files as input
files <- list.files(path = './SJ_files', full.names = TRUE,pattern = "_SJ.out.tab")  # 61 files
files2=as.data.frame(files)
files2$files=gsub("./SJ_files/","",files2$files)
files2$files=gsub("_SJ.out.tab","",files2$files)
files2$original = files
files2_new = files2[files2$files %in% merged_meta$Run,]

## Prepare promoter annotation----->>
gtf.file <-('gencode.v44.basic.annotation.gtf')
promoterAnnotation.gencode.v44.subset <- preparePromoterAnnotation(file = gtf.file,
                                                                   species = 'Homo_sapiens')

# Before running ProActiv make sure to match the order of samples in metafile and files2 
merged_meta = merged_meta[files2_new$files,]
allfilenames = files2_new$original  # ASSIGNING ALL THE ORIGINAL NAMES TO A NEW VARIABLE

# Define the conditions
conditions <- merged_meta$intrinsic_subtype
promoterAnnotation <- promoterAnnotation.gencode.v44.subset

## Run ProActiv
result <- proActiv(files = allfilenames, promoterAnnotation = promoterAnnotation, condition = conditions)
result_promoter <- rowData(result)
result_promoter = as.data.frame(result_promoter)

# For cleaner downstream analysis remove single-exon transcripts / promoters by eliminating promoter counts that are NA
result2 <- result[complete.cases(assays(result)$promoterCounts),]
result2_promoter = rowData(result2)
result2_promoter = as.data.frame(result2_promoter)
saveRDS(result2,"result2_35basal_vs_5others.rds")
#####################################################################################################
#####################################################################################################
