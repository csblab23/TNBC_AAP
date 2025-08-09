############## This script includes all the codes used for SUPPLEMENTARY S3 ################

# Load required libraries
library(dplyr)
library(tidyverse)
library(FactoMineR)
library(factoextra)


## Code For Supplementary figure S3a & S3c:---------------------->>

set.seed(456)

## Making PCA plots for genes:-->
# Load the FUSCC metafile
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
subset_meta <- metadata[,c('Run', 'batch','Tumor_Normal1','Intrinsic_Subtype')]   # let us just keep the required columns
subset_meta <- subset_meta[subset_meta$Tumor_Normal1=='Tumor',]

# Read the feature-counts matrix 
count <- as.data.frame(read_csv("Featurecounts_combined.csv"))
count=count[,-1] # remove extra s.no. column
count <- tibble::column_to_rownames(count, var = "Geneid")

# Match the sample number and order:
count <- count[,colnames(count) %in% subset_meta$Run,]
count <- count[,subset_meta$Run]

# Now calculate cpm  
cpm <- cpm(count)
is.exprs <- rowSums(cpm>1) >= 20 # Apply filtering criteria
counts2 <- count[is.exprs, ]

# Make sure sample order is same in both the files 
com_exp_mat=counts2[,subset_meta$Run]

# Prepare DGEList and voom normalize the matrix
x <- DGEList(counts=com_exp_mat)
x <- calcNormFactors(x, method='TMM')
v <- voom(x, plot=T)
vMat <- v$E

# Adjust batch by sva combat here
subset_meta$batch=as.factor(subset_meta$batch)
exprMat <- sva::ComBat(vMat, batch=subset_meta$batch)

# List of AAPs with names and promoter IDs
aap_lists <- list(
  fuscc = list(
    file = "nonsig_g_sig_p_fuscc.csv",
    id_col = "geneId",
    keep = NULL
  ),
  com9 = list(
    file = "common9aaps.csv",
    id_col = "promoterId",
    keep = c(10005,7033,2055,4814,30507,12842)
  )
)

for (aap_name in names(aap_lists)) {
  cat("\nProcessing:", aap_name, "\n")
  
  # Read data
  aap_df <- read_csv(aap_lists[[aap_name]]$file)
  
  # Filter for specific promoters if keep is set
  if (!is.null(aap_lists[[aap_name]]$keep)) {
    aap_df <- aap_df %>% filter(!!sym(aap_lists[[aap_name]]$id_col) %in% aap_lists[[aap_name]]$keep)
  }
  
  # Select rows in exprMat matching the gene IDs/promoter IDs
  gene_ids <- aap_df[[aap_lists[[aap_name]]$id_col]]
  exprMat_sub <- exprMat[rownames(exprMat) %in% gene_ids, ]
  
  # Transpose and merge with metadata
  t_expr <- t(exprMat_sub)
  df <- merge(subset_meta, t_expr, by = 0)
  df <- tibble::column_to_rownames(df, "Row.names")
  
  # Clean and factor subtype
  df$Intrinsic_Subtype <- gsub("Other", "Non_Basal", df$Intrinsic_Subtype)
  df$Intrinsic_Subtype <- factor(df$Intrinsic_Subtype, levels = c("Basal", "Non_Basal"))
  
  # PCA
  res.pca <- PCA(df[, 5:ncol(df)], graph = FALSE)
  
  # Plot
  p <- fviz_pca_ind(res.pca,
                    geom = "point",
                    pointshape = 19,
                    pointsize = 5,
                    alpha.ind = 0.8,
                    label = "none",
                    title = paste0(aap_name, " PCA"),
                    habillage = df$Intrinsic_Subtype,
                    palette = c("Basal" = "#4ecb8d", "Non_Basal" = "honeydew4"),
                    axes = c(1, 2),
                    addEllipses = FALSE,
                    mean.point = FALSE) +
    theme(axis.title.x = element_text(size = 32, face = "bold"),
          axis.title.y = element_text(size = 32, face = "bold"),
          axis.text.x = element_text(size = 30, colour = "black"),
          axis.text.y = element_text(size = 30, colour = "black"),
          plot.title = element_text(size = 35, hjust = 0.5, face = "bold"),
          legend.text = element_text(size = 32, colour = "black"),
          legend.title = element_blank())
  
  print(p)
  
  # Save plots
  ggsave(paste0(aap_name, "_pca_plot.png"), p, width = 8, height = 6, dpi = 1000)
}

#############################################################################
#############################################################################

## Code for Supplementary figure S3b & S3d:------------------->>

# PCA plots for promoters

# Read the result2 file:
result2 <- readRDS("result2_intrinsic_subtype.rds")

# Use Absolute promoter activity
abs_pa <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)  
colnames(abs_pa) <- gsub("_SJ.out","",colnames(abs_pa))
abs_pa <- abs_pa[,subset_meta$Run]

# Define list of AAP sets
aap_sets <- list(
  fuscc = list(
    file = "nonsig_g_sig_p_fuscc.csv",
    id_col = "promoterId",
    keep = NULL  # Keep all rows in the file
  ),
  validated = list(
    file = "common9aaps.csv",
    id_col = "promoterId",
    keep = c(10005,7033,2055,4814,30507,12842)
  )
)

# Loop through each dataset
for (aap_name in names(aap_sets)) {
  cat("\nProcessing:", aap_name, "\n")
  
  # Load the AAP file
  aap_info <- aap_sets[[aap_name]]
  aap_df <- read_csv(aap_info$file)
  
  # Filter if `keep` is specified
  if (!is.null(aap_info$keep)) {
    aap_df <- aap_df %>% filter(!!sym(aap_info$id_col) %in% aap_info$keep)
  }
  
  # Extract promoter IDs
  promoter_ids <- as.character(aap_df[[aap_info$id_col]])
  
  # Filter abs_pa
  abs_pa_filtered <- abs_pa[rownames(abs_pa) %in% promoter_ids, ]
  
  # Transpose and merge
  t_pa <- t(abs_pa_filtered)
  df <- merge(subset_meta, t_pa, by = 0)
  df <- tibble::column_to_rownames(df, "Row.names")
  
  # Format subtype column
  df$Intrinsic_Subtype <- gsub("Other", "Non_Basal", df$Intrinsic_Subtype)
  df$Intrinsic_Subtype <- factor(df$Intrinsic_Subtype, levels = c("Basal", "Non_Basal"))
  
  # Run PCA
  res.pca <- PCA(df[, 5:ncol(df)], graph = FALSE)
  
  # Plot
  p <- fviz_pca_ind(res.pca,
                    geom = "point",
                    pointshape = 19,
                    pointsize = 5,
                    alpha.ind = 0.8,
                    label = "none",
                    title = paste0("AAPs - ", aap_name),
                    habillage = df$Intrinsic_Subtype,
                    palette = c("Basal" = "#4ecb8d", "Non_Basal" = "honeydew4"),
                    axes = c(1, 2),
                    addEllipses = FALSE,
                    mean.point = FALSE) +
    theme(
      axis.title.x = element_text(size = 34, face = "bold"),
      axis.title.y = element_text(size = 34, face = "bold"),
      axis.text.x = element_text(size = 32, colour = "black"),
      axis.text.y = element_text(size = 32, colour = "black"),
      plot.title = element_text(size = 38, hjust = 0.5, face = "bold"),
      legend.text = element_text(size = 32, colour = "black"),
      legend.title = element_blank()
    )
  
  print(p)
  
  # Save plot
  ggsave(paste0("prom_pca_", aap_name, ".png"), p, width = 8, height = 6, dpi = 1000)
}

#############################################################################
#############################################################################

## Code For Supplementary figure S3e & S3f:---------------------->>

# Load the FUSCC metafile with TNBC subtype information
metadata <- as.data.frame(read_csv("fuscc_merged_meta.csv"))
metadata=tibble::column_to_rownames(metadata,var = '...1')
subset_metadata <- metadata

# Read the feature-counts matrix 
count <- as.data.frame(read_csv("Featurecounts_combined.csv"))
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

####-----

# Read the result2 file:
result2 <- readRDS("result2_updated_subype.rds")

# Load Absolute promoter activity
abs_pa <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)  
colnames(abs_pa) <- gsub("_SJ.out","",colnames(abs_pa))

# Match sample order
abs_pa <- abs_pa[,colnames(abs_pa)%in%subset_metadata$Run]
abs_pa <-abs_pa[,subset_metadata$Run]

# Fetch LSP1 promoters:
keep <- c('10530', '10532','10533')
abs_pa_sel <- abs_pa[rownames(abs_pa) %in% keep,]
abs_pa_sel_t <- t(abs_pa_sel)
colnames(abs_pa_sel_t)<-paste0("pr",colnames(abs_pa_sel_t))

# Merge with metadata & make boxplots:---------------------------------->>
finaldf <- merge(subset_metadata,abs_pa_sel_t,by = 0)
finaldf$bl1_vs_rest <- factor(finaldf$bl1_vs_rest, levels = c("BL1","REST"))

# Reshape data here:
finaldf <- finaldf %>%
  pivot_longer(cols = c(pr10530,pr10532,pr10533), 
               names_to = "promoters", 
               values_to = "expression")

# Define colors
my_colors <- c("BL1"= "darkorange1","REST"="orchid3")

# Boxplot
g <- ggboxplot(finaldf, x = "promoters", y = "expression", 
               color = "bl1_vs_rest",size = 1.5) +
  labs(y = "Absolute promoter activity") + 
  scale_color_manual(values = my_colors) + 
  theme(
    text = element_text(size = 25),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, face = "bold"),
    axis.text.x = element_text(size = 24, colour = "black"),
    axis.text.y = element_text(size = 24, colour = "black"),
    legend.text = element_text(size = 25),
    legend.title = element_blank()
  )

g<-g + stat_compare_means(aes(group = bl1_vs_rest), 
                          size=12,label = "p.signif",
                          label.y = 9,
                          label.x.npc = "middle") + scale_y_continuous(limits = c(0,10),expand = c(0,0))
g
ggsave("./fuscc/lsp1_pr.png", g, width = 6, height = 5.5, dpi = 1000)

######------

## For LSP1 gene now:
# LSP1: ENSG00000130592.17

g <- grep("ENSG00000130592.17",rownames(exprMat))
exprMat_sel <- as.data.frame(exprMat[g,])
mm <- merge(subset_metadata, exprMat_sel, by=0)
colnames(mm)[5]<- 'Condition'
colnames(mm)[7]<- 'Expression'

mm$Condition <- factor(mm$Condition , levels = c("BL1","REST"))
my_colors <- c("BL1"= "darkorange1","REST"="orchid3")

g <- ggboxplot(mm, x = "Condition", y = "Expression", 
               color = "Condition", size = 1.5) +
  labs(y = "Normalized gene expression") + 
  scale_color_manual(values = my_colors) + 
  theme(
    text = element_text(size = 25),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 23, face = "bold"),
    axis.text.x = element_text(size = 24, colour = "black"),
    axis.text.y = element_text(size = 24, colour = "black"),
    legend.text = element_text(size = 25),
    legend.title = element_blank()
  )

g<-g + stat_compare_means(aes(group = Condition), 
                          size=12,label = "p.signif",
                          label.y = 7.4,
                          label.x.npc = "middle") + stat_n_text(size = 8,
                                                                color = "black",
                                                                fontface = "bold") +
  scale_y_continuous(limits = c(0.4, 8),expand = c(0,0))

g
ggsave("./fuscc/lsp1_gene.png", g, width = 5.4, height = 5.6, dpi = 1000)

##########################################################################
##########################################################################

## Code For Supplementary figure S3g & S3h:---------------------->>

# Box plots for validation dataset/external dataset 1
# Load the metafile having TNBC subtype information  
metadata <- as.data.frame(read_csv("val_merged_meta.csv"))
metadata <- tibble::column_to_rownames(metadata,var = "...1")

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
is.exprs <- rowSums(cpm>1) >= 20 # Apply filtering criteria
counts2 <- count[is.exprs, ]

# Prepare DGEList and voom normalize the matrix
x <- DGEList(counts=counts2)
x <- calcNormFactors(x, method='TMM')
v <- voom(x, plot=T)
vMat <- v$E

####----

# Read the result2 file obtained by running proactiv on tnbc type4 (UNS removed)
result2 <- readRDS("result2_updated_subype_val.rds")

# Load Absolute promoter activity
abs_pa <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)  
colnames(abs_pa) <- gsub("_SJ.out","",colnames(abs_pa))

# Fetch LSP1 promoters:
keep <- c('10530', '10532','10533')
abs_pa_sel <- abs_pa[rownames(abs_pa) %in% keep,]
abs_pa_sel_t <- t(abs_pa_sel)
colnames(abs_pa_sel_t)<-paste0("pr",colnames(abs_pa_sel_t))

# Merge with metadata & make boxplots:---------------------------------->>
finaldf <- merge(metadata,abs_pa_sel_t,by = 0)
finaldf$bl1_vs_rest <- factor(finaldf$bl1_vs_rest, levels = c("BL1","REST"))

# Reshape data here:
finaldf <- finaldf %>%
  pivot_longer(cols = c(pr10530,pr10532,pr10533), 
               names_to = "promoters", 
               values_to = "expression")
# Define colors
my_colors <- c("BL1"= "darkorange1","REST"="orchid3")
my_colors2 <- c("BL1"= "#FFF2E6","REST"="#F4E1FF")

# Boxplot
g <- ggboxplot(finaldf, x = "promoters", y = "expression", 
               color = "bl1_vs_rest",
               fill = "bl1_vs_rest", size = 1.5,) +
  labs(y = "Absolute promoter activity") + 
  scale_color_manual(values = my_colors) + 
  scale_fill_manual(values = my_colors2) +
  theme(
    text = element_text(size = 25),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, face = "bold"),
    axis.text.x = element_text(size = 24, colour = "black"),
    axis.text.y = element_text(size = 24, colour = "black"),
    legend.text = element_text(size = 25),
    legend.title = element_blank()
  )

g<-g + stat_compare_means(aes(group = bl1_vs_rest), 
                          size=10,label = "p.signif",
                          label.y = 7.3,
                          label.x.npc = "middle") + scale_y_continuous(limits = c(0,8),expand = c(0,0))
g
ggsave("./val/lsp1_pr2.png", g, width = 6.2, height = 5.5, dpi = 1000)

####----

# For LSP1 gene now
# LSP1: ENSG00000130592.17

g <- grep("ENSG00000130592.17",rownames(vMat))
exprMat_sel <- as.data.frame(vMat[g,])
mm <- merge(metadata, exprMat_sel, by=0)
colnames(mm)[6]<- 'Condition'
colnames(mm)[7]<- 'Expression'

mm$Condition <- factor(mm$Condition , levels = c("BL1","REST"))
my_colors <- c("BL1"= "darkorange1","REST"="orchid3")
my_colors2 <- c("BL1"= "#FFF2E6","REST"="#F4E1FF")

g <- ggboxplot(mm, x = "Condition", y = "Expression", 
               color = "Condition", fill = "Condition",size = 1.5) +
  labs(y = "Normalized gene expression") + 
  scale_color_manual(values = my_colors) + 
  scale_fill_manual(values = my_colors2) +
  theme(
    text = element_text(size = 25),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(size = 24, colour = "black"),
    axis.text.y = element_text(size = 24, colour = "black"),
    legend.text = element_text(size = 25),
    legend.title = element_blank()
  )

g<-g + stat_compare_means(aes(group = Condition), 
                          size=12,label = "p.signif",
                          label.y = 7.6,
                          label.x.npc = "middle") + stat_n_text(size = 8,
                                                                color = "black",
                                                                fontface = "bold") +
  scale_y_continuous(limits = c(2.2, 8),expand = c(0,0))

g
ggsave("./val/lsp1_gene2.png", g, width = 5.4, height = 5.6, dpi = 1000)

#############################################################################
#############################################################################
