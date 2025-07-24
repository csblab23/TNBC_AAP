############## This script includes all the codes used for FIGURE PANEL 3 ################

# Load required libraries
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(readxl)
library(qtl2)
library(data.table)
library(proActiv)
library(EnvStats)
library(edgeR)



## Code For Figure 3a:---------------------->>
# Heatmap of the validated Basal-specific Active Alternative Promoters (AAPs)
# Load FUSCC metadata
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
metadata <- metadata[metadata$Tumor_Normal1=="Tumor",]  # just take tumor samples
subset_metadata <- metadata[,c('Run', 'batch','Tumor_Normal1','Intrinsic_Subtype')]   # just keep the required columns

# Read the result2 file obtained by running proactiv on FUSCC intrinsic subtype:
result2 <- readRDS("result2_intrinsic_subtype.rds")

# Use Absolute promoter activity:-------------------------------->>
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

# Keep just the 6 validated Basal-specific AAPs:
keep <- c("2055","4814","7033","10005","12842","30507")

abs_pa_filtered <- abs_pa[rownames(abs_pa) %in% keep,]
abs_pa_filtered <- abs_pa_filtered[keep,]

# Adding gene names corresponding to the selected promoters:
# pr_genes <- c('pr2055_MARK2','pr4814_PIEZO1',"pr7033_SLC29A1","pr10005_AKAP9","pr12842_SEC31A","pr30507_PHACTR4")
genes_df <- as.data.frame(read_csv('selected_6aaps_and_genes.csv'))
genes_df<-genes_df[,-c(1,3)]

abs_pa_filtered <- tibble::rownames_to_column(abs_pa_filtered,"promoterId")
abs_pa_filtered2 <- merge(genes_df, abs_pa_filtered, by="promoterId")
abs_pa_filtered2 <- abs_pa_filtered2[,-1]
abs_pa_filtered2 <- tibble::column_to_rownames(abs_pa_filtered2,"genes")

#####------

count_norm <- as.data.frame(t(abs_pa_filtered2))
count_norm <- tibble::rownames_to_column(count_norm,'Run')
subset_metadata$Intrinsic_Subtype <- gsub("Other","Non_Basal",subset_metadata$Intrinsic_Subtype)
all_info <- merge(subset_metadata,count_norm,by="Run")
rownames(all_info) <- all_info$Run
meta_req <- all_info[,1:4]
exp_only <- all_info[,5:dim(all_info)[2]]

exp_only <- exp_only[meta_req$Run,]
Samples=meta_req$Run
Sample_type=meta_req$Intrinsic_Subtype

# Create a data frame for column annotation
ann_df <- data.frame(Samples,Sample_type)

# Order the dataframe
ann_df_ordered <- ann_df[order(ann_df$Sample_type), ]
rownames(ann_df_ordered) <- NULL
ann_df_ordered <- tibble::column_to_rownames(ann_df_ordered,'Samples')

ann_colors <- list( "Sample_type" = c("Basal" = "#4ecb8d", "Non_Basal" = "honeydew4"))
col_fun = colorRamp2(c(-3,-1.5,-0.5,0, 0.5, 1.5,3), c("darkgoldenrod3","gold1","khaki","white" , "thistle2","hotpink1","darkmagenta"))

exp_only <- as.matrix(exp_only)
scaled <- scale(exp_only)
scaled <- t(scaled)
scaled <- scaled[,rownames(ann_df_ordered)]
scaled <- na.omit(scaled)

# Combine the heatmap and the annotation
h1 <- Heatmap(scaled,name="Scaled\nExpression",column_names_gp = grid::gpar(fontsize = 5),
              row_names_gp = grid::gpar(fontsize = 15), height = unit(12, "cm"),  width = unit(10,"cm"),
              col = col_fun,top_annotation = HeatmapAnnotation(df = ann_df_ordered ,col = ann_colors,which = 'col',
                                                               annotation_name_gp = grid::gpar(fontsize= 15,fontface = "bold")),
              show_row_names = T,show_column_names = F, cluster_rows = TRUE, cluster_columns = FALSE,
              heatmap_legend_param = list(legend_gp = gpar(fontsize = 16)))

save <- draw(h1,
             heatmap_legend_side = 'right',
             annotation_legend_side = 'right',
             row_sub_title_side = 'left',
             column_title_gp= gpar(fontsize = 15, fontface = 'bold'),
             column_title="Basal specific AAPs"
)
save
png("Heatmap_figure3A.png",width = 1500,height = 1000,res = 150)
save
dev.off()

############################################################################
############################################################################

## Code for Figure 3b,c,f,g Box plots:---->>

# For FUSCC:
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

# Take tumor samples only:
metadata <- metadata[metadata$Tumor_Normal1=='Tumor',]
subset_metadata <- metadata[,c('Run','batch','Intrinsic_Subtype')] # just keep the required columns

# Read the feature-counts matrix 
count <- as.data.frame(read_csv("Featurecounts_combined_fuscc.csv"))
count = count[,-1] # remove extra s.no. column
count <- tibble::column_to_rownames(count, var = "Geneid")

# Keep just same samples in count matrix as in meta:
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

# Read the result2 file:
result2 <- readRDS("result2_intrinsic_subtype.rds")  # (for boxplot)

# Box plots:--------->>
## Using the Absolute promoter activity by fetching from result2 file
abs_pa <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)  
colnames(abs_pa) <- gsub("_SJ.out","",colnames(abs_pa))
abs_pa <- abs_pa[,colnames(abs_pa)%in%subset_metadata$Run]   # keep same samples as in the metadata
abs_pa <-abs_pa[,subset_metadata$Run]

# Define a named list of genes and their promoter IDs
promoters_list <- list(
  AKAP9 = c("10004", "10005"),
  SEC31A = c("12841", "12842")
)

# Loop over each gene and its promoters
for (gene in names(promoters_list)) {
  
  keep <- promoters_list[[gene]]  # get the promoter IDs for this gene
  
  # Filter expression matrix by promoter IDs
  abs_pa_sel <- abs_pa[rownames(abs_pa) %in% keep, ]
  abs_pa_sel_t <- t(abs_pa_sel)
  colnames(abs_pa_sel_t) <- paste0("pr", colnames(abs_pa_sel_t))
  
  # Merge with metadata
  finaldf <- merge(subset_metadata, abs_pa_sel_t, by = 0)
  finaldf$Intrinsic_Subtype <- gsub("Other", "Non_Basal", finaldf$Intrinsic_Subtype)
  finaldf$Intrinsic_Subtype <- factor(finaldf$Intrinsic_Subtype, levels = c("Basal", "Non_Basal"))
  
  # Reshape to long format
  finaldf <- finaldf %>%
    pivot_longer(cols = starts_with("pr"), 
                 names_to = "promoters", 
                 values_to = "expression")
  
  # Define colors
  my_colors <- c("Basal" = "#4ecb8d", "Non_Basal" = "honeydew4")
  
  # Boxplot
  g <- ggboxplot(finaldf, x = "promoters", y = "expression", 
                 color = "Intrinsic_Subtype", size = 1.5) +
    labs(y = "Absolute promoter activity", title = gene) +
    scale_color_manual(values = my_colors) + 
    theme(
      text = element_text(size = 25),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 24, face = "bold"),
      axis.text.x = element_text(size = 24, colour = "black"),
      axis.text.y = element_text(size = 24, colour = "black"),
      legend.text = element_text(size = 25),
      legend.title = element_blank()
    ) +
    stat_compare_means(aes(group = Intrinsic_Subtype), 
                       size = 12, label = "p.signif",
                       label.y = 8.6, 
                       label.x.npc = "middle") +
    scale_y_continuous(limits = c(4, 9), expand = c(0, 0))
  
  # Print and save
  print(g)
  ggsave(paste0(gene, "_plot.png"), g, width = 5.5, height = 5.5, dpi = 1000)
}

##########--------------

## Box plots similarly for genes now:
# Define a named list of genes with Ensembl IDs
gene_list <- list(
  AKAP9 = "ENSG00000127914.19",
  SEC31A = "ENSG00000138674.17"
)

# Loop through each gene
for (gene_name in names(gene_list)) {
  
  ensg_id <- gene_list[[gene_name]]
  
  # Extract expression for that gene
  g_idx <- grep(ensg_id, rownames(exprMat))
  if (length(g_idx) == 0) {
    warning(paste("Gene", gene_name, "not found in exprMat. Skipping."))
    next
  }
  
  exprMat_sel <- as.data.frame(exprMat[g_idx, ])
  
  # Merge with metadata
  mm <- merge(subset_metadata, exprMat_sel, by = 0)
  colnames(mm)[4] <- "Condition"  
  colnames(mm)[5] <- "Expression"
  
  # Clean up labels
  mm$Condition <- gsub("Other", "Non_Basal", mm$Condition)
  mm$Condition <- factor(mm$Condition, levels = c("Basal", "Non_Basal"))
  
  # Define colors
  my_colors <- c("Basal" = "#4ecb8d", "Non_Basal" = "honeydew4")
  
  # Boxplot
  g <- ggboxplot(mm, x = "Condition", y = "Expression", 
                 color = "Condition", size = 1.5) +
    labs(y = "Normalized gene expression", title = gene_name) + 
    scale_color_manual(values = my_colors) +
    theme(
      text = element_text(size = 25),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 23, face = "bold"),
      axis.text.x = element_text(size = 24, colour = "black"),
      axis.text.y = element_text(size = 24, colour = "black"),
      legend.text = element_text(size = 25),
      legend.title = element_blank()
    ) +
    stat_compare_means(aes(group = Condition),
                       size = 12, label = "p.signif",
                       label.y = 9.5,
                       label.x.npc = "middle") +
    stat_n_text(size = 8, color = "black", fontface = "bold") +
    scale_y_continuous(limits = c(5.4, 10), expand = c(0, 0))
  
  # Print and save plot
  print(g)
  ggsave(paste0(gene_name, "_gene_plot.png"), g, width = 5.4, height = 5.6, dpi = 1000)
}

#############################################################################
#############################################################################

## Code for Figure 3d,e,h,i Box plots:---->>

# These are the box plots For Validation dataset/external dataset 1
# Load the metafile  
metadata <- as.data.frame(read_csv("TNBC_mergedmeta.csv"))
metadata=tibble::column_to_rownames(metadata,var = '...1')

# Read the feature-counts matrix 
count <- as.data.frame(read_csv("Raw_counts_61samples.csv"))
count=count[,-1] # remove extra s.no. column
count <- tibble::column_to_rownames(count, var = "Geneid")

# Keep same samples in count matrix as in metadata:
count <- count[,colnames(count) %in% metadata$Run]
count <- count[,metadata$Run]  # match order of samples 

##--

# Now calculate cpm  
cpm <- cpm(count)
is.exprs <- rowSums(cpm>1) >= 5 # Apply filtering criteria
counts2 <- count[is.exprs, ]

# Prepare DGEList and voom normalize the matrix
x <- DGEList(counts=counts2)
x <- calcNormFactors(x, method='TMM')
v <- voom(x, plot=T)
vMat <- v$E

# Read the result2 file obtained by running proactiv 
result2 <- readRDS("result2_35basal_vs_5others.rds")

# Use Absolute promoter activity:-------------------------------->>
abs_pa <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)  
colnames(abs_pa) <- gsub("_SJ.out","",colnames(abs_pa))
abs_pa <- abs_pa[,colnames(abs_pa)%in%metadata$Run]
abs_pa <-abs_pa[,metadata$Run]

# Define list of promoter sets (per gene)
promoters_list <- list(
  AKAP9 = c("10004", "10005"),
  SEC31A = c("12841", "12842")
)

# Define colors
my_colors  <- c("Basal" = "#4ecb8d", "Non_Basal" = "honeydew4")
my_colors2 <- c("Basal" = "#c8f4d9", "Non_Basal" = "snow2")

# Loop through each gene and promoter set
for (gene_name in names(promoters_list)) {
  
  keep <- promoters_list[[gene_name]]
  
  # Subset and transpose promoter activity matrix
  abs_pa_sel <- abs_pa[rownames(abs_pa) %in% keep, ]
  abs_pa_sel_t <- t(abs_pa_sel)
  colnames(abs_pa_sel_t) <- paste0("pr", colnames(abs_pa_sel_t))
  
  # Merge with metadata
  finaldf <- merge(metadata, abs_pa_sel_t, by = 0)
  
  # Factor subtype
  finaldf$intrinsic_subtype <- gsub("Other", "Non_Basal", finaldf$intrinsic_subtype)
  finaldf$intrinsic_subtype <- factor(finaldf$intrinsic_subtype, levels = c("Basal", "Non_Basal"))
  
  # Reshape to long format for ggplot
  finaldf <- finaldf %>%
    pivot_longer(cols = starts_with("pr"),
                 names_to = "promoters", 
                 values_to = "expression")
  
  # Boxplot
  g <- ggboxplot(finaldf, x = "promoters", y = "expression", 
                 color = "intrinsic_subtype",
                 fill = "intrinsic_subtype",
                 size = 1.5) +
    labs(title = gene_name, y = "Absolute promoter activity") + 
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
    ) +
    stat_compare_means(aes(group = intrinsic_subtype),
                       size = 10,
                       label = "p.signif",
                       label.y = 9.3,
                       label.x.npc = "middle") +
    scale_y_continuous(limits = c(0, 10), expand = c(0, 0))
  
  # Display and save
  print(g)
  ggsave(paste0(gene_name, "_plot.png"), g, width = 5.5, height = 5.5, dpi = 1000)
}

#####----

## Box plots similarly for gene now:
# Define list of genes: SYMBOL = ENSEMBL ID
gene_list <- list(
  AKAP9 = "ENSG00000127914.19",
  SEC31A = "ENSG00000138674.17"
)

# Color definitions
my_colors  <- c("Basal" = "#4ecb8d", "Non_Basal" = "honeydew4")
my_colors2 <- c("Basal" = "#c8f4d9", "Non_Basal" = "snow2")

# Loop through genes
for (gene_name in names(gene_list)) {
  
  ensg_id <- gene_list[[gene_name]]
  
  # Get expression row from vMat
  g_idx <- grep(ensg_id, rownames(vMat))
  if (length(g_idx) == 0) {
    warning(paste("Gene", gene_name, "not found in vMat. Skipping."))
    next
  }
  
  exprMat_sel <- as.data.frame(vMat[g_idx, ])
  
  # Merge with metadata
  mm <- merge(metadata, exprMat_sel, by = 0)
  
  # Adjust column names
  colnames(mm)[6]<- 'Condition'
  colnames(mm)[7]<- 'Expression'
  
  # Clean up condition labels
  mm$Condition <- gsub("Other", "Non_Basal", mm$Condition)
  mm$Condition <- factor(mm$Condition, levels = c("Basal", "Non_Basal"))
  
  # Boxplot
  g <- ggboxplot(mm, x = "Condition", y = "Expression", 
                 color = "Condition", fill = "Condition", size = 1.5) +
    labs(title = gene_name, y = "Normalized gene expression") + 
    scale_color_manual(values = my_colors) + 
    scale_fill_manual(values = my_colors2) +
    theme(
      text = element_text(size = 25),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 23, face = "bold"),
      axis.text.x = element_text(size = 24, colour = "black"),
      axis.text.y = element_text(size = 24, colour = "black"),
      legend.text = element_text(size = 25),
      legend.title = element_blank()
    ) +
    stat_compare_means(aes(group = Condition), 
                       size = 12, label = "p.signif",
                       label.y = 12.5,
                       label.x.npc = "middle") +
    stat_n_text(size = 8, color = "black", fontface = "bold") +
    scale_y_continuous(limits = c(7, 13), expand = c(0, 0))
  
  # Print and save plot
  print(g)
  ggsave(paste0(gene_name, "_gene_plot.png"), g, width = 5.4, height = 5.6, dpi = 1000)
}

#########################################################################
#########################################################################

# Code for Figure 3j & k:----------->>

# Load the FUSCC metafile with subtype information  
metadata <- as.data.frame(read_csv("fuscc_merged_meta.csv"))
metadata=tibble::column_to_rownames(metadata,var = '...1')
subset_metadata <- metadata

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

# Read the result2 file:
result2 <- readRDS("result2_updated_subype_fuscc.rds")

# Box plots:--------->>
# Load Absolute promoter activity:--------------------->>
abs_pa <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)  
colnames(abs_pa) <- gsub("_SJ.out","",colnames(abs_pa))

# Match the number and order of samples:
abs_pa <- abs_pa[,colnames(abs_pa)%in%subset_metadata$Run]
abs_pa <-abs_pa[,subset_metadata$Run]

# Now fetch required promoters:
keep <- c('10530', '10532','10533')  # for LSP1
abs_pa_sel <- abs_pa[rownames(abs_pa) %in% keep,]
abs_pa_sel_t <- t(abs_pa_sel)
colnames(abs_pa_sel_t)<-paste0("pr",colnames(abs_pa_sel_t))

# Merge with metadata & make boxplots:---->>
finaldf <- merge(subset_metadata,abs_pa_sel_t,by = 0)
finaldf$bl1_vs_rest <- factor(finaldf$bl1_vs_rest, levels = c("BL1","REST"))

# Reshape data here:
finaldf <- finaldf %>%
  pivot_longer(cols = c(pr10530,pr10532,pr10533), 
               names_to = "promoters", 
               values_to = "expression")

# Assign colours
my_colors <- c("BL1"= "darkorange1","REST"="orchid3")

# Plot:
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
ggsave("lsp1_pr_plot.png", g, width = 6, height = 5.5, dpi = 1000)

####----

## For LSP1 gene now
# LSP1: ENSG00000130592.17

g <- grep("ENSG00000130592.17",rownames(exprMat))
exprMat_sel <- as.data.frame(exprMat[g,])
mm <- merge(subset_metadata, exprMat_sel, by=0)
colnames(mm)[5]<- 'Condition'
colnames(mm)[7]<- 'Expression'

# Define levels and assign colors
mm$Condition <- factor(mm$Condition , levels = c("BL1","REST"))
my_colors <- c("BL1"= "darkorange1","REST"="orchid3")

# Boxplot
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
ggsave("lsp1_gene_plot.png", g, width = 5.4, height = 5.6, dpi = 1000)

##############################################################################
##############################################################################

# Code for Figure 3l & m:-------------------->>

# Load the validation dataset metafile with subtype information
metadata <- as.data.frame(read_csv("val_merged_meta.csv"))
metadata <- tibble::column_to_rownames(metadata,var = "...1")

# Read the feature-counts matrix 
count <- as.data.frame(read_csv("Raw_counts_61samples.csv"))
count=count[,-1] # remove extra s.no. column
count <- tibble::column_to_rownames(count, var = "Geneid")

# Keep just same samples in count matrix as in meta:
count <- count[,colnames(count) %in% metadata$Run]
count <- count[,metadata$Run]  # make order of samples same 

##--

# Now calculate cpm  
cpm <- cpm(count)
is.exprs <- rowSums(cpm>1) >= 20 # Apply filtering criteria
counts2 <- count[is.exprs, ]

# Prepare DGEList and voom normalize the matrix
x <- DGEList(counts=counts2)
x <- calcNormFactors(x, method='TMM')
v <- voom(x, plot=T)
vMat <- v$E

##--

# Read the result2 file obtained by running proactiv on tnbc type4 (UNS removed)
result2 <- readRDS("result2_updated_subype_val.rds")

# Load Absolute promoter activity:-------------------------------->>
abs_pa <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)  
colnames(abs_pa) <- gsub("_SJ.out","",colnames(abs_pa))

# Fetch LSP1 promoters:
keep <- c('10530', '10532','10533')  
abs_pa_sel <- abs_pa[rownames(abs_pa) %in% keep,]
abs_pa_sel_t <- t(abs_pa_sel)
colnames(abs_pa_sel_t)<-paste0("pr",colnames(abs_pa_sel_t))

# Merge with metadata & make boxplots:----->>
finaldf <- merge(metadata,abs_pa_sel_t,by = 0)
finaldf$bl1_vs_rest <- factor(finaldf$bl1_vs_rest, levels = c("BL1","REST"))

# Reshape data here:
finaldf <- finaldf %>%
  pivot_longer(cols = c(pr10530,pr10532,pr10533), 
               names_to = "promoters", 
               values_to = "expression")

# Define levels and assign colors:
finaldf$bl1_vs_rest <- factor(finaldf$bl1_vs_rest, levels = c("BL1","REST"))
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
g <- g + stat_compare_means(aes(group = bl1_vs_rest), 
                          size=10,label = "p.signif",
                          label.y = 7.3,
                          label.x.npc = "middle") + scale_y_continuous(limits = c(0,8),expand = c(0,0))
g
ggsave("lsp1_pr2.png", g, width = 6.2, height = 5.5, dpi = 1000)

####----

## For LSP1 gene
# LSP1: ENSG00000130592.17

g <- grep("ENSG00000130592.17",rownames(vMat))
exprMat_sel <- as.data.frame(vMat[g,])
mm <- merge(metadata, exprMat_sel, by=0)
colnames(mm)[6]<- 'Condition'
colnames(mm)[7]<- 'Expression'

# Define levels and assign colors
mm$Condition <- factor(mm$Condition , levels = c("BL1","REST"))
my_colors <- c("BL1"= "darkorange1","REST"="orchid3")
my_colors2 <- c("BL1"= "#FFF2E6","REST"="#F4E1FF")

# Boxplot:
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
g <- g + stat_compare_means(aes(group = Condition), 
                          size=12,label = "p.signif",
                          label.y = 7.6,
                          label.x.npc = "middle") + stat_n_text(size = 8,
                                                                color = "black",
                                                                fontface = "bold") +
  scale_y_continuous(limits = c(2.2, 8),expand = c(0,0))
g
ggsave("lsp1_gene.png", g, width = 5.4, height = 5.6, dpi = 1000)

#############################################################################
#############################################################################