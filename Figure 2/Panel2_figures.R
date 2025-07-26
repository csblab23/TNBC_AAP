############## This script includes all the codes used for FIGURE PANEL 2 ################

# Load required libraries
library(dplyr)
library(tidyverse)
library(data.table)
library(ggpubr)
library(proActiv)
library(EnvStats)
library(ggvenn)
library(eulerr)
library(edgeR)

## For Figure 2a:------------>>

# Venn diagram for genes having both upregulated and downregulated DRPs:
# Load the DEG wilcoxon result and keep just the nonsignificant genes:
degs <- as.data.frame(read_csv('DEG_wilcoxon_tum_vs_nor.csv'))
colnames(degs)[1] <- 'gene_id'
degs_nonsig <- degs[degs$pValues>0.05,]

# Fetch significant promoters for them which are both upregulated and downregulated
# NOTE: considering only absolute promoter activity here
# Load the upregulated DRP dataframe (pvalue<0.05)
up_drp <- as.data.frame(read_csv("upreg_drp_absPA_bypvalue.csv"))  
colnames(up_drp)[1] <- 'promoterId'  # renaming the column 1

# Load the downregulated DRP dataframe (pvalue<0.05)
down_drp <- as.data.frame(read_csv("down_drp_absPA_bypvalue.csv"))
colnames(down_drp)[1] <- 'promoterId'

## Add the gene ids in the above dataframes:
# Read the result2 file:
result2 <- readRDS("result2_tum_vs_adjnormal.rds")

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

# First keep only the nonsignificant genes in this:
filtered_df_nonsig <- filtered_df[filtered_df$geneId%in% degs_nonsig$gene_id,]

## Preparing separate dataframes for up drp and down drp:
# For upreg drp:
merged_up <- merge(up_drp,filtered_df_nonsig,by='promoterId')
# For downreg drp:
merged_down <- merge(down_drp,filtered_df_nonsig,by='promoterId')
# Now find common genes:
common_genes <- intersect(merged_up$geneId, merged_down$geneId)

# Prepare proportionate venn diagram
up_drps <- unique(merged_up$geneId) 
down_drps <- unique(merged_down$geneId)
venn <- list(UP_DRPs = up_drps, DOWN_DRPs = down_drps)
euler_data <- euler(venn)
p1 <- plot(euler_data,fills = list(fill = c("salmon3", "bisque2"),alpha = 0.7),
           labels = list(col = "black",cex=2.1,hjust=0.5),quantities = list(cex = 2.2,hjust=0.5),cex = 2.5)
p1
ggsave("Fig_both_up_down_drps_venn.png", p1, width = 6, height = 4, dpi = 1000)

#######################################################################################
#######################################################################################

## For Figure 2b & c:------------>>

## Make box plots for the HDAC9 gene and promoter:
# Load the FUSCC metafile
metadata <- as.data.frame(read_csv("meta_with_batchinfo.csv"))
metadata <- metadata[,-1]
metadata$Tumor_Normal1=metadata$Tumor_Normal
# Cleaning the tissue type column
metadata =metadata %>% 
  mutate(Tumor_Normal1 = case_when(Tumor_Normal=="TumorTissue"  ~ "Tumor",
                                   Tumor_Normal=="TumorTissue.2018"  ~ "Tumor",
                                   Tumor_Normal=="TumorTissue.rep"  ~ "Tumor",
                                   Tumor_Normal=="PairedNormalTissue"  ~ "AdjNormal",
                                   Tumor_Normal=="PairedNormalTissue.rep"  ~ "AdjNormal"))
rownames(metadata) <- metadata$Run  # making the sample names as rownames

# Subset to just keep the required columns
subset_meta <- metadata[,c('Run', 'batch','Tumor_Normal1')]   
merged_meta <- subset_meta

# Read the feature-counts matrix 
count <- as.data.frame(read_csv("Featurecounts_combined_fuscc.csv"))
count=count[,-1] # remove extra s.no. column
count <- tibble::column_to_rownames(count, var = "Geneid")

# Now calculate cpm  
cpm <- cpm(count)
is.exprs <- rowSums(cpm>1) >= 20 # Apply filtering criteria
counts2 <- count[is.exprs, ]

# Make samples order same in both the files 
com_exp_mat=counts2[,merged_meta$Run]

# Prepare DGEList and voom normalize the matrix
x <- DGEList(counts=com_exp_mat)
x <- calcNormFactors(x, method='TMM')
v <- voom(x, plot=T)
vMat <- v$E

# Adjust batch by sva combat here
merged_meta$batch=as.factor(merged_meta$batch)
exprMat <- sva::ComBat(vMat, batch=merged_meta$batch)

## Box plot for HDAC9 gene expression ---->>
# HDAC9: ENSG00000048052.25
# Retain the gene of interest in the exprMat:
g <- grep("ENSG00000048052.25",rownames(exprMat))
exprMat_sel <- as.data.frame(exprMat[g,])
mm <- merge(merged_meta, exprMat_sel, by=0)
colnames(mm)[4]<- 'Condition'
colnames(mm)[5]<- 'Expression'
mm$Condition <- gsub("Tumor","TNBC",mm$Condition)
my_colors <- c("AdjNormal"= "steelblue1","TNBC"= "tomato2")
mm$Condition <- factor(mm$Condition, levels = c("TNBC", "AdjNormal"))

# Box plot
g <- ggboxplot(mm, x = "Condition", y = "Expression", 
               color = "Condition", 
               fill = NA, size = 1.6) +
  labs(y = "Normalized gene expression") + 
  scale_color_manual(values = my_colors) + 
  theme(
    text = element_text(size = 25),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text.x = element_text(size = 24, colour = "black"),
    axis.text.y = element_text(size = 24, colour = "black"),
    legend.text = element_text(size = 25),
    legend.title = element_blank()
  )

g<-g + stat_compare_means(aes(group = Condition), 
                          size=12,label = "p.signif",
                          label.y = 7.5,
                          label.x.npc = "middle")+ stat_n_text(size = 7,
                                                               color = "black",
                                                               fontface = "bold")+ scale_y_continuous(limits = c(2,8),expand = c(0,0))

g
ggsave("HDAC9_gene.png", g, width = 7, height = 7, dpi = 1000)

# Box plot for HDAC9 promoters:----->>
# Read the result2 file:
result2 <- readRDS("result2_tum_vs_adjnormal.rds")

# Fetch the absolute promoter activity 
abs_pa <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)  
colnames(abs_pa) <- gsub("_SJ.out","",colnames(abs_pa))

# Match the order of samples
abs_pa <- abs_pa[,colnames(abs_pa)%in%merged_meta$Run]
abs_pa <-abs_pa[,merged_meta$Run]

# Fetch HDAC9 promoters
keep <- c('1077','1079')

# Subset in promoter activity dataframe
abs_pa_sel <- abs_pa[rownames(abs_pa) %in% keep,]
abs_pa_sel_t <- t(abs_pa_sel)
colnames(abs_pa_sel_t)<-paste0("pr",colnames(abs_pa_sel_t))

# Merge with metadata & make boxplot:----->>
finaldf <- merge(merged_meta,abs_pa_sel_t,by = 0)
finaldf$Tumor_Normal1 <- gsub("Tumor","TNBC",finaldf$Tumor_Normal1)
finaldf$Tumor_Normal1 <- factor(finaldf$Tumor_Normal1,levels = c("TNBC","AdjNormal"))

# Reshape data here:
finaldf <- finaldf %>%
  pivot_longer(cols = c(pr22914,pr22915), 
               names_to = "promoters", 
               values_to = "expression")
my_colors <- c("TNBC"= "tomato2","AdjNormal"= "steelblue1")

g <- ggboxplot(finaldf, x = "promoters", y = "expression", 
               color = "Tumor_Normal1", 
               fill = NA,size = 1.6) +
  labs(y = "Absolute promoter activity") + 
  scale_color_manual(values = my_colors) + 
  theme(
    text = element_text(size = 25),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text.x = element_text(size = 24, colour = "black"),
    axis.text.y = element_text(size = 24, colour = "black"),
    legend.text = element_text(size = 25),
    legend.title = element_blank()
  )

g<-g + stat_compare_means(aes(group = Tumor_Normal1), 
                          size=12,
                          label = "p.signif",
                          label.y = 7,
                          label.x.npc = "middle") + scale_y_continuous(limits = c(-0.2,8),expand = c(0,0))
g
ggsave("HDAC9_pr_fuscc.png", g, width = 6, height = 6, dpi = 1000)

##########################################################################
##########################################################################

## For Figure 2d & e:------------>>

# Similarly making HDAC9 box plots for the external dataset 1
# Load the metafile
metadata <- as.data.frame(read_csv("meta_sel_61samples.csv"))
metadata <- tibble::column_to_rownames(metadata,var = "...1")

# Read the feature counts matrix 
count <- as.data.frame(read_csv("Raw_counts_61samples.csv"))
count <- count[,-1]
count <- tibble::column_to_rownames(count,'Geneid')

# Now calculate cpm  
cpm <- cpm(count)
is.exprs <- rowSums(cpm>1) >= 5 # Apply filtering criteria
counts <- count[is.exprs, ]

# Make sample order same in both the files 
com_exp_mat=counts[,metadata$Run]

# Prepare DGEList and voom normalize the matrix
x <- DGEList(counts=com_exp_mat)
x <- calcNormFactors(x, method='TMM')
v <- voom(x, plot=T)
vMat <- v$E

## Box plot for HDAC9 gene expression
# HDAC9: ENSG00000048052.25
# Keep gene of interest in vMat
g <- grep("ENSG00000048052.25",rownames(vMat))
exprMat_sel <- as.data.frame(vMat[g,])
rownames(metadata) <- metadata$Run
mm <- merge(metadata, exprMat_sel, by=0)
colnames(mm)[4]<- 'Condition'
colnames(mm)[5]<- 'Expression'
mm$Condition <- gsub("Normal","AdjNormal", mm$Condition)

# Specify colors
my_colors1 <- c("TNBC"= "tomato2","AdjNormal"= "steelblue1")
my_colors2 <- c("TNBC"= "#FFE6E5","AdjNormal"= "#EAFAFF")
mm$Condition <- factor(mm$Condition, levels = c("TNBC", "AdjNormal"))
g <- ggboxplot(mm, x = "Condition", y = "Expression", 
               color = "Condition", 
               fill = "Condition", size = 1.6) +
  labs(y = "Normalized gene expression") + 
  scale_color_manual(values = my_colors1) +
  scale_fill_manual(values = my_colors2) +
  theme(
    text = element_text(size = 25),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text.x = element_text(size = 24, colour = "black"),
    axis.text.y = element_text(size = 24, colour = "black"),
    legend.text = element_text(size = 25),
    legend.title = element_blank()
  )

g<-g + stat_compare_means(aes(group = Condition), 
                          size=14,label = "p.signif",
                          label.y = 8.6,
                          label.x.npc = "middle")+ stat_n_text(size = 7,
                                                               color = "black",
                                                               fontface = "bold") 
g
ggsave("HDAC9_gene_new.png", g, width = 5.5, height = 6, dpi = 1000)

##--

## Box plot for HDAC9 promoters
# Read the result2 file:
result2 <- readRDS("result2_40tnbc_vs_21normal.rds")

# Fetch the Absolute promoter activity:---->>
abs_pa <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)  
colnames(abs_pa) <- gsub("_SJ.out","",colnames(abs_pa))
abs_pa <- abs_pa_filtered  

# Now keep HDAC9 promoters
keep <- c('1077','1079')

# Retain the same promoters in absolute promoter activity dataframe
abs_pa_sel <- abs_pa[rownames(abs_pa) %in% keep,]
abs_pa_sel_t <- t(abs_pa_sel)
colnames(abs_pa_sel_t)<-paste0("pr",colnames(abs_pa_sel_t))

# Merge with metadata & make boxplots:----->>
finaldf <- merge(metadata,abs_pa_sel_t,by = 0)
finaldf$tumor_normal <- gsub("Normal","AdjNormal",finaldf$tumor_normal)
finaldf$tumor_normal <- factor(finaldf$tumor_normal,levels = c("TNBC","AdjNormal"))

# Reshape data here:
finaldf <- finaldf %>%
  pivot_longer(cols = c(pr68402,pr68403), 
               names_to = "promoters", 
               values_to = "expression")

my_colors1 <- c("TNBC"= "tomato2","AdjNormal"= "steelblue1")
my_colors2 <- c("TNBC"= "#FEE1E0","AdjNormal"= "#E5F8FF")
g <- ggboxplot(finaldf, x = "promoters", y = "expression", 
               color = "tumor_normal", 
               fill = "tumor_normal", size = 1.6) +
  labs(y = "Absolute promoter activity") + 
  scale_color_manual(values = my_colors1) +
  scale_fill_manual(values = my_colors2) +
  theme(
    text = element_text(size = 25),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text.x = element_text(size = 24, colour = "black"),
    axis.text.y = element_text(size = 24, colour = "black"),
    legend.text = element_text(size = 25),
    legend.title = element_blank()
  )

g<-g + stat_compare_means(aes(group = tumor_normal), 
                          size=14,label = "p.signif",
                          label.y = 6,
                          label.x.npc = "middle") + scale_y_continuous(limits = c(-0.6, 7),expand = c(0,0))
g
ggsave("HDAC9_pr_val.png", g, width = 6, height = 6, dpi = 1000)

#############################################################################
#############################################################################

## For Figure 2g & h:------------>>

# Box plots for the HDAC9 pr1077 & pr1079 transcript expression in FUSCC
# Load the metafile  
metadata <- as.data.frame(read_csv("meta_with_batchinfo.csv"))
metadata <- metadata[,-1]
metadata$Tumor_Normal1=metadata$Tumor_Normal
# Adjust the tissue type column
metadata =metadata %>% 
  mutate(Tumor_Normal1 = case_when(Tumor_Normal=="TumorTissue"  ~ "Tumor",
                                   Tumor_Normal=="TumorTissue.2018"  ~ "Tumor",
                                   Tumor_Normal=="TumorTissue.rep"  ~ "Tumor",
                                   Tumor_Normal=="PairedNormalTissue"  ~ "AdjNormal",
                                   Tumor_Normal=="PairedNormalTissue.rep"  ~ "AdjNormal"))

rownames(metadata) <- metadata$Run  # making the sample names as rownames
subset_meta <- metadata[,c('Project_ID','Run', 'batch','Tumor_Normal1')]   # just keep the required columns

# Read the result2 file for all 448 samples obtained through proActiv
result2 <- readRDS("./tum_nor_alternate_prom_result/result2_tum_vs_adjnormal.rds")
result2_promoter=rowData(result2)
result2_promoter=as.data.frame(result2_promoter)

# Subset HDAC9
result2_sel <- result2_promoter[result2_promoter$geneId=="ENSG00000048052.25",]
all_hdac9_tr_flat <- unlist(result2_sel$txId, use.names = FALSE)

## Load the transcript FPKM df:
tr_fpkm <- as.data.frame(read_csv('DF_FPKM_transcript.csv'))
tr_fpkm <- tr_fpkm[,-1]
colnames(tr_fpkm)<-gsub('.isoforms.results','',colnames(tr_fpkm))
tr_fpkm<-tibble::column_to_rownames(tr_fpkm,var = 'transcript_id')

# Log2(FPKM+1)
tr_fpkm <- log2(tr_fpkm + 1)

# Keep same samples as in metadata:
tr_fpkm_sel <- tr_fpkm[,colnames(tr_fpkm)%in%subset_meta$Run]
tr_fpkm_sel<-tr_fpkm_sel[,subset_meta$Run]

# Adjust batch by sva combat here----------------
subset_meta$batch=as.factor(subset_meta$batch)
exprMat <- sva::ComBat(tr_fpkm_sel, batch=subset_meta$batch)

# Merge with metadata
merged_all <- merge(exprMat,subset_meta,by=0)
merged_all <- tibble::column_to_rownames(merged_all,"Row.names")

# For pr1079 transcripts:
merged_all_sel <- merged_all[,c(9,10,19:22)]
merged_all_sel <- as.data.frame(merged_all_sel)

## Box plot:
# Pivot to long format
long_df <- merged_all_sel %>%
  pivot_longer(
    cols = starts_with("ENST"),
    names_to = "Transcript_ID",
    values_to = "Expression"
  )
long_df$Tumor_Normal1 <- gsub("Tumor","TNBC",long_df$Tumor_Normal1)
my_colors <- c("TNBC"= "tomato2","AdjNormal"= "steelblue1")
long_df$Tumor_Normal1 <- factor(long_df$Tumor_Normal1,levels = c("TNBC","AdjNormal"))

# Boxplot of expression per transcript, grouped by Tumor/Normal
g1 <- ggplot(long_df, aes(x = Tumor_Normal1, y = Expression, color = Tumor_Normal1)) +
  geom_boxplot(outlier.size = 0.5,size = 1.6) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ Transcript_ID) +
  theme_bw() +
  theme(
    text = element_text(size = 25),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text.x = element_text(size = 24, colour = "black"),
    axis.text.y = element_text(size = 24, colour = "black"),
    strip.text = element_text(size = 22,colour = "black",face ="bold"),
    legend.text = element_blank(),
    legend.title = element_blank()) +
  labs(x = "Sample Type", y = "Transcript Expression")  +
  coord_cartesian(ylim = c(-0.6, 5)) +
  stat_compare_means(aes(group = Tumor_Normal1),
                     size=16,
                     label = "p.signif",
                     label.y = 4.1,label.x.npc = "middle") + stat_n_text(size = 6,color = "black",
                                                                         fontface = "bold")
g1
ggsave("HDAC9_expression_boxplot_pr1079.png", plot = g1, width = 8, height = 6, units = "in",
       bg = "white")

# For pr1077 transcripts:
merged_all_sel <- merged_all[,c(5,15,19:22)]  
merged_all_sel <- as.data.frame(merged_all_sel)

## Box plot:
# Pivot to long format
long_df <- merged_all_sel %>%
  pivot_longer(
    cols = starts_with("ENST"),
    names_to = "Transcript_ID",
    values_to = "Expression"
  )

long_df$Tumor_Normal1 <- gsub("Tumor","TNBC",long_df$Tumor_Normal1)
my_colors <- c("TNBC"= "tomato2","AdjNormal"= "steelblue1")
long_df$Tumor_Normal1 <- factor(long_df$Tumor_Normal1,levels = c("TNBC","AdjNormal"))

# Boxplot of expression per transcript, grouped by Tumor/Normal
g1 <- ggplot(long_df, aes(x = Tumor_Normal1, y = Expression, color = Tumor_Normal1)) +
  geom_boxplot(outlier.size = 0.5,size = 1.6) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ Transcript_ID) +
  theme_bw() +
  theme(
    text = element_text(size = 25),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text.x = element_text(size = 24, colour = "black"),
    axis.text.y = element_text(size = 24, colour = "black"),
    strip.text = element_text(size = 22,colour = "black",face ="bold"),
    legend.text = element_blank(),
    legend.title = element_blank()) +
  labs(x = "Sample Type", y = "Transcript Expression")  +
  coord_cartesian(ylim = c(-0.6, 5)) +
  stat_compare_means(aes(group = Tumor_Normal1),
                     size=16,
                     label = "p.signif",
                     label.y = 4.1,label.x.npc = "middle") + stat_n_text(size = 6,color = "black",
                                                                         fontface = "bold")
g1
ggsave("HDAC9_expression_boxplot_pr1077.png", plot = g1, width = 8, height = 6, units = "in",
       bg = "white")

#############################################################################
#############################################################################

## For Figure 2i & j:------------>>

# Box plots for the HDAC9 pr1077 & pr1079 transcript expression in TCGA
# Read TCGA metadata:
tcga_metadata <- read_csv("C:/Users/Simran/Desktop/BRCA/tcga_data_new/MEGENA/Megena_on_TCGA/all_csv_files/Tum_plus_nor_nodupl_survivalinfo.csv")
tcga_metadata <- tcga_metadata[,-1]
tcga_metadata <- as.data.frame(tcga_metadata)

# Take just TNBC & NORMAL samples
keep <- c("YES","Normal")
tcga_metadata <- tcga_metadata[tcga_metadata$TNBC %in% keep,]

## Read HDAC9 pr1079 transcripts FPKM (log2 fpkm + 0.001)
# https://xenabrowser.net/datapages/?dataset=tcga_RSEM_isoform_fpkm&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
t1 <- read_tsv("ENST00000456174.6_xenaDownload.tsv")
t1 <- as.data.frame(t1)
t1_sel <- t1[t1$sample%in%tcga_metadata$trimmed_sample_submitter_id,]

t2 <- read_tsv("ENST00000524023.1_xenaDownload.tsv")
t2 <- as.data.frame(t2)
t2_sel <- t2[t2$sample%in%tcga_metadata$trimmed_sample_submitter_id,]

# Merge them
mm <- merge(t1_sel,t2_sel,by= "sample")
# Now merge with TCGA metadata
mm1 <- merge(mm,tcga_metadata,by.x = "sample",by.y = "trimmed_sample_submitter_id")

## Box plot:
# Just keeping required columns:
mm1 <- mm1[,c("sample","ENST00000456174.6","ENST00000524023.1","TNBC")]
# Pivot to long format
long_df <- mm1 %>%
  pivot_longer(
    cols = starts_with("ENST"),
    names_to = "Transcript_ID",
    values_to = "Expression"
  )

long_df$TNBC <- gsub("YES","TNBC",long_df$TNBC)
long_df$TNBC <- gsub("Normal","AdjNormal",long_df$TNBC)
my_colors <- c("TNBC"= "tomato2","AdjNormal"= "steelblue1")
my_colors2 <- c("TNBC"= "rosybrown1","AdjNormal"= "lightskyblue")

long_df$TNBC <- factor(long_df$TNBC,levels = c("TNBC","AdjNormal"))
long_df <- as.data.frame(long_df)

## Boxplot of expression per transcript, grouped by Tumor/Normal
g <- ggplot(long_df, aes(x = TNBC, y = Expression, color = TNBC)) +
  geom_boxplot(outlier.size = 0.5,size = 1.6) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ Transcript_ID) +
  theme_bw() +
  labs(x = "Sample Type", y = "Transcript Expression")  +
  theme(
    text = element_text(size = 25),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text.x = element_text(size = 24, colour = "black"),
    strip.text = element_text(size = 22,colour = "black",face ="bold"),
    axis.text.y = element_text(size = 24, colour = "black"),
    legend.text = element_blank(),
    legend.title = element_blank()
  ) +
  coord_cartesian(ylim = c(-11.5, 4.5)) +
  stat_compare_means(aes(group = TNBC),
                     size=12,
                     label = "p.signif",
                     label.y = 3,label.x.npc = "middle") + stat_n_text(size = 6,color = "black",
                                                                       fontface = "bold")
g
ggsave("TCGA_pr1079_transcript_expression_boxplot.png", plot = g, width = 8, height = 6, units = "in",
       bg = "white")

##--

# Now reading HDAC9 pr1077 transcript:
t3 <- read_tsv("ENST00000417496.6_xenaDownload.tsv")
t3 <- as.data.frame(t3)
t3_sel <- t3[t3$sample%in%tcga_metadata$trimmed_sample_submitter_id,]

# Merge with TCGA metadata
mm2 <- merge(t3_sel,tcga_metadata,by.x = "sample",by.y = "trimmed_sample_submitter_id")

## Box plot:
# Pivot to long format
long_df <- mm2 %>%
  pivot_longer(
    cols = starts_with("ENST"),
    names_to = "Transcript_ID",
    values_to = "Expression"
  )

long_df$TNBC <- gsub("YES","TNBC",long_df$TNBC)
long_df$TNBC <- gsub("Normal","AdjNormal",long_df$TNBC)
my_colors <- c("TNBC"= "tomato2","AdjNormal"= "steelblue1")

long_df$TNBC <- factor(long_df$TNBC,levels = c("TNBC","AdjNormal"))
long_df <- as.data.frame(long_df)

## Boxplot of expression per transcript, grouped by Tumor/Normal
g <- ggplot(long_df, aes(x = TNBC, y = Expression, color = TNBC)) +
  geom_boxplot(outlier.size = 0.5,size = 1.6) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ Transcript_ID) +
  theme_bw() +
  labs(x = "Sample Type", y = "Transcript Expression")  +
  theme(
    text = element_text(size = 25),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text.x = element_text(size = 24, colour = "black"),
    strip.text = element_text(size = 22,colour = "black",face ="bold"),
    axis.text.y = element_text(size = 24, colour = "black"),
    legend.text = element_blank(),
    legend.title = element_blank()
  ) +
  coord_cartesian(ylim = c(-11, 4.5)) +
  stat_compare_means(aes(group = TNBC),
                     size=12,
                     label = "p.signif",
                     label.y = 3,label.x.npc = "middle") + stat_n_text(size = 6,color = "black",
                                                                       fontface = "bold")
g
ggsave("TCGA_pr1077_transcript_expression_boxplot.png", plot = g, width = 5, height = 6, units = "in",
       bg = "white")

################################################################################################
################################################################################################

## For Figure 2k,l,m,n:------------>>

# KM PLOTS FOR HDAC9 PROMOTERS:----------->>
# Load the metadata
metadata <- as.data.frame(read_csv("meta_with_batchinfo.csv"))
# Data preprocessing
metadata <- metadata[, -1]  # Remove extra column
metadata <- metadata %>%
  mutate(Tumor_Normal1 = case_when(
    Tumor_Normal %in% c("TumorTissue", "TumorTissue.2018", "TumorTissue.rep") ~ "Tumor",
    Tumor_Normal %in% c("PairedNormalTissue", "PairedNormalTissue.rep") ~ "AdjNormal"
  ))
rownames(metadata) <- metadata$Run

# Take tumor samples
tumor_metadata <- subset(metadata, Tumor_Normal1 == 'Tumor')
tumor_metadata <- tumor_metadata[,c(1,2,3,85,87)]

## Read Overall Survival information:------------------------------>>
survival <- read_csv("survival_info_fuscc.csv")
colnames(survival) <- survival[1,]
survival <- survival[-1,]
survival <- survival[,c(1,43,45)]  # keeping: Project_ID   OS    OS_time_Months
survival <- as.data.frame(survival)

# Combine metadata
combined_meta <- merge(tumor_metadata,survival,by="Project_ID")

# For using Relapse-Free Survival (RFS)
test_surv_df <- data.frame(Run = tumor_metadata$Run, Event = tumor_metadata$RFS_Status,
                           Time = tumor_metadata$RFS_time_Months)

# Load result2 file
result2 <- readRDS('result2_v2_updated.rds')
result2_promoter = rowData(result2) 
filtered_df <- as.data.frame(result2_promoter)

# Not taking internal promoters' genes for downstream analysis
filtered_df <- filtered_df %>% 
  filter(internalPromoter == 'FALSE')

# Taking Absolute Promoter Activity
Ab_counts <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)
colnames(Ab_counts) <- gsub("_SJ.out", "", colnames(Ab_counts))
keep <-c("1077","1079")

Ab_counts <- Ab_counts[rownames(Ab_counts)%in%keep,]
Ab_counts <- t(Ab_counts)
colnames(Ab_counts) <- paste0("pr",colnames(Ab_counts))

test_surv_df <- as.data.frame(test_surv_df)
test_surv_df<-tibble::column_to_rownames(test_surv_df,"Run")

# Merge both dataframes
Ab_counts <- merge(test_surv_df, Ab_counts, by=0)
rownames(Ab_counts) <- Ab_counts$Row.names
Ab_counts <- Ab_counts[,-1]

Ab_counts <- tibble::rownames_to_column(Ab_counts,"Run")
merged_df <- Ab_counts

merged_df$Event <- as.numeric(merged_df$Event)
merged_df$Time <- as.numeric(merged_df$Time)

# Define promoter variables
prom <- c("pr1079", "pr1077")

for (pr in prom) {
  # Determine optimal cutoff using maxstat
  mod_maxstat <- maxstat.test(Surv(Time, Event) ~ merged_df[[pr]], 
                              data = merged_df, smethod = "LogRank")
  cutoff_value <- mod_maxstat$estimate
  
  # Binarize based on cutoff
  merged_df$binary <- ifelse(merged_df[[pr]] > cutoff_value, "high", "low")
  merged_df$binary <- as.factor(merged_df$binary)
  
  # Survival fit for RFS
  sfit <- survfit(Surv(Time, Event) ~ binary, data = merged_df)
  
  # Dynamically compute legend labels
  n_counts <- table(merged_df$binary)
  legend_labels <- paste0(toupper(names(n_counts)), " (n=", n_counts, ")")
  
  # Create KM plot for RFS
  p <- ggsurvplot(
    sfit,
    pval = TRUE,
    pval.size = 10,
    title = pr,
    legend.labs = legend_labels,
    legend.title = "",
    palette = c("#008080", "#FFA500"),
    risk.table = TRUE,
    xlim = c(0, 110),
    font.tickslab = 16,
    font.legend = 16,
    size = 1,
    risk.table.height = 0.25,
    risk.table.fontsize = 4.5,
    xlab = "Time (months)",
    ylab = "RFS"  
  )
  
  # Customize plot appearance
  p$plot <- p$plot + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title.x = element_text(face = "bold", size = 18),
      axis.title.y = element_text(face = "bold", size = 18),
      axis.text = element_text(face = "bold", size = 16),
      legend.position = c(0.7, 0.55),
      legend.background = element_rect(fill = "white", color = "white"),
      panel.border = element_rect(color = "black", fill = NA, size = 1.2),
      text = element_text(family = "Arial", face = 'bold')
    )
  
  # Print and save plot
  print(p)
  ggsave(paste0(pr, "_RFS.png"), p$plot, width = 4.4, height = 4, dpi = 600, bg = "white")
}

####----

# For using Overall Survival (OS)
test_surv_df <- data.frame(Run = combined_meta$Run, Event = combined_meta$OS,
                           Time = combined_meta$OS_time_Months)

# Load result2 file
result2 <- readRDS('result2_v2_updated.rds')
result2_promoter = rowData(result2) 
filtered_df <- as.data.frame(result2_promoter)

# Not taking internal promoters' genes for downstream analysis
filtered_df <- filtered_df %>% 
  filter(internalPromoter == 'FALSE')

# Taking Absolute Promoter Activity
Ab_counts <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)
colnames(Ab_counts) <- gsub("_SJ.out", "", colnames(Ab_counts))
keep <-c("1077","1079")

Ab_counts <- Ab_counts[rownames(Ab_counts)%in%keep,]
Ab_counts <- t(Ab_counts)
colnames(Ab_counts) <- paste0("pr",colnames(Ab_counts))

test_surv_df <- as.data.frame(test_surv_df)
test_surv_df<-tibble::column_to_rownames(test_surv_df,"Run")

# Merge both dataframes
Ab_counts <- merge(test_surv_df, Ab_counts, by=0)
rownames(Ab_counts) <- Ab_counts$Row.names
Ab_counts <- Ab_counts[,-1]

Ab_counts <- tibble::rownames_to_column(Ab_counts,"Run")
merged_df <- Ab_counts
merged_df$Event <- as.numeric(merged_df$Event)
merged_df$Time <- as.numeric(merged_df$Time)

# List of promoters to loop through
prom <- c("pr1079", "pr1077")

for (pr in prom) {
  # Determine optimal cutoff using maxstat
  mod_maxstat <- maxstat.test(Surv(Time, Event) ~ merged_df[[pr]], 
                              data = merged_df, smethod = "LogRank")
  cutoff_value <- mod_maxstat$estimate
  
  # Binarize based on cutoff
  merged_df$binary <- ifelse(merged_df[[pr]] > cutoff_value, "high", "low")
  merged_df$binary <- as.factor(merged_df$binary)
  
  # Survival fit
  sfit <- survfit(Surv(Time, Event) ~ binary, data = merged_df)
  
  # Dynamically compute the legend labels
  n_counts <- table(merged_df$binary)
  legend_labels <- paste0(toupper(names(n_counts)), " (n=", n_counts, ")")
  
  # Create KM plot
  p <- ggsurvplot(
    sfit,
    pval = TRUE,
    pval.size = 10,
    title = pr,
    legend.labs = legend_labels,
    legend.title = "",
    palette = c("#008080", "#FFA500"),
    risk.table = TRUE,
    xlim = c(0, 110),
    font.tickslab = 16,
    font.legend = 16,
    size = 1,
    risk.table.height = 0.25,
    risk.table.fontsize = 4.5,
    xlab = "Time (months)",
    ylab = "OS"
  )
  
  # Further customization
  p$plot <- p$plot + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title.x = element_text(face = "bold", size = 18),
      axis.title.y = element_text(face = "bold", size = 18),
      axis.text = element_text(face = "bold", size = 16),
      legend.position = c(0.7, 0.55),
      legend.background = element_rect(fill = "white", color = "white"),
      panel.border = element_rect(color = "black", fill = NA, size = 1.2),
      text = element_text(family = "Arial", face = 'bold')
    )
  
  # Print and save plot
  print(p)
  ggsave(paste0(pr, "_OS.png"), p$plot, width = 4.4, height = 4, dpi = 600, bg = "white")
}
#################################################################################
#################################################################################
