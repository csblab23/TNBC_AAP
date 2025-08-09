############## This script includes all the codes used for FIGURE PANEL 1 ################

# Load required libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(webr)
library(dplyr)
library(proActiv)

## For Figure 1a:--------------->>

# Read the result2 file for all 448 FUSCC samples obtained after proActiv
result2 <- readRDS("result2_tum_vs_adjnormal.rds")
result2_promoter=rowData(result2)
result2_promoter=as.data.frame(result2_promoter)

# Not considering internal promoters for downstream analysis
result_prom <- result2_promoter[result2_promoter$internalPromoter == 'FALSE',]

# Check for single-promoter and multi-promoter genes
data_plot <- as.data.frame(table(unname(table(result_prom$geneId))))
colnames(data_plot) <- c('No_of_promoters', "Frequency")
data_plot$No_of_promoters <- as.numeric(data_plot$No_of_promoters) 
data_plot_modified2 <- data_plot %>%
  mutate(No_of_promoters = ifelse(No_of_promoters >= 4, 4, No_of_promoters))

# Obtain the total number of genes
total_genes <- sum(data_plot_modified2$Frequency)
data_plot_modified2$percent_of_genes <- (data_plot_modified2$Frequency / total_genes) * 100

# Compile multi-promoter genes
data <- data_plot_modified2[1:4,]
sum(data_plot_modified2$percent_of_genes[4:15]) 
data$percent_of_genes[4] <- sum(data_plot_modified2$percent_of_genes[4:15])
sum(data_plot_modified2$Frequency[4:15])  
data$Frequency[4] <- sum(data_plot_modified2$Frequency[4:15]) 

str(data)
data$No_of_promoters <- as.character(data$No_of_promoters)
data$No_of_promoters <- gsub("4",">=4",data$No_of_promoters)
data$No_of_promoters <- factor(data$No_of_promoters, levels = c('1', '2', '3', '>=4'))
data$final_freq <- round(data$percent_of_genes, 1)

# Calculate the difference to make the sum 100
diff <- 100 - sum(data$final_freq)
# Add the difference to the smallest positive value
data$final_freq[which.min(data$final_freq)] <- data$final_freq[which.min(data$final_freq)] + diff

## Take single promoters & multipromoters: just 2 categories to show on a pie donut plot:
data_sel <- as.data.frame(data[,c(1,4)])
data_sel$No_of_promoters <- if_else(data_sel$No_of_promoters=="1", "Single_promoter_genes","Multi-promoter_genes")
data_sel <- data_sel %>% 
  group_by(No_of_promoters) %>%
  summarize(total = sum(final_freq), .groups = 'drop')

# Plot:
plot <- ggplot(data_sel, aes(x = "", y = total, fill = No_of_promoters)) +
  geom_bar(width = 0.8, stat = "identity",color = "black",linewidth=0.6) +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_manual(values = c("Multi-promoter_genes" = "peru", "Single_promoter_genes" = "peachpuff2"),
                    labels = c("Multi promoter genes", "Single promoter genes")) +
  theme(legend.position = "bottom") +
  geom_text(aes(label = paste0(total,"%")), position = position_stack(vjust = 0.5),size = 12) +
  annotate("text", x = 0, y = 0, label = "37909 Genes\n47657 Promoters", size = 5.8, fontface = "bold") +
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_blank(),
    plot.margin = margin(5, 10, 10, 5),
    legend.direction = "vertical",
    legend.position = c(0.5, 0.01)
  )

print(plot)
ggsave("Fig1A_piedonut.png", plot, width = 6, height = 10, dpi = 1000)

#########################################################################
#########################################################################

## For Figure 1b:----------->>

# Load the result2 file for FUSCC:
result2 <- readRDS("result2_tum_vs_adjnormal.rds")
result2_promoter = rowData(result2)
result2_promoter = as.data.frame(result2_promoter)

# Keep non-internal promoters only:
result2_prom <- result2_promoter[result2_promoter$internalPromoter == 'FALSE', ]
# Check the distribution of the promoter types:
filtered_df <- result2_prom

# Remove inactive promoters from here
# we are considering only active promoters here
filtered_df_sel <- filtered_df[!filtered_df$Tumor.class=='Inactive',]
filtered_df_sel <- filtered_df_sel[,c('promoterId', 'geneId', 'Tumor.class')]

# Group by gene id and find number of genes with 1,2,3,>=4 promoters
gene_freq <- filtered_df_sel %>%
  group_by(geneId) %>%
  summarise(frequency = n())

df_with_freq <- filtered_df_sel %>%
  inner_join(gene_freq, by = "geneId")

df_with_freq <- df_with_freq %>%
  mutate(frequency_group = ifelse(frequency >= 4, ">=4", as.character(frequency))) %>%
  mutate(frequency_group = factor(frequency_group, levels = c("1", "2", "3", ">=4")))

# Identify the distribution of major/minor promoters for each category
df_prop <- df_with_freq %>%
  group_by(frequency_group, Tumor.class) %>%
  summarise(count = n()) %>%
  mutate(prop = (count / sum(count)) * 100)

# Stacked bar plot:
g<-ggplot(df_prop, aes(x = factor(frequency_group), y = prop, fill = Tumor.class)) +
  geom_bar(stat = "identity",color = "black",linewidth = 0.5) +
  geom_text(aes(label = paste0(round(prop, 1), "%")),
            position = position_stack(vjust = 0.5), color = "black", fontface = "bold", family = "Arial", size = 10) +
  scale_fill_manual(values = c("Major" = "cadetblue", "Minor" = '#CCC591')) +
  scale_y_continuous(limits = c(0,100),expand = c(0,0)) +
  labs(x = "No. of Promoters", y = "Proportion of Promoter Types") +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 28, family = "Arial", colour = 'black'),
    axis.title = element_text(size = 30, face = "bold", family = "Arial"),
    legend.text = element_text(size = 29, family = "Arial"),
    axis.text.x = element_text(angle = 90,size = 28 ,hjust = 1, vjust = 1, family = "Arial",colour = "black"),
    legend.title = element_blank(),
    legend.key.size = unit(1,"cm"),
    plot.margin = margin(10, 10, 10, 10),
  )
g
ggsave("fig1B_stackedbar.png", g, width = 9, height = 8, dpi = 1000)

######################################################################
######################################################################

## For Figure 1c:---->>

# Read the result2 file for FUSCC dataset:
result2 <- readRDS("result2_tum_vs_adjnormal.rds")
result2_promoter = rowData(result2)
result2_promoter = as.data.frame(result2_promoter)
result2_prom <- result2_promoter[result2_promoter$internalPromoter == 'FALSE', ]

# Keep genes that have more than or equal to two promoters (multi promoter genes):
filtered_df <- result2_prom %>%
  as.data.frame() %>%
  group_by(geneId) %>%
  filter(n() >= 2) 

pdata2 <- as_tibble(filtered_df) %>%
  mutate(promoterPosition = ifelse(promoterPosition > 5, 5, promoterPosition)) %>%
  filter(Tumor.class %in% c('Minor', 'Major'))

################### Chunk for knowing proportions ############################
table_output <- table(pdata2$geneId, pdata2$promoterPosition)
df_output <- as.data.frame.matrix(table_output)
colnames(df_output) <- c("position_1", "position_2", "position_3", "position_4", "position_5")
# Count the number of genes with promoters at each position
colSums(df_output)
prop.table(colSums(df_output))
sum(prop.table(colSums(df_output)))

####----

# Define the proportions
proportions <- c(0.46988645, 0.36179388, 0.10133678, 0.03564755, 0.03133535)
# Calculate the cumulative positions
positions <- cumsum(c(0, proportions))
# Summarize the data to get counts for each promoterPosition
pdata2_summary <- pdata2 %>%
  group_by(promoterPosition, Tumor.class) %>%
  summarize(count = n(), .groups = 'drop')

# Add widths and center positions
pdata2_summary <- pdata2_summary %>%
  mutate(
    width = proportions[as.integer(promoterPosition)],
    center = (positions[as.integer(promoterPosition)] + positions[as.integer(promoterPosition) + 1]) / 2
  )

pdata2_summary$Tumor.class = as.factor(pdata2_summary$Tumor.class)
pdata2_summary$Tumor.class <- relevel(pdata2_summary$Tumor.class, ref = "Major")

# Plot
gp<-ggplot(pdata2_summary) +
  geom_col(aes(x = center, y = count / sum(count), fill = Tumor.class, color = 'black'), 
           width = pdata2_summary$width, position = 'fill') +
  labs(x = expression(bold(Promoter ~ Position ~ "(5'" %->% " 3')"))) +
  ylab('Proportion of promoter types') + 
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = paste0(seq(0, 100, 25), '%'),expand = c(0,0)) +
  scale_x_continuous(breaks = (positions[-length(positions)] + positions[-1]) / 2, labels = c('1', '2', '3', '4', '>=5'),expand = c(0,0)) +
  scale_color_manual(values = c("black")) +
  scale_fill_manual(values = c("Minor" = "#CCC591", "Major" = "cadetblue")) +
  theme_classic() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 28, family = "Arial", colour = 'black'),
        axis.title = element_text(size = 30, face = "bold", family = "Arial", color = 'black'),
        legend.text = element_blank(),
        axis.text.x = element_text(size = 20, angle = 90, hjust = 1, colour = "black",family = "Arial", margin = margin(t = 5)),
        axis.title.y = element_text(size = 30, face = 'bold',colour = "black"),
        legend.title = element_blank(),
        plot.margin = unit(c(0.2, 0.1, 0.1, 0.1),"inches")) +
  
  guides(fill = guide_legend(title = "Promoter Category"))
gp
ggsave("fig1C_positionbarplot.png", gp, width = 8, height = 8, dpi = 1200)

#############################################################################
#############################################################################

## For Figure 1d:---->>

# Load result2 file for FUSCC
result2 <- readRDS("result2_tum_vs_adjnormal.rds")
result2_promoter = rowData(result2)
result2_promoter = as.data.frame(result2_promoter)

# Taking only genes that have internal promoter == False
result2_prom <- result2_promoter[result2_promoter$internalPromoter == 'FALSE', ]

# Select only the required columns:
df <- result2_prom[,c(1,2,13)]
df2 <- df %>% group_by(geneId)
# Just keep single promoter genes:
df2 <- df %>% group_by(geneId) %>% filter(n()<=1)

# Remove inactive promoters because we only want single prom genes with single active promoter:
df2 <-df2[!df2$Tumor.class=='Inactive',]
single_prom_genes <-df2

# Now for multipromoter genes:
df3 <- df %>% group_by(geneId) %>%
  filter(n()>=2)

# Segregate the single_active and multiple_active promoters from these:
df4 <- df3 %>%
  group_by(geneId) %>%
  mutate(prom_sum = sum(Tumor.class %in% c('Major','Minor')),
         
         class = case_when(prom_sum == 1 ~ 'single_active_pr',               
                           TRUE ~ NA_character_  # use NA_character_ for character columns
         )
  )

df4 <- df4 %>%
  mutate(
    class = case_when(
      prom_sum > 1 ~ 'multiple_active_pr',
      TRUE ~ class  # preserve existing values in the 'class' column
    )
  )

## Remove the columns with NA as sum (these have all inactive promoters)
df4_sel <- na.omit(df4)
df4_sel<-df4_sel[!df4_sel$Tumor.class=='Inactive',]

# Find out the unique single active and multiple active :
single_active_g <-df4_sel[df4_sel$class=='single_active_pr',]
multi_active_g<-df4_sel[df4_sel$class=='multiple_active_pr',]

## Stacked bar plot now:
df_prom <- data.frame(Prom_type = c('Multiple_active_pr(multi_pr_genes)', 'Single_active_pr(multi_pr_genes)',
                                    'Single_active_promoter(single_pr_genes)'),
                      Number = c(1998, 2488, 11788),
                      pr = c(1,1,1))

# Calculate the percentage:
total <- sum(df_prom$Number)
df_prom$per <- round(df_prom$Number/total *100,1)

# Convert to factor
df_prom$Prom_type <- factor(df_prom$Prom_type,levels = c('Multiple_active_pr(multi_pr_genes)',
                                                         'Single_active_pr(multi_pr_genes)',
                                                         'Single_active_promoter(single_pr_genes)'))

# Plot:
p<-ggplot(df_prom, aes(x = pr, y = per, fill = Prom_type)) +
  geom_bar(stat = "identity", position = "stack" ,width = 0.9,color="black") +
  geom_text(aes(label = paste0(round(per, 1), "%")),
            position = position_stack(vjust = 0.5),
            fontface = "bold",
            color = "black",
            family = "Arial",
            size = 8,
            angle =90) +
  scale_fill_manual(values = c("Multiple_active_pr(multi_pr_genes)" = "#e6c2c1", 
                               "Single_active_pr(multi_pr_genes)" = '#f2e4e4', 
                               "Single_active_promoter(single_pr_genes)" = "#966f6e"),
                    labels = c("Multiple Active Promoters\n(multi promoter genes)", 
                               "Single Active Promoters\n(multi promoter genes)", 
                               "Single Active Promoters\n(single promoter genes)")) +
  ylab('% of Active Promoters') + 
  labs(fill = 'Promoter Category') +
  scale_y_continuous(limits = c(0, 100),  
                     breaks = c(0,20,40,60,80,100),
                     labels = c("0","","","","","100"),
                     expand = expansion(mult = c(0, 0))) +  # remove extra space at the ends
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 30, family = "Arial", colour = 'black'),
    legend.text = element_text(size = 24, family = "Arial",color = "black",hjust = 0.5),
    axis.title.x =element_blank(),
    axis.text.y = element_text(angle = 90),
    axis.title.y = element_text(size = 30, face = 'bold',vjust = -2.5,hjust = 0.4),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 1),
    axis.line.x = element_blank(),
    legend.title = element_blank(),
    legend.spacing.y = unit(2.0,"cm"),
    legend.key.size = unit(2.5, "cm"),
    plot.margin = margin(20, 10, 20, 10))  # Add margin on the left for y-axis label
p
ggsave("Fig_categoriesofpromoters3.png", p, width = 8, height = 8, dpi = 1000)

##--

## For the line plot:
# Load result2 file for FUSCC:
result2 <- readRDS("result2_tum_vs_adjnormal.rds")

# Load the absolute promoter activity:
abs_pa <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)  
colnames(abs_pa) <- gsub("_SJ.out","",colnames(abs_pa))
abs_pa <- tibble::rownames_to_column(abs_pa,var = "promoterId")

# Read promoter categories file:
promcategoriesdf <- as.data.frame(read_csv('prom_categories.csv'))
promcategoriesdf<-promcategoriesdf[,-1]

# Inactive are not required here:
promcategoriesdf2<-promcategoriesdf[!promcategoriesdf$Tumor.class=='Inactive',]
# Merge with result2_prom_sel df:
merged_df <- merge(promcategoriesdf2, abs_pa, by="promoterId")

# Specify the columns to be used:
cols_sel <- colnames(merged_df[5:ncol(merged_df)])

# Apply summarise with across so that we get output as expected
merged_df2 <- merged_df %>% 
  group_by(geneId) %>%
  summarise(across(all_of(cols_sel), sum, .names = "{col}_total"), .groups = 'drop')

# Calculate the rowmeans
gene_expr_df <- merged_df2 %>%
  mutate(avg_gene_expr = rowMeans(.[, 2:ncol(.)], na.rm = TRUE)) %>%
  relocate(avg_gene_expr)

# Just keep required columns:
gene_expr_df_sel <- gene_expr_df[,c('avg_gene_expr','geneId')]

# Calculating average major Promoter Activity:
merged_df_major <- merged_df[merged_df$Tumor.class=="Major",]

# Calculate the Average expression of major promoter across all the samples:
merged_df_major_sel <- merged_df_major %>%
  mutate(avg_major_PA = rowMeans(.[, 5:ncol(.)], na.rm = TRUE)) %>%
  relocate(avg_major_PA)

# Subset & keep the required columns:
major_pa_df <- merged_df_major_sel[,c('avg_major_PA','promoterId','geneId','Tumor.class','type')]

####----

# Merge the two dfs:
combined <- merge(gene_expr_df_sel, major_pa_df, by='geneId',all.x = T, all.y =T)
# Final merge to compile all information
combined_final <- merge(combined, promcategoriesdf2, by='geneId')
combined_final <- combined_final[,-c(4:6)]
colnames(combined_final)[5] <-'Tumor_class'
colnames(combined_final)[6] <-'gene_type'
combined_final$gene_type <- factor(combined_final$gene_type,
                                   levels = c("single_promoter_gene",
                                              "multi_promoter_gene_multiple_active_pr",
                                              "multi_promoter_gene_single_active_pr"))

# Plot:
g<-ggplot(data = combined_final, aes(x = avg_gene_expr, y = avg_major_PA, color = gene_type)) +
  geom_point(size=2,alpha = 0.9) +                           
  scale_color_manual(values = c("single_promoter_gene"="#966f6e",
                                "multi_promoter_gene_multiple_active_pr" = "#e6c2c1",
                                "multi_promoter_gene_single_active_pr" = "#f2e4e4"),
                     labels = c("Single Promoter Genes","Multiple Active Promoters\n(multi promoter genes)", "Single Active Promoters\n(multi promoter genes)"))+
  annotate("segment",x = 0, y = 0, xend = 30, yend = 30, 
           linetype = "dashed", color = "gold3", size = 1) +
  scale_x_continuous(limits = c(0,30)) +
  labs(x = "Average Gene Expression",     
       y = "Average Major Promoter Activity") +
  theme_classic() +
  theme(
    axis.title.x = element_text(colour = "black", size = 23,face = "bold"),
    axis.text.x = element_text(colour = "black", size = 22),
    axis.title.y = element_text(colour = "black", size = 23, face = "bold"),
    axis.text.y = element_text(colour = "black", size = 22),
    legend.text = element_text(colour = "black", size = 20),
    legend.title = element_blank(),
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "right",
    legend.spacing.x = unit(3.0,"cm"),
    legend.key.spacing.y = unit(0.6,"cm"),
    legend.box.background = element_rect(color = "black", linewidth = 1),
    panel.border = element_rect(fill = NA, color = "black", size = 1)
  )

g
ggsave("Fig_average_gexpr_pa_updated3.png", g, width = 10, height = 7, dpi = 500)

#####################################################################################
#####################################################################################

## For figure 1e ---->>

# We will use the differentially regulated promoters (DRPs) obtained from Absolute Promoter Activity:
# Load the upregulated DRPs
up_drp <- read_csv("upreg_drp_absPA_bypvalue.csv")
up_drp<- up_drp[,-1]

# Load the downregulated DRPs
down_drp <- read_csv("down_drp_absPA_bypvalue.csv")
down_drp<-down_drp[,-1]

# Combine all drps
all_drps <- rbind(up_drp, down_drp)

# Rank the Promoters by log2FC normalized promoter activity
all_drps_ordered <- all_drps[order(all_drps$log2foldChange),]
to_be_removed <- c(Inf,-Inf) # remove infinite fold changes
all_drps_ordered <- all_drps_ordered[!all_drps_ordered$log2foldChange %in% to_be_removed,]

## Now we need gene names for these promoter ids:
# Fetch the differentially expressed genes (DEG) results file
deg_all <- read_csv("DEG_wilcoxon_tum_vs_nor.csv")
colnames(deg_all)[1]<- "geneId"

# Important to fetch only the non-significant genes:
deg_nonsig <- deg_all[deg_all$pValues>0.05,] 

# Load filtered_df file wherein --> internal promoters are absent; and genes with >=2 promoters are present
# This file was prepared from the result2 object
filtered_df <- read_csv("filtered_df.csv")
filtered_df<-filtered_df[,-1]
filtered_df_subset <- filtered_df[,c('promoterId','geneId')]

# Merge to have both prom id and gene id:
df1 <- merge(deg_nonsig, filtered_df_subset, by="geneId")

# Move back to the promoter df:
subset_all_Drps_ordered <- as.data.frame(all_drps_ordered[,1])  # just taking prom id column
df2 <- merge(df1,subset_all_Drps_ordered, by.x='promoterId', by.y='prom_id')

# Make sure to remove the genes with foldchnge between -1 & 1 only:
df2_sel <- df2[df2$log2foldChange<1 & df2$log2foldChange>-1,]

# Keep just these promoters only in the above ordered dataframe:
all_drps_ordered_sel <- all_drps_ordered[all_drps_ordered$prom_id %in% df2_sel$promoterId,]
colnames(all_drps_ordered_sel)[2]<- 'promoter_log2FC'

# Making a subset of just required columns:
all_drps_ordered_sel <- all_drps_ordered_sel[,c('prom_id','promoter_log2FC')]

# Adjusting the column names
colnames(df2_sel)[3]<- 'gene_log2FC'
df2_sel <- df2_sel[,c('promoterId', 'geneId', 'gene_log2FC')]

# Final merge and order by log2FC
merged_final <- merge(all_drps_ordered_sel, df2_sel, by.x='prom_id', by.y = 'promoterId')
merged_final <- merged_final[order(merged_final$promoter_log2FC),]  # order by prom log2fc

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

# Subset required columns
subset_meta <- metadata[,c('Run', 'batch','Tumor_Normal1')]   # let us just keep the required columns

# Read the result2 file:
result2 <- readRDS("tum_nor_alternate_prom_result/result2_tum_vs_adjnormal.rds")

# Load Absolute promoter activity
abs_pa <- as.data.frame(result2@assays@data@listData$absolutePromoterActivity)  
colnames(abs_pa) <- gsub("_SJ.out","",colnames(abs_pa))

# we have to do filtering now so let us fetch the rowData from result2:
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

# Keep the same promoter ids in the absolute counts matrix
abs_pa_filtered <- abs_pa[rownames(abs_pa) %in% filtered_df$promoterId,]
count_norm <- abs_pa_filtered

# Heatmap
t_sel=t(count_norm)  # transpose so rows are samples and columns correspond to genes
ti=merge(subset_meta,t_sel,by=0)
ti_all=ti[order(ti$Tumor_Normal1),,drop=F]
rownames(ti_all)=ti_all$Row.names
ti_all$Tumor_Normal1 <- gsub('AdjNormal','NORMAL', ti_all$Tumor_Normal1)
ti_all$Tumor_Normal1 <- gsub('Tumor',"TNBC",ti_all$Tumor_Normal1)
meta1=ti_all[,1:4]
exp=ti_all[,5:dim(ti_all)[2]]

## keep the above filtered promoters only here:
exp_sel <- exp[,colnames(exp) %in% merged_final$prom_id]
exp_sel <- as.matrix(exp_sel)
merged_final$prom_id <- as.character(merged_final$prom_id)
exp_sel<- exp_sel[,merged_final$prom_id]  # making order of promoter ids same in both files
keep <- colSums(exp_sel>1) >=60  # filtering here so that heatmap appears better as there are many 0 zero values
exp_sel2 <- exp_sel[,keep]

library(ComplexHeatmap)
library("circlize")

png("heatmap_fig.png", width = 1200, height = 1200,res = 100 )
meta1=meta1[order(meta1$Tumor_Normal1),]
exp_sel2=exp_sel2[rownames(meta1),]
meta11=meta1[,c(3:4)]

# Assign appropriate colours
cell_colors <- c("NORMAL"= "steelblue1","TNBC"= "tomato2")
col_fun = colorRamp2(c(-1.5,-0.7,0,0.7,1.5), c("bisque2","bisque", "papayawhip","salmon1","salmon3"))

# Convert to factor
meta11$Tumor_Normal1=as.factor(meta11$Tumor_Normal1)
h1 <- Heatmap(meta11$Tumor_Normal1,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              width = unit(2, "cm"),
              col = cell_colors,
              show_heatmap_legend = FALSE
)

h2=Heatmap(scale(exp_sel2),name = " ", cluster_rows = TRUE, cluster_columns = FALSE,
           show_row_names = FALSE, show_column_names = FALSE,
           col = col_fun,heatmap_legend_param = list(title_gp = gpar(fontsize = 10,fontface = 'bold'), 
                                                     labels_gp = gpar(fontsize = 30,fontface = "bold")))

h1+h2
dev.off()
combined_heatmap <- h1+h2

# Extract the order of column names
column_order <- combined_heatmap@ht_list[[" "]]@column_names_param[["labels"]]

##########

# Keep same promoter ids
merged_final_new <- merged_final[merged_final$prom_id %in% column_order,]
rownames(merged_final_new) <- merged_final_new$prom_id
merged_final_new <- merged_final_new[column_order,]

# Create a line plot of log2foldchange for promoter
g1 <- ggplot(merged_final_new, aes(x = seq_along(prom_id), y = promoter_log2FC)) +
  geom_line(color="midnightblue",linewidth = 1.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  
  labs(x = "Ranked Promoters", y = "Promoter Activity \n(log2FC)") +          
  theme_minimal() + 
  theme(axis.text = element_text(size = 21,face = 'bold',color = 'black'),
        axis.title.x = element_text(size = 22,face = 'bold'),
        axis.title.y = element_text(size = 22,face = 'bold'))

# Create a scatter plot of log2foldchange for gene
g2 <- ggplot(merged_final_new, aes(x = seq_along(geneId), y = gene_log2FC)) +
  geom_point(color="seagreen",size=3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  
  labs(x = "Ranked Promoters", y = "Gene Expression \n(log2FC)") +        
  theme_minimal() + 
  theme(axis.text = element_text(size = 21,face = 'bold',color = 'black'),
        axis.title.x = element_text(size = 22,face = 'bold'),
        axis.title.y = element_text(size = 22,face = 'bold')) + scale_y_continuous(limits = c(-1,1))

# combine the 2 plots vertically:
library(gridExtra)
combined_plots <- grid.arrange(g1, g2, nrow = 2)
ggsave("Fig_heatmap_bottomplots.png", combined_plots, width = 10, height = 7, dpi = 500)

###############################################################
###############################################################

## For Figure 1f & g ---->>

# 1. Venn between upregulated DEG & upregulated DRP genes (by pvalue<0.05):
up_deg <- as.data.frame(read_csv("upreg_2443deg_bypvalue.csv"))
up_deg<-up_deg[,-1]

up_drp <- as.data.frame(read_csv("upreg_drp_absPA_bypvalue.csv"))
up_drp<-up_drp[,-1]
colnames(up_drp)[1]<-'promoterId'

## We have to add the gene information here in up_drp df:
# Read the result2 file:
result2 <- readRDS("./tum_nor_alternate_prom_result/result2_tum_vs_adjnormal.rds")

# We have to do filtering now so fetch the rowData from result2:
result2_promoter = rowData(result2) 
filtered_df <- as.data.frame(result2_promoter)
filtered_df <- filtered_df[, c('promoterId','geneId','internalPromoter')]

# Filter 1: keep internal promoters = FALSE ------------------------------
table(filtered_df$internalPromoter)
filtered_df <- filtered_df[filtered_df$internalPromoter == "FALSE",]

# Filter 2: keep the genes which have >=2 prom ---------------------------
# for this we will fetch the rowData
filtered_df <- filtered_df %>%
  group_by(geneId) %>%
  filter(n() >= 2) %>%
  mutate(ID = paste('pr', promoterId, '_', geneId, sep='')) 

merged_df <- merge(up_drp,filtered_df,by='promoterId')
dim(merged_df)
length(unique(merged_df$geneId))

# Now find common genes:
common_genes <- intersect(up_deg$...2, merged_df$geneId)
write.csv(common_genes,"common_48genes_UP_pvalue.csv")

# Plot venn diagram
# Specify values to use in venn diagram
library(eulerr)
up_degss <- unique(up_deg$...2) 
up_drpgs <- unique(merged_df$geneId) 

venn <- list(UP_DEGs = up_degss, UP_DRPGs = up_drpgs)
euler_data <- euler(venn)

# Create proportionate venn diagram:
p1<-plot(euler_data,  fills = list(fill = c("burlywood4","salmon3"), alpha = 0.7),
         labels = list(cex=2.2,hjust=0.5), quantities = list(cex = 2,hjust=0.5),cex = 2.5)
p1
ggsave("Fig_venn_upreg.png", p1, width = 8, height = 7, dpi = 1000)

#######

# 2. Venn between downregulated DEGs & downregulated DRP genes:
down_deg <- as.data.frame(read_csv("down_2039deg_bypvalue.csv"))
down_deg<-down_deg[,-1]

down_drp <- as.data.frame(read_csv("down_drp_absPA_bypvalue.csv"))
down_drp<-down_drp[,-1]
colnames(down_drp)[1] <- 'promoterId'

## Add the gene id information here in up_drp df:
# Read the result2 file:
result2 <- readRDS("tum_nor_alternate_prom_result/result2_tum_vs_adjnormal.rds")

# Apply the filtering
result2_promoter = rowData(result2) 
filtered_df <- as.data.frame(result2_promoter)
filtered_df <- filtered_df[, c('promoterId','geneId','internalPromoter')]

# Filter 1: keep internal promoters = FALSE ------------------------------
table(filtered_df$internalPromoter)
filtered_df <- filtered_df[filtered_df$internalPromoter == "FALSE",]

# Filter 2: keep the genes which have >=2 prom ---------------------------
# for this we will fetch the rowData
filtered_df <- filtered_df %>%
  group_by(geneId) %>%
  filter(n() >= 2) %>%
  mutate(ID = paste('pr', promoterId, '_', geneId, sep='')) 

merged_df <- merge(down_drp,filtered_df,by='promoterId')

# Now find common genes:
common_genes <- intersect(down_deg$...2, merged_df$geneId)
write.csv(common_genes,"common_151genes_DOWN_pvalue.csv")

# Specify the values for venn:
library(eulerr)
down_degss <- unique(down_deg$...2)
down_drpgs <- unique(merged_df$geneId)
venn <- list(DOWN_DEGs = down_degss, DOWN_DRPGs = down_drpgs)
euler_data <- euler(venn)

# Create proportionate venn diagram:
p2<-plot(euler_data,  fills = list(fill = c("darkseagreen4","bisque2"), alpha = 0.7),
         labels = list(cex=2.2,hjust=0.5), quantities = list(cex = 2.5,hjust=0.5),cex = 2.5)
p2
ggsave("Fig_venn_downreg_feb5.png", p2, width = 8, height = 7, dpi = 1000)

#############################################################################
#############################################################################

## For Figure 1h & i ---->>

# Read all the gmt files:
gmt_file1 <- upload_GMT_file("c6.all.v2024.1.Hs.symbols (1).gmt") 
gmt_file2 <- upload_GMT_file("c2.all.v2024.1.Hs.symbols.gmt") 
gmt_file3 <- upload_GMT_file("c2.reactome.gmt") 
gmt_file4 <- upload_GMT_file("c5.biological.process.gmt") 
gmt_file5 <- upload_GMT_file("c5.cellular.component.gmt") 
gmt_file6 <- upload_GMT_file("c5.molecular.function.gmt") 
gmt_file7 <- upload_GMT_file("hallmark.gmt") 
gmt_file8 <- upload_GMT_file("c7.all.v2024.1.Hs.symbols.gmt") 
gmt_file9 <- upload_GMT_file("c2.biocarta.gmt")
gmt_file10 <- upload_GMT_file("c2.kegg_legacy.gmt")  

# Make a list of all gmt files:
gmtfiles <- list(gmt_file1, gmt_file2, gmt_file3, gmt_file4, gmt_file5, 
                 gmt_file6, gmt_file7, gmt_file8, gmt_file9, gmt_file10)

## 1. ANALYSIS FOR UPREGULATED GENES:------------------------->>
# between 2443 upregulated DEG & 692 upregulated DRP genes:
up_deg <- as.data.frame(read_csv("upreg_2443deg_bypvalue.csv"))
up_deg<-up_deg[,-1]

up_drp <- as.data.frame(read_csv("upreg_drp_absPA_bypvalue.csv"))
up_drp<-up_drp[,-1]
colnames(up_drp)[1]<-'promoterId'

## Add the gene information here in up_drp df:
# Read the result2 file:
result2 <- readRDS("result2_tum_vs_adjnormal.rds")

# Filtering
result2_promoter = rowData(result2) 
filtered_df <- as.data.frame(result2_promoter)
filtered_df <- filtered_df[, c('promoterId','geneId','internalPromoter')]

# Filter 1: keep internal promoters = FALSE ------------------------------
table(filtered_df$internalPromoter)
filtered_df <- filtered_df[filtered_df$internalPromoter == "FALSE",]

# Filter 2: keep the genes which have >=2 prom ---------------------------
# for this we will fetch the rowData
filtered_df <- filtered_df %>%
  group_by(geneId) %>%
  filter(n() >= 2) %>%
  mutate(ID = paste('pr', promoterId, '_', geneId, sep='')) 

merged_df <- merge(up_drp,filtered_df,by='promoterId')

# Now find common genes:
common_genes <- intersect(up_deg$...2, merged_df$geneId)
write.csv(common_genes,"common_48genes_UP_pvalue.csv")

# Removing these common genes from deg list:
deg_2395 <- up_deg[!up_deg$...2%in%common_genes,]
deg_2395 <- deg_2395$...2
write.csv(deg_2395,'deg_up_2395.csv')

# Removing these common genes from DRPG list:
drpg_664 <- merged_df[!merged_df$geneId%in%common_genes,]
drpg_664<-drpg_664$geneId
drpg_664<-unique(drpg_664)
write.csv(drpg_664,'drpg_up_664.csv')

##--

# Read the gene lists of interest ------DEG
all_gene_ids <- read_csv("deg_up_2395.csv")
all_gene_ids <- all_gene_ids$x

# Remove digits after dots
all_gene_ids_2 <- gsub("\\..*","",all_gene_ids)

# Convert the ensembl ids to gene names:
genes <- gprofiler2::gconvert(all_gene_ids_2,organism = 'hsapiens', target = 'ENSG')
gene_names=genes[,c(2,5)]  # fetching just the id and gene name columns

up_deg <- gene_names$name
up_deg <- unique(up_deg)

# Similarly read the second gene list ----- DRPGs
all_gene_ids <- read_csv("drpg_up_664.csv")
all_gene_ids <- all_gene_ids$x

# Remove digits after dots:
all_gene_ids_2 <- gsub("\\..*","",all_gene_ids)

# Convert ensembl ids to gene names:
genes <- gprofiler2::gconvert(all_gene_ids_2,organism = 'hsapiens', target = 'ENSG')
gene_names=genes[,c(2,5)]  # fetching just the id and gene name columns

up_drpg <- gene_names$name
up_drpg<-unique(up_drpg)

# Read the third input list ------- overlapping genes
all_gene_ids <- read_csv("common_48genes_UP_pvalue.csv")
all_gene_ids <- all_gene_ids$x

# Remove digits after dots:
all_gene_ids_2 <- gsub("\\..*","",all_gene_ids)

# Convert ensembl ids to gene names:
genes <- gprofiler2::gconvert(all_gene_ids_2,organism = 'hsapiens', target = 'ENSG')
gene_names=genes[,c(2,5)]  # fetching just the id and gene name columns
head(gene_names)
up_common <- gene_names$name
input_lists <- list(DEGs = up_deg, DRPGs = up_drpg, COMMON = up_common)
str(input_lists)

## Perform pathway analysis using gost in a loop:------->>
# Initialize an empty dataframe to store results
result <- data.frame()

for (i in seq_along(gmtfiles)){
  
  gostres <- gost(query = input_lists,  organism = gmtfiles[[i]][[1]],
                  ordered_query = FALSE, multi_query = FALSE, significant = TRUE,
                  exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE,
                  user_threshold = 0.05, correction_method = "fdr", domain_scope = "annotated",
                  custom_bg = NULL, numeric_ns = "", sources = NULL, as_short_link = FALSE)
  
  # Extract the results:
  p_c2 <- gostres$result
  table(p_c2$query)
  
  if ("query" %in% colnames(p_c2)) {
    p_c2_updated <- p_c2 %>%
      group_by(query) %>%
      arrange(p_value) %>%
      slice_head(n = 10)  # fetching top 10 pathways obtained
    
    # Append the current result to the main dataframe
    result <- bind_rows(result, p_c2_updated)
  }
}

# Save the results:
result <- as.data.frame(result)
table(result$query)
result$query <- factor(result$query, levels = c("DEGs", "DRPGs", "COMMON"))

# Sort the data frame based on the custom order
result_sorted <- result[order(result$query), ]
rownames(result_sorted) <- NULL
result_sorted <- as.data.frame(result_sorted)

# Save the results:
writexl::write_xlsx(result_sorted,"upreg_pathway_result.xlsx")

# Remove duplicates and keep unique pathways:
upreg <- as.data.frame(read_excel('upreg_pathway_result.xlsx'))

# verifying that there are no duplicate pathways for each query:
upreg2 <- upreg %>%
  group_by(query) %>%
  distinct(term_id,.keep_all = TRUE) %>%
  ungroup()
table(upreg2$query)

## Fetch the unique pathways:
com_up <- upreg2[upreg2$query=='COMMON',]
deg_up <- upreg2[upreg2$query=='DEGs',]
drpg_up <- upreg2[upreg2$query=='DRPGs',]

##-----

# remove these overlapping pathways inorder to fetch out the unique ones:
overlap_deg_drpgs <- intersect(deg_up$term_id,drpg_up$term_id)
overlap_deg_common <- intersect(deg_up$term_id,com_up$term_id)
overlap_drpg_common <- intersect(drpg_up$term_id,com_up$term_id)
remove_pathways <- c(overlap_deg_drpgs,overlap_deg_common,overlap_drpg_common)

# Creating separate dataframes unique to each category
uniq_common <- com_up[!com_up$term_id %in% remove_pathways,]
uniq_common_sorted <- uniq_common[order(uniq_common$p_value),]  # sort by pvalue
writexl::write_xlsx(uniq_common_sorted,'unique_common_UP.xlsx')

uniq_degs <- deg_up[!deg_up$term_id %in% remove_pathways,]
uniq_degs_sorted <- uniq_degs[order(uniq_degs$p_value),]
writexl::write_xlsx(uniq_degs_sorted,'unique_degs_UP.xlsx')

uniq_drpgs <- drpg_up[!drpg_up$term_id %in% remove_pathways,]
uniq_drpgs_sorted <- uniq_drpgs[order(uniq_drpgs$p_value),]
writexl::write_xlsx(uniq_drpgs_sorted,'unique_drpgs_UP.xlsx')


# Bubbleplots:
combined_df <- rbind(uniq_common_sorted, uniq_degs_sorted, uniq_drpgs_sorted)
writexl::write_xlsx(combined_df,"FINAL_upreg_pathway_result.xlsx")

####----

# Just keep the selected relevant pathways:
pathways_sel <- c('REACTOME_ACTIVATED_PKN1_STIMULATES_TRANSCRIPTION_OF_AR_ANDROGEN_RECEPTOR_REGULATED_GENES_KLK2_AND_KLK3',
                  'REACTOME_RHO_GTPASE_EFFECTORS', 'REACTOME_HDACS_DEACETYLATE_HISTONES',
                  'REACTOME_SIRT1_NEGATIVELY_REGULATES_RRNA_EXPRESSION', 
                  'SMID_BREAST_CANCER_BASAL_UP','REACTOME_DISEASES_OF_DNA_REPAIR',
                  'GOBP_DNA_RECOMBINATION','GOBP_MITOTIC_CELL_CYCLE',
                  'GOBP_REGULATION_OF_SMALL_GTPASE_MEDIATED_SIGNAL_TRANSDUCTION',
                  'GOMF_UBIQUITIN_LIKE_PROTEIN_LIGASE_BINDING',
                  'GOMF_G_PROTEIN_COUPLED_CHEMOATTRACTANT_RECEPTOR_ACTIVITY')

combined_df_sel <- combined_df[combined_df$term_id %in% pathways_sel,]
combined_df_sel<-as.data.frame(combined_df_sel)

# Cleaning the pathway names
combined_df_sel$term_id <- gsub("REACTOME_","",combined_df_sel$term_id)
combined_df_sel$term_id <- gsub("GOBP_","",combined_df_sel$term_id)
combined_df_sel$term_id <- gsub("GOMF_","",combined_df_sel$term_id)
combined_df_sel$term_id <- gsub("HALLMARK_","",combined_df_sel$term_id)
combined_df_sel$term_id <- gsub("ACTIVATED_PKN1_STIMULATES_TRANSCRIPTION_OF_AR_ANDROGEN_RECEPTOR_REGULATED_GENES_KLK2_AND_KLK3",
                                "ACTIVATED_PKN1_STIMULATES\nANDROGEN_RECEPTOR_REGULATED\nGENES_KLK2_AND_KLK3",combined_df_sel$term_id)
combined_df_sel$term_id <- gsub("SIRT1_NEGATIVELY_REGULATES_RRNA_EXPRESSION",
                                "SIRT1_NEGATIVELY_REGULATES\nRRNA_EXPRESSION",combined_df_sel$term_id)
combined_df_sel$term_id <- gsub("G_PROTEIN_COUPLED_CHEMOATTRACTANT_RECEPTOR_ACTIVITY",
                                "G_PROTEIN_COUPLED_CHEMOATTRACTANT\nRECEPTOR_ACTIVITY",combined_df_sel$term_id)
combined_df_sel$term_id <- gsub("REGULATION_OF_SMALL_GTPASE_MEDIATED_SIGNAL_TRANSDUCTION",
                                "REGULATION_OF_SMALL_GTPASE\nMEDIATED_SIGNAL_TRANSDUCTION",combined_df_sel$term_id)
combined_df_sel$term_id <- gsub("UBIQUITIN_LIKE_PROTEIN_LIGASE_BINDING",
                                "UBIQUITIN_LIKE_PROTEIN\nLIGASE_BINDING",combined_df_sel$term_id)

# Ensure that term_id is a factor with levels ordered by query_mapped
combined_df_sel$term_id <- factor(combined_df_sel$term_id,
                                  levels = unique(combined_df_sel$term_id[order(combined_df_sel$query)]))

# Calculate -log10 of pvalue:
combined_df_sel$neglog10_pvalue <- -log10(combined_df_sel$p_value)

# Bar plot for pathways obtained:
p <- ggplot(combined_df_sel, aes(x = term_id, y = neglog10_pvalue, fill = query)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("COMMON" = "peachpuff3", "DEGs" = "burlywood4", "DRPGs"="salmon3"),
                    breaks = c("DRPGs", "DEGs", "COMMON")) +
  labs(x = "Pathway", 
       y = "-log10(Pvalue)") +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 1,colour = "black",size = 34),
        axis.text.y = element_text(colour = "black",size = 33),
        axis.line = element_line(colour = "black"),
        axis.ticks.y = element_line(colour = "black"),
        axis.title.x = element_text(colour = "black",face = "bold",size = 38),
        axis.title.y = element_text(colour = "black",face = "bold",size = 38),
        legend.title = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
        legend.text = element_text(colour = "black",size = 35,face = "bold"),
        panel.grid.major = element_line(color = "honeydew3", linewidth = 0.8, linetype = "dashed"),
        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.8, linetype = "dashed"),
        legend.spacing.y = unit(2, "cm"))

p
ggsave("Selected_upreg_pathways_FINAL.png", p, width = 22, height = 18, dpi = 1000)

#######################################################################################
#######################################################################################

## 2. ANALYSIS FOR DOWNREGULATED GENES:------------------------->>
# between 2024 downregulated DEG & 1531 downregulated DRP genes:
down_deg <- as.data.frame(read_csv("down_2039deg_bypvalue.csv"))
down_deg<-down_deg[,-1]

down_drp <- as.data.frame(read_csv("down_drp_absPA_bypvalue.csv"))
down_drp<-down_drp[,-1]
colnames(down_drp)[1] <- 'promoterId'

## Add the gene information here in up_drp df:
# Read the result2 file:
result2 <- readRDS("result2_tum_vs_adjnormal.rds")

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

merged_df <- merge(down_drp,filtered_df,by='promoterId')

# Now find common genes:
common_genes <- intersect(down_deg$...2, merged_df$geneId)
write.csv(common_genes,"common_151genes_DOWN_pvalue.csv")

# Remove common genes from deg dataframe:
deg_1888 <- down_deg[!down_deg$...2%in%common_genes,]
deg_1888 <- deg_1888$...2
write.csv(deg_1888,'deg_down_1888.csv')

# Remove common genes from drpg dataframe:
drpg_1225 <- merged_df[!merged_df$geneId%in%common_genes,]
drpg_1225<-drpg_1225$geneId
drpg_1225<-unique(drpg_1225)
write.csv(drpg_1225,'drpg_down_1225.csv')

##--

# Read the gene lists of interest ------DEG
all_gene_ids <- read_csv("deg_down_1888.csv")
all_gene_ids <- all_gene_ids$x

# Remove digits after dots
all_gene_ids_2 <- gsub("\\..*","",all_gene_ids)

# Convert the ensembl ids to gene names:
genes <- gprofiler2::gconvert(all_gene_ids_2,organism = 'hsapiens', target = 'ENSG')
gene_names=genes[,c(2,5)]  # fetching just the id and gene name columns
down_deg <- gene_names$name
down_deg <- unique(down_deg)

# Similarly read the second gene list ----- DRPGs
all_gene_ids <- read_csv("drpg_down_1225.csv")
all_gene_ids <- all_gene_ids$x

# Remove digits after dots:
all_gene_ids_2 <- gsub("\\..*","",all_gene_ids)

# Convert ensembl ids to gene names:
genes <- gprofiler2::gconvert(all_gene_ids_2,organism = 'hsapiens', target = 'ENSG')
gene_names=genes[,c(2,5)]  # fetching just the id and gene name columns
down_drpg <- gene_names$name
down_drpg<-unique(down_drpg)

# Read the third input list ------- overlapping genes
all_gene_ids <- read_csv("common_151genes_DOWN_pvalue.csv")
all_gene_ids <- all_gene_ids$x
head(all_gene_ids)

# Remove digits after dots:
all_gene_ids_2 <- gsub("\\..*","",all_gene_ids)

# Convert ensembl ids to gene names:
genes <- gprofiler2::gconvert(all_gene_ids_2,organism = 'hsapiens', target = 'ENSG')
gene_names=genes[,c(2,5)]  # fetching just the id and gene name columns
head(gene_names)
down_common <- gene_names$name
down_common <- unique(down_common)

input_lists <- list(DEGs = down_deg, DRPGs = down_drpg, COMMON = down_common)
str(input_lists)

## Perform pathway analysis using gost in a loop:------------>>
# Initialize an empty dataframe to store results
result <- data.frame()
for (i in seq_along(gmtfiles)){
  
  gostres <- gost(query = input_lists,  organism = gmtfiles[[i]][[1]],
                  ordered_query = FALSE, multi_query = FALSE, significant = TRUE,
                  exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE,
                  user_threshold = 0.05, correction_method = "fdr", domain_scope = "annotated",
                  custom_bg = NULL, numeric_ns = "", sources = NULL, as_short_link = FALSE)
  
  # Extract the results:
  p_c2 <- gostres$result
  table(p_c2$query)
  
  if ("query" %in% colnames(p_c2)) {
    p_c2_updated <- p_c2 %>%
      group_by(query) %>%
      arrange(p_value) %>%
      slice_head(n = 10)  # fetching top 10 pathways obtained
    
    # Append the current result to the main dataframe
    result <- bind_rows(result, p_c2_updated)
  }
}

# Save the results:
result <- as.data.frame(result)
result$query <- factor(result$query, levels = c("DEGs", "DRPGs", "COMMON"))

# Sort the data frame based on the custom order
result_sorted <- result[order(result$query), ]
rownames(result_sorted) <- NULL
result_sorted <- as.data.frame(result_sorted)

# Save the results:
writexl::write_xlsx(result_sorted,"downreg_pathway_result.xlsx")

####

# Remove duplicate pathways and keep unique:
downreg <- as.data.frame(read_excel('downreg_pathway_result.xlsx'))

# Verifying that there are no duplicate pathways for each query:
downreg2 <- downreg %>%
  group_by(query) %>%
  distinct(term_id,.keep_all = TRUE) %>%
  ungroup()

## Fetch the unique pathways:
com_down <- downreg2[downreg2$query=='COMMON',]
deg_down <- downreg2[downreg2$query=='DEGs',]
drpg_down <- downreg2[downreg2$query=='DRPGs',]

##-----

overlap_deg_drpgs <- intersect(deg_down$term_id,drpg_down$term_id)
overlap_deg_common <- intersect(deg_down$term_id,com_down$term_id)
overlap_drpg_common <- intersect(drpg_down$term_id,com_down$term_id)

# Remove these overlapping pathways inorder to fetch out the unique ones:
remove_pathways <- c(overlap_deg_drpgs,overlap_deg_common,overlap_drpg_common)

# Creating separate dataframes for each category
uniq_common <- com_down[!com_down$term_id %in% remove_pathways,]
uniq_common_sorted <- uniq_common[order(uniq_common$p_value),]  # sort by pvalue
writexl::write_xlsx(uniq_common_sorted,'unique_common_DOWN.xlsx')

uniq_degs <- deg_down[!deg_down$term_id %in% remove_pathways,]
uniq_degs_sorted <- uniq_degs[order(uniq_degs$p_value),]
writexl::write_xlsx(uniq_degs_sorted,'unique_degs_DOWN.xlsx')

uniq_drpgs <- drpg_down[!drpg_down$term_id %in% remove_pathways,]
uniq_drpgs_sorted <- uniq_drpgs[order(uniq_drpgs$p_value),]
writexl::write_xlsx(uniq_drpgs_sorted,'unique_drpgs_DOWN.xlsx')

# Bubbleplots for them:
combined_df <- rbind(uniq_common_sorted, uniq_degs_sorted, uniq_drpgs_sorted)
writexl::write_xlsx(combined_df,"FINAL_downreg_pathway_result.xlsx")

####----

# Just keep the selected pathways:
pathways_sel <- c('REACTOME_GPCR_LIGAND_BINDING', 'KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION',
                  'HALLMARK_KRAS_SIGNALING_DN','REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION','GOBP_CELL_JUNCTION_ORGANIZATION',
                  'LIM_MAMMARY_STEM_CELL_UP','GOBP_CELL_MORPHOGENESIS','GOMF_ACTIN_BINDING','KEGG_PATHWAYS_IN_CANCER',
                  'GOBP_POSTSYNAPTIC_MEMBRANE_ASSEMBLY','GOBP_NEURON_DEVELOPMENT','REACTOME_ACTIVATION_OF_TRKA_RECEPTORS')

combined_df_sel <- combined_df[combined_df$term_id %in% pathways_sel,]

# Cleaning the pathway name column 
combined_df_sel$term_id <- gsub("REACTOME_","",combined_df_sel$term_id)
combined_df_sel$term_id <- gsub("GOBP_","",combined_df_sel$term_id)
combined_df_sel$term_id <- gsub("GOMF_","",combined_df_sel$term_id)
combined_df_sel$term_id <- gsub("HALLMARK_","",combined_df_sel$term_id)
combined_df_sel$term_id <- gsub("KEGG_","",combined_df_sel$term_id)
combined_df_sel$term_id <- gsub("ACTIVATION_OF_TRKA_RECEPTORS",
                                "ACTIVATION_OF_RECEPTOR\nTYROSINE_KINASES",combined_df_sel$term_id)
combined_df_sel$term_id <- gsub("CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
                                "CYTOKINE_CYTOKINE\nRECEPTOR_INTERACTION",combined_df_sel$term_id)
combined_df_sel$term_id <- gsub("EXTRACELLULAR_MATRIX_ORGANIZATION",
                                "EXTRACELLULAR_MATRIX\nORGANIZATION",combined_df_sel$term_id )
combined_df_sel$term_id <- gsub("POSTSYNAPTIC_MEMBRANE_ASSEMBLY",
                                "POSTSYNAPTIC_MEMBRANE\nASSEMBLY",combined_df_sel$term_id)

# Ensure that term_id is a factor with levels ordered by query_mapped
combined_df_sel$term_id <- factor(combined_df_sel$term_id,
                                  levels = unique(combined_df_sel$term_id[order(combined_df_sel$query)]))

# Calculate -log10 of pvalue:
combined_df_sel$neglog10_pvalue <- -log10(combined_df_sel$p_value)

# Bar plot for the obtained pathways:
p <- ggplot(combined_df_sel, aes(x = term_id, y = neglog10_pvalue, fill = query)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("COMMON" = "#ABB090", "DEGs" = "darkseagreen4", "DRPGs"="bisque2"),
                    breaks = c("DRPGs", "DEGs", "COMMON")) +
  labs(x = "Pathway", 
       y = "-log10(Pvalue)") +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 1,colour = "black",size = 34),
        axis.text.y = element_text(colour = "black",size = 33),
        axis.line = element_line(colour = "black"),
        axis.ticks.y = element_line(colour = "black"),
        axis.title.x = element_text(colour = "black",face = "bold",size = 38),
        axis.title.y = element_text(colour = "black",face = "bold",size = 38),
        legend.title = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
        legend.text = element_text(colour = "black",size = 35,face = "bold"),
        panel.grid.major = element_line(color = "honeydew3", linewidth = 0.8, linetype = "dashed"),
        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.8, linetype = "dashed"),
        legend.spacing.y = unit(2, "cm"))

p
ggsave("Selected_downreg_pathways_FINAL.png", p, width = 19, height = 18, dpi = 1000)

#####################################################################################
#####################################################################################
