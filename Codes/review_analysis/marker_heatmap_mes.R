source('Codes/Functions.R')
source('~/XSpecies/Codes/convert_human_to_ortholog_functions.R')
library(Hmisc)
library(RColorBrewer)
library(viridis)
library(reshape2)


################################
##### Generating average expression based heatmap

get_scaled_by_gene <- function(x){
  y <- scale(t(x), center = TRUE, scale = TRUE)
  y[is.nan(y)] <- 0
  y2 <- t(y)
  y2 <- as.data.frame(y2)
  return(y2)
}


########################################################
############ Importing the nuc-seq data  #########

## importing the gene expression data
merged_samples <- readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_mesenchymal_subclusters_merged_data_res2.5_cluster24Only.rds')
head(merged_samples)

merged_samples$cluster <-  as.character(merged_samples$SCT_snn_res.1)
table(merged_samples$cluster)

# gene(FB-associated)
Set1 = c("Dpt", "Serpinf1", "Gfh", "Vkorc1", "Cgln1", "Mgst1", "Clu", "Rpsa", "Ftl1", "Igfbp4",
         "Mgst1", "Rpl11", "Rps15", "Rps2", "Rps23", "Rps5", "Rpsa", "Rpl24")
# gene (HSC-associated): 
#Set2 = c("Colec11", "Col4a1", "Prelp", "Ecm1", "Colec10", "Hgf", "Lrat", "Pth1r", "Calcrl", "Rspo3", "Tgfbi", "Notch1", "Igfbp7", "BGN")
Set2 = c("Reln", "Pth1r", "Rxra","Vipr1","Colec11", "Col4a1", "Col3a1", "Prelp", "Ecm1", "Colec10", "Hgf", 
  "Lrat", "Pth1r", "Calcrl", "Rspo3", "Tgfbi", "Notch1", "Igfbp7", "Bgn")
# gene (VSMC-associated):
Set3 =  c("Acta2", "Gadd45g", "Id2", "Errfil", "Rasd1", "H3f3a", "Net1")

markers_list = c(Set1,Set2, Set3) 
markers_list[!markers_list %in% row.names(merged_samples)] #"Gfh"    "Cgln1"  "BGN"    "Errfil" not inncluded
# "Gfh"    "Cgln1"  "Rxra"   "Vipr1"  "Errfil"

markers_list = markers_list[markers_list %in% row.names(merged_samples)]

#merged_samples_subset2 <- CreateAssayObject(GetAssayData(merged_samples_subset)[rownames(merged_samples_subset) %in% markers_list,])
merged_samples_subset <- merged_samples[rownames(merged_samples) %in% markers_list,]

dim(merged_samples_subset)
#merged_samples_subset2$manual_annotation = merged_samples_subset$manual_annotation
cluster_names = as.character(merged_samples_subset$cluster )
cluster_names_types = names(table(cluster_names))

### dividing the expression matrix based on the clusters
cluster_expression <- sapply(1:length(cluster_names_types), function(i){
  a_cluster_name = cluster_names_types[i]
  GetAssayData(merged_samples_subset, 'data')[,cluster_names == a_cluster_name] 
}, simplify = F)

names(cluster_expression) = cluster_names_types
lapply(cluster_expression, head)


## calculate the average expression of each gene in each cluster
cluster_average_exp <- lapply(cluster_expression, function(x){
  df = data.frame(average=rowSums(x)/ncol(x))
  return(df)
})
lapply(cluster_average_exp, dim)

## Concatenate all the clusters together to make a matrix
cluster_average_exp_df = do.call(cbind,cluster_average_exp)
colnames(cluster_average_exp_df) = paste0(names(cluster_average_exp))
head(cluster_average_exp_df)


## scale and center all the genes in the matrix
cluster_average_exp_df_scaled <- get_scaled_by_gene(cluster_average_exp_df) ## scaling over the clusters of interest
colnames(cluster_average_exp_df_scaled)
rownames(cluster_average_exp_df_scaled)
cluster_average_exp_df_scaled = cluster_average_exp_df_scaled[markers_list,]

pheatmap(t(cluster_average_exp_df_scaled), cluster_rows = FALSE,
         cluster_cols = FALSE, 
         #annotation_row = data.frame(anno_df),
         #annotation_col = data.frame(anno_df_genes),
         color= inferno(20, direction = +1),  
         border_color= FALSE, 
         annotation_names_row=FALSE,
         annotation_names_col=FALSE,
         #annotation_colors = list(CellType=anno_df_col_list),
         annotation_legend = TRUE
)
  