source('Codes/Functions.R')
Initialize()

get_scaled_by_gene <- function(x){
  y <- scale(t(x), center = TRUE, scale = TRUE)
  y[is.nan(y)] <- 0
  y2 <- t(y)
  y2 <- as.data.frame(y2)
  return(y2)
}


########################################################################
################. Human - Andrews et al. mesenchymal clusters ################
########################################################################

# integrated_data = readRDS('Objects/mesenchymal/Integrated_with_Subannotations.rds')
hep_data = readRDS('~/rat_sham_sn_data/Tallulah_Paper_RDS/hep_harmony_annotated.RDS')

hep_data <- FindVariableFeatures(hep_data)
#human_HVGs <- VariableFeatures(hep_data)
#hep_data = hep_data[rownames(hep_data) %in% human_HVGs,] 
head(hep_data)
table(hep_data$sub_annotations)
table(hep_data$All_Integrated_Manual)
cluster_names_types = names(table(hep_data$sub_annotations))
cluster_names = hep_data$sub_annotations

### dividing the expression matrix based on the clusters
cluster_expression <- sapply(1:length(cluster_names_types), function(i){
  a_cluster_name = cluster_names_types[i]
  GetAssayData(hep_data, 'data')[,cluster_names == a_cluster_name] 
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
#colnames(cluster_average_exp_df) = paste0('cluster_',names(cluster_average_exp))
colnames(cluster_average_exp_df) = names(cluster_average_exp)
head(cluster_average_exp_df)

## scale and center all the genes in the matrix
cluster_average_exp_df <- get_scaled_by_gene(cluster_average_exp_df) ## scaling over the clusters of interest
cluster_average_exp_df$human_ID=rownames(cluster_average_exp_df)
head(cluster_average_exp_df)

rat_to_human_genes <- readRDS('~/RatLiver/Results/rat_to_human_genes_all.rds')
rat_to_human_genes = rat_to_human_genes[rat_to_human_genes$hsapiens_homolog_orthology_type=='ortholog_one2one',]
rat_to_human_genes <- rat_to_human_genes[,c('hsapiens_homolog_associated_gene_name', 'symbol')]
colnames(rat_to_human_genes) = c('human_genes', 'rat_genes')

head(rat_to_human_genes)
head(cluster_average_exp_df)

human_data = merge(cluster_average_exp_df, rat_to_human_genes, by.x='human_ID', by.y='human_genes')
head(human_data)
dim(human_data)

########################################################################
################ Rat - sinle nuc seq data ################
########################################################################
## importing the gene expression data
merged_samples = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')

########################################################################
############### Analysis based on resolution 0.6 ##############
Resolution = 0.6
merged_samples <- FindClusters(merged_samples, resolution = Resolution, verbose = FALSE)
table(merged_samples$SCT_snn_res.0.6)

rat_HVGs <- VariableFeatures(merged_samples)
head(merged_samples)
table(merged_samples$annot_IM_g,merged_samples$SCT_snn_res.0.6)
table(merged_samples$SCT_snn_res.0.6[merged_samples$annot_IM_g == 'Hep'])
nucseq_heps = as.character(c(0:4, 6, 9:11))

########################################################################
############### Analysis based on resolution 2.5 ##############

Resolution = 2.5
merged_samples <- FindClusters(merged_samples, resolution = Resolution, verbose = FALSE)
table(merged_samples$SCT_snn_res.2.5)
rat_HVGs <- VariableFeatures(merged_samples2)


#### defining the mega clusters based on 0.6 resolution
Hep0 = as.character(c(0, 4, 7, 10, 16, 25, 27)) # 27 is the little tail
Hep1 = as.character(c(1, 6, 15, 17, 12, 31)) # 31 is the tail
Hep2 = as.character(c(2, 5, 9, 21)) #
Hep3 = as.character(c(3, 8,13, 23, 32 ))


merged_samples$clusters = as.character(merged_samples$SCT_snn_res.2.5)
merged_samples$Hep_clusters = ifelse(merged_samples$clusters %in% as.character(Hep0), 'Hep0', merged_samples$clusters) 
merged_samples$Hep_clusters = ifelse(merged_samples$Hep_clusters %in% as.character(Hep1), 'Hep1', merged_samples$Hep_clusters) 
merged_samples$Hep_clusters = ifelse(merged_samples$Hep_clusters %in% as.character(Hep2), 'Hep2', merged_samples$Hep_clusters) 
merged_samples$Hep_clusters = ifelse(merged_samples$Hep_clusters %in% as.character(Hep3), 'Hep3', merged_samples$Hep_clusters) 
table(merged_samples$Hep_clusters)

merged_samples = merged_samples[, merged_samples$Hep_clusters %in% c('Hep0', 'Hep1', 'Hep2', 'Hep3')]

####################################################



##### annotate the hepatocytes based on the mouse zonation layers

merged_samples$cluster <- as.character(merged_samples$Hep_clusters)
cluster_names = merged_samples$cluster 
cluster_names_types = names(table(cluster_names))

### dividing the expression matrix based on the clusters
cluster_expression <- sapply(1:length(cluster_names_types), function(i){
  a_cluster_name = cluster_names_types[i]
  GetAssayData(merged_samples, 'data')[,cluster_names == a_cluster_name] 
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
colnames(cluster_average_exp_df) = paste0(names(cluster_average_exp)) #'cluster_'
head(cluster_average_exp_df)

### QC-refined old samples (unified QC and MT removed)
clusters_to_check <- paste0('cluster_',nucseq_heps)
Hep_cluster_average_exp <- get_scaled_by_gene(cluster_average_exp_df[,colnames(cluster_average_exp_df) %in% clusters_to_check]) ## scaling over the clusters of interest

## scale and center all the genes in the matrix
Hep_cluster_average_exp <- get_scaled_by_gene(cluster_average_exp_df)
Hep_cluster_average_exp$rat_ID=rownames(Hep_cluster_average_exp)
head(Hep_cluster_average_exp)
dim(Hep_cluster_average_exp)
Hep_cluster_average_exp <- Hep_cluster_average_exp[Hep_cluster_average_exp$rat_ID %in% rat_HVGs,]

rat_cluster_average.df = Hep_cluster_average_exp
head(rat_cluster_average.df)
head(human_data)

############################################################
merge_human_rat = merge(human_data, rat_cluster_average.df, by.x='rat_genes', by.y='rat_ID')
head(merge_human_rat)
dim(merge_human_rat)



merge_human_rat = merge_human_rat[-c(1, 2)] #10
num_human_clusters = ncol(human_data) - 2
pheatmap::pheatmap(cor(merge_human_rat)[1:num_human_clusters, (num_human_clusters+1):ncol(merge_human_rat)])

############################################################
merge_mouse_rat = merge(mouse_data, rat_cluster_average.df, by.x='rat_genes', by.y='rat_ID')
head(merge_mouse_rat)
dim(merge_mouse_rat)

merge_mouse_rat = merge_mouse_rat[-c(1, 2)] #10
pheatmap::pheatmap(cor(merge_mouse_rat)[1:3, 4:ncol(merge_mouse_rat)])








