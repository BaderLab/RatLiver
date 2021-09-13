source('Codes/Functions.R')
Initialize()

old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)
set1_data <- your_scRNAseq_data_object
set1_data@meta.data$cluster = as.character(sCVdata_list$res.0.6@Clusters)
set1_data$strain = unlist(lapply(str_split(set1_data$orig.ident, '_'), '[[', 2))

##### filtering based on highly variable genes
set1_data <- FindVariableFeatures(set1_data, nfeatures=500)
set1_HVGs <- VariableFeatures(set1_data,)
set1_data = set1_data[set1_HVGs,]
######### calculating the cluster average expression values

strain_names_types = names(table(set1_data$strain))
### dividing the expression matrix based on the clusters
strain_expression <- sapply(1:length(strain_names_types), function(i){
  a_strain_name = strain_names_types[i]
  set1_data[, set1_data$strain== a_strain_name]
}, simplify = F)

names(strain_expression) = strain_names_types
lapply(strain_expression, dim)

strain_average_exp <- lapply(strain_expression, function(x){
  ## calculate the average expression of each gene in each cluster
  df = data.frame(average=rowSums(x)/ncol(x))
  return(df)
})
## Concatenate all the clusters together to make a matrix
strain_average_df = do.call(cbind,strain_average_exp)
colnames(strain_average_df) = names(strain_average_exp)
head(strain_average_df)
strain_average_df$rat_ID = rownames(strain_average_df)


num_PCs = 35
#### loading varimax factors
varimax_res_set1 <- readRDS('Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds')
score_set1 <- data.frame(varimax_res_set1$rotScores)
colnames(score_set1) = paste0('Varimax_', 1:ncol(score_set1))
loading_set1 <- data.frame(varimax_res_set1$rotLoadings[,1:ncol(varimax_res_set1$rotLoadings)])
loading_set1$genes = rownames(loading_set1)

merged_old = merge(strain_average_df,loading_set1[,c(1:num_PCs,ncol(loading_set1))], 
                   by.x='rat_ID', by.y='genes')
dim(merged_old)
head(merged_old)
merged_old_cor = cor(merged_old[,-1])
strain_num = 2
pheatmap(merged_old_cor[1:strain_num, (strain_num+1):ncol(merged_old_cor)],
         cluster_cols = F, 
         main = paste0('set1-HVGs strain-varimax correlation')) #2000HGV



varimax_res_set2 <- readRDS('Results/new_samples/varimax_rotated_object_new.rds')
score_set2 <- data.frame(varimax_res_set2$rotScores)
colnames(score_set2) = paste0('Varimax_', 1:ncol(score_set2))
loading_set2 <- data.frame(varimax_res_set2$rotLoadings[,1:ncol(varimax_res_set2$rotLoadings)])
loading_set2$genes = rownames(loading_set2)


merged_new = merge(rat_new_cluster_average_exp_df,loading_set2[,c(1:num_PCs,ncol(loading_set2))], by.x='rat_ID', by.y='genes')
dim(merged_new)
merged_new_cor = cor(scale(merged_new[,-1]))
cluster_num = ncol(rat_new_cluster_average_exp_df) -1
pheatmap(merged_new_cor[1:cluster_num, (cluster_num+1):ncol(merged_new_cor)],
         cluster_cols = F, 
         main = paste0('set2-Top ',num_DE_genes,' marker cluster-varimax correlation'))





