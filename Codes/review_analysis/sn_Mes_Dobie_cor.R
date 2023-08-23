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
################ Mouse - Dobie et al -  mesenchymal clusters ################
########################################################################

####################
Dobie_data = readRDS('Objects/mesenchymal/Dobie/Dobie_annotated.RDS')
Dobie_data = Dobie_data[rownames(Dobie_data) %in% VariableFeatures(Dobie_data),] 
cluster_names_types = names(table(Dobie_data$Manualanno_Sample))
cluster_names = Dobie_data$Manualanno_Sample

### dividing the expression matrix based on the clusters
cluster_expression <- sapply(1:length(cluster_names_types), function(i){
  a_cluster_name = cluster_names_types[i]
  GetAssayData(Dobie_data, 'data')[,cluster_names == a_cluster_name] 
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
cluster_average_exp_df$mouse_ID=rownames(cluster_average_exp_df)
head(cluster_average_exp_df)

rat_to_mouse_genes <- readRDS('~/XSpecies/rat_to_mouse_genes.rds')
rat_to_mouse_genes = rat_to_mouse_genes[rat_to_mouse_genes$mmusculus_homolog_orthology_type=='ortholog_one2one',]
rat_to_mouse_genes <- rat_to_mouse_genes[,c('mmusculus_homolog_associated_gene_name', 'symbol')]
colnames(rat_to_mouse_genes) = c('mouse_genes', 'rat_genes')

head(rat_to_mouse_genes)
head(cluster_average_exp_df)

mouse_data = merge(cluster_average_exp_df, rat_to_mouse_genes, by.x='mouse_ID', by.y='mouse_genes')
head(mouse_data)
dim(cluster_average_exp_df)
dim(mouse_data) ## more than 800 genes have been removed




#################################################################################################
#################### subcluster average expression of the mesenchymal subpopulation. #########
#################################################################################################
mes_data = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_mesenchymal_subclusters_merged_data_res2.5_cluster24Only.rds')

mes_data$cluster <- as.character(mes_data$SCT_snn_res.1)
cluster_names = mes_data$cluster 
cluster_names_types = names(table(cluster_names))

merged_samples = mes_data
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
colnames(cluster_average_exp_df) = paste0('cluster_',names(cluster_average_exp))
head(cluster_average_exp_df)

### QC-refined old samples (unified QC and MT removed)
## scale and center all the genes in the matrix
mes_cluster_average_exp <- get_scaled_by_gene(cluster_average_exp_df) ## scaling over the clusters of interest
mes_cluster_average_exp$rat_ID=rownames(mes_cluster_average_exp)
head(mes_cluster_average_exp)
dim(mes_cluster_average_exp)
#mes_cluster_average_exp <- mes_cluster_average_exp[mes_cluster_average_exp$rat_ID %in% rat_HVGs,]

rat_cluster_average.df = mes_cluster_average_exp
head(rat_cluster_average.df)
head(mouse_data)

rat_HVGs = VariableFeatures(mes_data)
############################################################
merge_mouse_rat = merge(mouse_data, rat_cluster_average.df, by.x='rat_genes', by.y='rat_ID')
head(merge_mouse_rat)
dim(merge_mouse_rat)
merge_mouse_rat_c = merge_mouse_rat[merge_mouse_rat$rat_genes %in% rat_HVGs,]
dim(merge_mouse_rat_c)
merge_mouse_rat_c = merge_mouse_rat_c[merge_mouse_rat_c$mouse_ID %in%  VariableFeatures(Dobie_data),]
dim(merge_mouse_rat_c)

head(merge_mouse_rat_c)
merge_mouse_rat_c = merge_mouse_rat_c[,-c(1, 2)] #10
dim(merge_mouse_rat_c)
head(merge_mouse_rat_c)
num_mouse_clusters = ncol(mouse_data) - 2
cor_mat = cor(merge_mouse_rat_c)[1:num_mouse_clusters, (num_mouse_clusters+1):ncol(merge_mouse_rat_c)]
colnames(cor_mat) = c('snMes.0', 'snMes.1', 'snMes.2')
rownames(cor_mat)
pheatmap::pheatmap(cor_mat, fontsize = 15) #main=paste0('#genes: ', nrow(merge_mouse_rat_c))





