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
stellate_data = readRDS('Objects/mesenchymal/Stellate_annotated_V3.RDS')
stellate_seur = GetAssayData(stellate_data)

stellate_data <- FindVariableFeatures(stellate_data)
human_HVGs <- VariableFeatures(stellate_data)
stellate_data = stellate_data[rownames(stellate_data) %in% human_HVGs,] 
cluster_names_types = names(table(stellate_data$sub_annotations_new))
cluster_names = stellate_data$sub_annotations_new

### dividing the expression matrix based on the clusters
cluster_expression <- sapply(1:length(cluster_names_types), function(i){
  a_cluster_name = cluster_names_types[i]
  GetAssayData(stellate_data, 'data')[,cluster_names == a_cluster_name] 
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

########################################################################
################ Rat - Set1 map -  mesenchymal clusters ################
########################################################################

rat_cluster_average.df <- readRDS('Results/rat_new_cluster_average_exp_all.rds')
rat_cluster_average.df <- readRDS('Results/rat_old_cluster_average_exp_all.rds')

colnames_set2 = c("pDC (17)", "Naive T cell (10)", "Erythroid (5)", "Hep (0)" , "Hep (1)", "Hep (14)",
                  "Hep (15)", "Hep (2)", "Hep (3)", "Hep (9)", "Inflammatory Mac (11)", "LSEC (4)",
                  "LSEC (6)", "Non-Inflammatory Mac (8)", "Mature B cell (12)", "gd T cell (7)", 
                  "Non-Inflammatory Mac (13)", "Stellate (16)", "rat_ID" )
colnames(rat_cluster_average.df) = colnames_set2

mes_sub_seur = readRDS('Results/old_samples/Mesenchymal_subclusters.rds')
mes_sub_seur <- FindVariableFeatures(mes_sub_seur)
rat_mes_HVGs <- VariableFeatures(mes_sub_seur)
mes_sub_seur = mes_sub_seur[rownames(mes_sub_seur) %in% rat_mes_HVGs,] 


umap_df = getUmapDF(mes_sub_seur)
umap_df$gcluster = mes_sub_seur$final_cluster
umap_df$res.0.2 = mes_sub_seur$SCT_snn_res.0.2
umap_df$res.0.4 = mes_sub_seur$SCT_snn_res.0.4
ggplot(umap_df,aes(UMAP_1, UMAP_2, color=gcluster))+geom_point()+theme_bw()

mes_sub_seur$mes_subclust = paste0('rat_mes_', mes_sub_seur$SCT_snn_res.0.4)

cluster_names_types = names(table(mes_sub_seur$mes_subclust ))
cluster_names = mes_sub_seur$mes_subclust

### dividing the expression matrix based on the clusters
cluster_expression <- sapply(1:length(cluster_names_types), function(i){
  a_cluster_name = cluster_names_types[i]
  GetAssayData(mes_sub_seur, 'data')[,cluster_names == a_cluster_name] 
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
head(cluster_average_exp_df)
cluster_average_exp_df$rat_ID = rownames(cluster_average_exp_df)
rat_cluster_average.df = cluster_average_exp_dfs
head(rat_cluster_average.df)

############################################################
merge_human_rat = merge(human_data, rat_cluster_average.df, by.x='rat_genes', by.y='rat_ID')
head(merge_human_rat)
dim(merge_human_rat)

merge_human_rat = merge_human_rat[-c(1, 2)] #10
pheatmap::pheatmap(cor(merge_human_rat)[1:7, 8:ncol(merge_human_rat)])

############################################################
merge_mouse_rat = merge(mouse_data, rat_cluster_average.df, by.x='rat_genes', by.y='rat_ID')
head(merge_mouse_rat)
dim(merge_mouse_rat)

merge_mouse_rat = merge_mouse_rat[-c(1, 2)] #10
pheatmap::pheatmap(cor(merge_mouse_rat)[1:3, 4:ncol(merge_mouse_rat)])


