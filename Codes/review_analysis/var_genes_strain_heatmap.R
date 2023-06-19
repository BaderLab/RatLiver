##### select cluster nuumber 5 - marco+ cells of the single cell RNAseq data TLH
#### make a heatmap out of their expression values


############ set-1 map
old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)
cluster_df = data.frame(umi=colnames(your_scRNAseq_data_object),
                        clusters=as.character(sCVdata_list$res.0.6@Clusters))
set1_info <- read.csv('figure_panel/set-1-final-info.csv')
colnames(set1_info)[1] = 'clusters'
set1_info$clusters = as.character(set1_info$clusters)
cluster_df_set1 = merge(cluster_df, set1_info[,1:3], by.x='clusters', by.y='clusters', all.x=T,order=F)
#### re-ordering the rows based on UMIs of the initial data
cluster_df_set1 <- cluster_df_set1[match(colnames(your_scRNAseq_data_object),cluster_df_set1$umi),]
cluster_df_set1$umi == colnames(your_scRNAseq_data_object)
cluster_df_set1$sample_name = unlist(lapply(str_split(cluster_df_set1$umi, pattern = '_'), 
                                            function(x) paste0(x[-length(x)],collapse = '_')))
head(cluster_df_set1)
table(cluster_df_set1$sample_name)



merged_samples <- your_scRNAseq_data_object
merged_samples@meta.data = cbind(merged_samples@meta.data, cluster_df_set1)
table(cluster_df_set1$label)
macrophage_ids = cluster_df_set1$umi[cluster_df_set1$label %in% c('Non-Inflammatory Mac (5)', 'Non-Inflammatory Mac (10)')] #,

merged_samples$strain = ifelse(merged_samples$sample_name %in% c('rat_Lew_01', 'rat_Lew_02'), 'LEW', 'DA')
table(merged_samples$strain, merged_samples$sample_name)

merged_samples <- merged_samples[,colnames(merged_samples)%in%macrophage_ids]
dim(merged_samples)
sample_data = list(count=GetAssayData(merged_samples@assays$RNA))

strain_names_types = names(table(merged_samples$strain))
### dividing the expression matrix based on the clusters
strain_expression <- sapply(1:length(strain_names_types), function(i){
  a_strain_name = strain_names_types[i]
  GetAssayData(merged_samples, 'data')[,merged_samples$strain == a_strain_name] 
}, simplify = F)

names(strain_expression) = strain_names_types
lapply(strain_expression, head)


## calculate the average expression of each gene in each cluster
strain_average_exp <- lapply(strain_expression, function(x){
  df = data.frame(average=rowSums(x)/ncol(x))
  return(df)
})
lapply(strain_average_exp, dim)

## Concatenate all the clusters together to make a matrix
strain_average_exp_df = do.call(cbind,strain_average_exp)
colnames(strain_average_exp_df) = names(strain_average_exp)
head(strain_average_exp_df)
strain_average_exp_df$gene = rownames(strain_average_exp_df)





###################################################
##############   Not being used for now ###########
###################################################
rot_data <- readRDS('Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds') ## MT-removed
scores <- data.frame(rot_data$rotScores)
Loadings <- rot_data$rotLoadings
num_genes <- nrow(Loadings)
num_cells <- nrow(scores)
colnames(scores) = paste0('PC_', 1:ncol(scores))
scores$UMI <- rownames(scores)

merged_info = merge(re.df, strain_average_exp_df, by.x='gene', by.y='gene', sort=F)
merged_info_ord = merged_info[match(re.df_ord$gene, merged_info$gene),]
head(merged_info_ord)
merged_info_ord_sub = merged_info_ord[1:25,c("gene", 'DA', 'LEW')]
row.names(merged_info_ord_sub) <- merged_info_ord_sub$gene
###################################################


genes_to_check = c('Ly6al', 'Cd163', 'Hmox1', 'Siglec5', 'Itgal', 'Il18', 'Ccl3', "Timp2")

merged_info_sub = strain_average_exp_df[strain_average_exp_df$gene %in% genes_to_check,]
merged_info_ord_sub = strain_average_exp_df[genes_to_check,]


merged_info_ord_sub <- merged_info_ord_sub[,-3]
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks((as.matrix(merged_info_ord_sub)), n = 50)

pheatmap(merged_info_ord_sub, cluster_rows = F, cluster_cols = F, color=viridis(100, direction = -1), main='c5 and c10')

pheatmap(merged_info_ord_sub, cluster_rows = T, cluster_cols = F,  
         color = inferno(length(mat_breaks) - 1),
         breaks= mat_breaks)


pheatmap(merged_info_ord_sub, cluster_rows = T, cluster_cols = F,
         color=colorRampPalette(c("navy", "white", "red"))(50))
