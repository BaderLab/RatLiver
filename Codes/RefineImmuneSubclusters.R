source('Codes/Functions.R')
Initialize()
new_data_scCLustViz_object_Immune <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_ImmuneSub.RData"
load(new_data_scCLustViz_object_Immune)
old_immune_data = your_scRNAseq_data_object
old_immune_data$subclusters = paste0('old-', as.character(sCVdata_list$RNA_snn_res.1@Clusters))
  
new_data_scCLustViz_object_Immune <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_ImmuneSub_c17Included.RData"
load(new_data_scCLustViz_object_Immune)
new_immune_data = seur
new_immune_data$subclusters = as.character(new_immune_data$SCT_snn_res.1)

df_old <- data.frame(umi=colnames(old_immune_data), old_subclust=old_immune_data$subclusters) 
df_new <- data.frame(umi=colnames(new_immune_data), new_subclust=new_immune_data$subclusters, Cluster=new_immune_data$final_cluster)
merged_df = merge(df_new, df_old, by.x='umi', by.y='umi', all.x=T)
pheatmap(table(paste0('new-',merged_df$new_subclust), merged_df$old_subclust))
pheatmap(table(paste0('new-',merged_df$new_subclust), merged_df$Cluster))
