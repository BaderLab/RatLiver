source('Codes/Functions.R')
#source('Codes/FactorAnalysisUtils.R')
Initialize()
library(plyr)
library(Polychrome)


mes_sub_seur = readRDS('Results/old_samples/Mesenchymal_subclusters.rds')

cluster_names <- names(table(mes_sub_seur$SCT_snn_res.0.4))
merged_samples = mes_sub_seur
### finding the markers while removing the mito-genes
Cluster_markers <- sapply(1:length(cluster_names), 
                          function(i) FindMarkers(merged_samples, 
                                                  ident.1=cluster_names[i],
                                                  logfc.threshold = 0,
                                                  min.pct=0,
                                                  min.cells.feature = 1,
                                                  min.cells.group = 1
                          ), 
                          simplify = FALSE)
names(Cluster_markers) <- cluster_names


P_value_thr = 0.05
Cluster_markers_final <- sapply(1:length(Cluster_markers), function(i) {
  
  ## selecting the cluster of interest's marker dataframe (x)
  x = Cluster_markers[[i]]
  a_cluster_name <- names(Cluster_markers)[i]
  
  ## sort rows of x based on log-fold change
  x = x[order(x$avg_log2FC, decreasing = T),]
  
  ## sort based on adj p-value and sign of logFC  
  # x$ranking_score=-log10(x$p_val_adj+.Machine$double.xmin)*sign(x$avg_log2FC)
  # x = x[order(x$ranking_score, decreasing = T),]
  
  ## add the average expression of each gene as a column
  selected_cells = Idents(merged_samples) == a_cluster_name
  data <- GetAssayData(merged_samples)[rownames(x),]
  x$avg_exp <- rowSums(data[,selected_cells])/sum(selected_cells)
  
  ## filtering genes with adj-p-value higher than 0.05
  #x = x[x$p_val_adj<P_value_thr,]
  
  return(x)
}, simplify = F)

names(Cluster_markers_final) <- names(Cluster_markers)
lapply(Cluster_markers_final, head)
#lapply(Cluster_markers_final, dim)

################################################################################
Resolution = 0.4
saveRDS(Cluster_markers_final, paste0('Results/old_samples/TLH_mesenchymal_markers_res',Resolution,'.rds'))

lapply(Cluster_markers_final, head)
#### saving markers #####
marker_dir = paste0('Results/old_samples/TLH_mesenchymal_markers_res',Resolution,'/')
dir.create(marker_dir)

for(i in 1:length(Cluster_markers_final)){
  #df <- data.frame(genes=rownames(Cluster_markers[[i]]),
  #                 score=Cluster_markers[[i]]$ranking_score)
  df <- data.frame(Cluster_markers_final[[i]])
  print(head(df, 25))
  file_name <- names(Cluster_markers_final)[i]
  write.csv(df, paste0(marker_dir, file_name,'.csv'), 
            col.names = T, row.names = T, quote = F)
}
