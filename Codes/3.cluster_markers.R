source('Codes/Functions.R')
source('Codes/FactorAnalysisUtils.R')
Initialize()
library(plyr)


map = 'new' #'old' 
#map='subcluster'

output_file = 'mt_removed_immuneSub_c17Inc'
if(map == 'old') output_file = 'mt40_lib1500_MTremoved'

merged_samples = NULL
if(map=='old') merged_samples = readRDS('Objects/merged_samples_oldSamples_mt40_lib1500_MTremoved.rds')
if(map=='new') merged_samples = readRDS('Objects/merged_samples_newSamples_MT-removed.rds')
table(merged_samples$sample_name)

merged_samples <- readRDS('Results/new_samples/Immune_subclusters.rds') 
merged_samples <- readRDS('Results/new_samples/endothelial_subclusters.rds') 

#new_data_scCLustViz_object_Immune <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_ImmuneSub.RData"
new_data_scCLustViz_object_Immune <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_ImmuneSub_c17Included.RData"

new_data_scCLustViz_object_endothelial <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_EndothelialSub.RData"

load(new_data_scCLustViz_object_Immune)
merged_samples = your_scRNAseq_data_object
merged_samples = seur
#### finding markers for the clusters ####

#### final clusters - old samples: resolution 0.6, 17 clusters >> running on screen: runner4
### Heps: any clusters other than 13, 9, 5, 7, 14, 3 


#### final clusters - new samples: resolution 0.6, 18 clusters
### Heps: 0, 1, 2, 3, 5, 9, 14, 15

if(map=='new') merged_samples$cluster = as.character(merged_samples$res.0.6)
if(map=='old') merged_samples$cluster = as.character(merged_samples$res.0.6)
if(map=='subcluster') merged_samples$cluster = as.character(sCVdata_list$RNA_snn_res.1@Clusters)
if(map=='subcluster') merged_samples$cluster = as.character(sCVdata_list$res.1@Clusters)

Idents(merged_samples) <- paste0('cluster_', as.character(merged_samples$cluster))
cluster_names <-  levels(merged_samples)
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
#lapply(Cluster_markers_final, head)
#lapply(Cluster_markers_final, dim)

#Cluster_markers <- readRDS('Results/old_samples/Cluster_markers_mergedOldSamples.rds')
#saveRDS(Cluster_markers_final, 'Results/new_samples/Cluster_markers_final.rds')
saveRDS(Cluster_markers_final, paste0('Results/',map,'_samples/Cluster_markers_', output_file, '.rds'))


#### saving markers #####
dir.create(paste0('Results/',map,'_samples/final_markers_',output_file, '/'))
for(i in 1:length(Cluster_markers_final)){
  #df <- data.frame(genes=rownames(Cluster_markers[[i]]),
  #                 score=Cluster_markers[[i]]$ranking_score)
  df <- data.frame(Cluster_markers_final[[i]])
  print(head(df, 25))
  file_name <- names(Cluster_markers_final)[i]
  write.csv(df, paste0('Results/',map,'_samples/final_markers_',output_file, '/', file_name,'.txt'), 
            col.names = T, row.names = T, quote = F)
}



### check the data using UMAP - sanity check
df <- data.frame(Embeddings(merged_samples, 'umap'), 
                 sample=merged_samples$sample_name, 
                 cluster=as.character(merged_samples$res.1) )
df$isCluster = df$cluster == '22'                       
ggplot(df, aes(UMAP_1, UMAP_2, color=isCluster))+geom_point()+theme_classic()
