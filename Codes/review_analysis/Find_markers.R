source('Codes/Functions.R')
#source('Codes/FactorAnalysisUtils.R')
Initialize()
library(plyr)
library(Polychrome)




#merged_data <- readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_standardQC.rds')
merged_samples = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
#merged_samples <- FindNeighbors(merged_samples,reduction="harmony",verbose=T)
Resolution = 2.5
merged_samples <- FindClusters(merged_samples, resolution = Resolution, verbose = FALSE)
#merged_samples$cluster = as.character(merged_samples$SCT_snn_res.0.6)
merged_samples$cluster = as.character(merged_samples$SCT_snn_res.2.5) ## 34 clusters
table(merged_samples$cluster)

Idents(merged_samples) <- paste0('cluster_', as.character(merged_samples$cluster))
cluster_names <-  levels(merged_samples)


df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap_h')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap_h')[,2], 
                      #cluster=merged_samples$cluster, 
                      cluster=merged_samples$cluster, 
                      Alb=GetAssayData(merged_samples)['Alb',], 
                      sample_name = merged_samples$sample_name, 
                      strain = merged_samples$strain, 
                      umi=colnames(merged_samples))

# creating color palette with many distinct colors
P34 = createPalette(length(cluster_names),  c("#ff0000", "#00ff00", "#0000ff"))
swatch(P34)
color_pallete = as.vector(P34)



ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cluster))+geom_point(size=1,alpha=0.6)+theme_classic()+
  scale_color_manual(values=color_pallete)

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
#saveRDS(Cluster_markers_final, paste0('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_cluster_markers_res',Resolution,'.rds'))
Cluster_markers_final <- readRDS(paste0('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_cluster_markers_res',Resolution,'.rds'))
lapply(Cluster_markers_final, head)
#### saving markers #####
marker_dir = paste0('~/rat_sham_sn_data/standardQC_results/cluster_markers_res',Resolution,'/')
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




