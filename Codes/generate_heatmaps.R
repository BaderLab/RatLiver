### generating heatmap based on each cluster markers and the seurat default heatmap function

## importing the markers
Cluster_markers_merged <- readRDS('~/XSpecies/objects/Cluster_markers_merged.rds')

## importing the gene expression data
merged_samples <- readRDS('Results/preproc_rats/merged/merged_rat_samples_2.rds')
DefaultAssay(your_scRNAseq_data_object) <- 'SCT'


num_genes <- 35

pdf('plots/top_marker_heatmaps_mergedRats.pdf', height = 7, width = 13)
for(i in 1:length(Cluster_markers_merged)){
  
  cluster_name = names(Cluster_markers_merged)[i]
  
  p=DoHeatmap(object = your_scRNAseq_data_object, 
            features = Cluster_markers_merged[[i]]$V2[1:num_genes], 
            label = TRUE, 
            angle = 90, size=3)+
    ggtitle(paste0(cluster_name, ' Heatmap'))+
    theme(plot.title = element_text(color = "blue", size = 12, face = "bold", hjust=0.23))
  print(p)
}

dev.off()




###### part of jeff's code related to the heatmap generation

#heatmap of top 15 genes per cluster compared to all others
deG<-get(load(paste0(timePoint,"_precalc_",gsub(".","",names(minDEgenes)[f3],fixed=T),"_deGall.RData")))
heatGenes <- lapply(deG,function(X) rownames(X[order(X$fwer),])[1:15])
clustMeans <- sapply(CGS,function(X) X$MTC[X$genes %in% unique(unlist(heatGenes))])
rownames(clustMeans) <- CGS[[1]]$genes[CGS[[1]]$genes %in% unique(unlist(heatGenes))]
hG <- hclust(dist(clustMeans),"complete")
hC <- hclust(dist(t(clustMeans)),"single")
tempLabCol <- paste(paste0("#",seq_along(deG)),
                    paste(sapply(deG,nrow),"DE"),sep=": ")
clustMeans2 <- sapply(CGS,function(X) X$MTC[X$genes %in% unique(markers[,2])])
rownames(clustMeans2) <- CGS[[1]]$genes[CGS[[1]]$genes %in% unique(markers[,2])]
colnames(clustMeans2) <- paste0("Clust #",seq_along(deG))
ord1<-hclust(dist(t(clustMeans2)))
clustMeans3<-clustMeans2[,ord1$order]
Labels3<-cbind(rep("Col Labels",ncol(clustMeans3)),colnames(clustMeans3))
heatmap.2(clustMeans,Rowv=as.dendrogram(hG),Colv=as.dendrogram(hC),scale="row",
          col="viridis",trace="none",
          ColSideColors=clustCols,labCol=tempLabCol)
readline(prompt="Press Enter to continue: ")

