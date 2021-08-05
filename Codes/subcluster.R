source('Codes/Functions.R')
Initialize()

#merged_samples <- readRDS('~/RatLiver/Objects/merged_samples_newSamples.rds')
merged_samples <- readRDS('Objects/merged_samples_newSamples_MT-removed.rds')

cell_types <- 'Immune cells'#'endothelial cells' #
### visualize the selected population
endothelial_clusters <- c(4, 6, 16)
Immune_clusters <- c(7, 8, 10, 11, 12, 13, 17)
cluster_name = as.character(Immune_clusters) #Immune_clusters
merged_samples$final_cluster = as.character(merged_samples$res.0.6)
umap_df <- data.frame(Embeddings(merged_samples, 'umap'))
umap_df$clusters = merged_samples$final_cluster
colnames(umap_df)[1:2] = c('UMAP_1', 'UMAP_2')
umap_df$cluster_interest <- ifelse(umap_df$clusters %in% cluster_name,umap_df$clusters , 'other')
ggplot(umap_df,aes(x= UMAP_1, y= UMAP_2, color=clusters))+geom_point()+theme_classic()+ggtitle('all clusters')
ggplot(umap_df,aes(x= UMAP_1, y= UMAP_2, color=cluster_interest))+geom_point()+theme_classic()+ggtitle(cell_types)

gene_name <- 'Ptprc'
umap_df$gene <- GetAssayData(merged_samples)[gene_name,]
ggplot(umap_df,aes(x= UMAP_1, y= UMAP_2, color=gene))+geom_point()+theme_classic()+
  ggtitle('Immune cells')+scale_color_viridis(direction = -1)

cluster3_sub_df$umi <- rownames(cluster3_sub_df)
head(cluster3_sub_df)
umap_df$subclusters <- ifelse(rownames(umap_df) %in% cluster3_sub_df$umi ,cluster3_sub_df$subcluster , 'other')
ggplot(umap_df,aes(x= UMAP_1, y= UMAP_2, color=subclusters))+geom_point()+
  theme_classic()+ggtitle(cluster_name)





### constructing a subset the input seurat object
UMIs <- colnames(merged_samples[,merged_samples$final_cluster %in% cluster_name])
seur = merged_samples[,colnames(merged_samples) %in% UMIs]

seur@assays$SCT <- NULL
seur <- SetAssayData(
  object = seur,
  assay.type = 'RNA',
  new.data = GetAssayData(seur)[,colnames(seur) %in% UMIs],
  slot = 'data'
)

colnames(seur@meta.data)
seur@meta.data <- seur@meta.data[,colnames(seur@meta.data) %in% c("orig.ident","nCount_RNA","nFeature_RNA","mito_perc",
                                                                  "nCount_SCT","nFeature_SCT","cell_type","sample_type",
                                                                  "strain_type", 'sample_name', 'strain', 'final_cluster')]
seur@reductions$pca <- NULL
seur@reductions$tsne <- NULL
seur@reductions$umap <- NULL
seur@reductions$harmony <- NULL

seur <- SCTransform(seur,conserve.memory=F,verbose=T,
                    return.only.var.genes=F,variable.features.n = nrow(seur))

DefaultAssay(seur) = 'SCT'
seur <- RunPCA(seur,verbose=T, features=rownames(seur))

## PCA

plot(100 * seur@reductions$pca@stdev^2 / seur@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")
PC_NUMBER = 20

### Harmony
#seur <- RunHarmony(seur, "sample_type",assay.use="RNA")
seur <- RunHarmony(seur, "sample_name",assay.use="SCT")
seur <- RunUMAP(seur,dims=1:PC_NUMBER, reduction="harmony")
seur <- RunTSNE(seur,dims=1:PC_NUMBER, reduction="harmony")


max_seurat_resolution <- 2.5 ## change this to higher values
FDRthresh <- 0.01 # FDR threshold for statistical tests
min_num_DE <- 1
seurat_resolution <- 0 # Starting resolution is this plus the jump value below.
seurat_resolution_jump <- 0.2

DefaultAssay(seur) = 'RNA' # why? shouldn't the data be normalized first? 
seur <- FindNeighbors(seur,reduction="harmony",dims=1:PC_NUMBER,verbose=F)


sCVdata_list <- list()
DE_bw_clust <- TRUE
while(DE_bw_clust) {
  if (seurat_resolution >= max_seurat_resolution) { break }
  seurat_resolution <- seurat_resolution + seurat_resolution_jump 
  # ^ iteratively incrementing resolution parameter 
  
  seur <- FindClusters(seur,resolution=seurat_resolution,verbose=F)
  
  message(" ")
  message("------------------------------------------------------")
  message(paste0("--------  res.",seurat_resolution," with ",
                 length(levels(Idents(seur)))," clusters --------"))
  message("------------------------------------------------------")
  
  if (length(levels(Idents(seur))) <= 1) { 
    message("Only one cluster found, skipping analysis.")
    next 
  } 
  # ^ Only one cluster was found, need to bump up the resolution!
  
  if (length(sCVdata_list) >= 1) {
    temp_cl <- length(levels(Clusters(sCVdata_list[[length(sCVdata_list)]])))
    if (temp_cl == length(levels(Idents(seur)))) { 
      temp_cli <- length(levels(interaction(
        Clusters(sCVdata_list[[length(sCVdata_list)]]),
        Idents(seur),
        drop=T
      )))
      if (temp_cli == length(levels(Idents(seur)))) { 
        message("Clusters unchanged from previous, skipping analysis.")
        next 
      }
    }
  }
  
  curr_sCVdata <- CalcSCV(
    inD=seur,
    assayType="RNA",
    assaySlot="counts",
    cl=Idents(seur), 
    # ^ your most recent clustering results get stored in the Seurat "ident" slot
    exponent=NA, 
    # ^ going to use the corrected counts from SCTransform
    pseudocount=NA,
    DRthresh=0.1,
    DRforClust="harmony",
    calcSil=T,
    calcDEvsRest=T,
    calcDEcombn=T
  )
  
  DE_bw_NN <- sapply(DEneighb(curr_sCVdata,FDRthresh),nrow)
  # ^ counts # of DE genes between neighbouring clusters at your selected FDR threshold
  message(paste("Number of DE genes between nearest neighbours:",min(DE_bw_NN)))
  
  if (min(DE_bw_NN) < min_num_DE) { DE_bw_clust <- FALSE }
  # ^ If no DE genes between nearest neighbours, don't loop again.
  
  sCVdata_list[[paste0("res.",seurat_resolution)]] <- curr_sCVdata
}

saveRDS(seur, 'Results/new_samples/Immune_subclusters_c17Included.rds')
new_data_scCLustViz_object_Immune <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_ImmuneSub_c17Included.RData"
save(seur,sCVdata_list,file=new_data_scCLustViz_object_Immune)

# seur@meta.data <- cbind(seur@meta.data, data.frame(lapply(sCVdata_list, function(x) x@Clusters)))

sCVdata_list$res.0.6@Clusters
umap_df$umi = rownames(umap_df)
subcluster_df <- data.frame(subclust=sCVdata_list$res.0.6@Clusters, 
                            umi=names(sCVdata_list$res.0.6@Clusters))


seur_2 <- readRDS('Results/new_samples/Immune_subclusters.rds') ### Immune_ , endothelial_
seur_umap <- data.frame(getEmb(seur,'umap'), cluster=as.character(seur$SCT_snn_res.1))
seur_umap <- data.frame(getEmb(seur,'umap'), cluster=as.character(seur$final_cluster))

ggplot(seur_umap, aes(UMAP_1, UMAP_2, color=cluster))+geom_point()+theme_classic()+scale_color_manual(values=colorPalatte)



head(subcluster_df)
head(umap_df)
umap_df <- merge(umap_df, subcluster_df, by.x='umi', by.y='umi', all.x=T)
head(umap_df)
ggplot(umap_df,aes(x= UMAP_1, y= UMAP_2, color=subclust))+geom_point()+
  theme_classic()+ggtitle(cluster_name)

markers_df_list = lapply(sCVdata_list$res.0.6@DEvsRest, function(markers_df){
  markers_df_2 <- markers_df[order(markers_df$logGER*(-log10(markers_df$FDR)),decreasing = T),]
  Rp_genes_index <- grep(pattern = 'Rp', x = rownames(markers_df_2))
  Mt_genes_index <- grep(pattern = 'Mt-', x = rownames(markers_df_2))
  markers_df_2 <- markers_df_2[-c(Rp_genes_index, Mt_genes_index),]
  return(markers_df_2)})

names(markers_df_list) = paste0('subcluster_', names(markers_df_list))
lapply(markers_df_list, head, 15)

merged_samples_2 <- merged_samples
subcluster.df.total = data.frame(umi=colnames(merged_samples_2) ,cluster=Idents(merged_samples_2))

subcluster.df.oneClust = data.frame(umi=names(sCVdata_list$res.0.56@Clusters) ,
                                    cluster=paste0('subcluster_', as.character(sCVdata_list$res.0.56@Clusters)))
subcluster.df.total <- merge(subcluster.df.total, subcluster.df.oneClust, by.x='umi', by.y='umi', all.x=T)
colnames(subcluster.df.total)[2:3] = c('cluster', 'subcluster')
head(subcluster.df.total)
is_other_cluster <- is.na(subcluster.df.total$subcluster)
subcluster.df.total$subcluster[is_other_cluster] <- as.character(subcluster.df.total$cluster[is_other_cluster])
head(subcluster.df.total)

Idents(merged_samples_2) = subcluster.df.total$subcluster

subcluster_index = 5
features <- unique(as.character(rownames(markers_df_list[[subcluster_index]])[1:60]))

DotPlot(merged_samples_2, features = features) + RotatedAxis() + 
  ggtitle(paste0(cluster_name, ' ', names(markers_df_list)[subcluster_index], ' markers'))


dir.create('~/RatLiver/Results/subclusters')
dir.create('~/RatLiver/Results/subclusters/new_samples/')
saveRDS(sCVdata_list, paste0('~/RatLiver/Results/subclusters/new_samples/',cluster_name,'_subclust.rds'))



dir.create('~/XSpecies/objects/merged_subclusters')
saveRDS(sCVdata_list, '~/XSpecies/objects/merged_subclusters/cluster_3_subclust.rds')

sCVdata_list <- readRDS('~/XSpecies/objects/merged_subclusters/cluster_3_subclust.rds')
lapply(sCVdata_list$res.0.02@Clusters, head)

cluster3_sub_df <- data.frame(subcluster=sCVdata_list$res.0.02@Clusters)
head(cluster3_sub_df)



sCVdata_list <- readRDS('~/XSpecies/objects/merged_subclusters/cluster_5_subclust.rds')
table(sCVdata_list$res.0.1@Clusters)
lapply(sCVdata_list$res.0.1@DEvsRest, head)

#### checking the top markers of sub-clusters 
seur_genes_df <- mapper
Cluster_markers <- sCVdata_list$res.0.17@DEvsRest
Cluster_markers_merged <- sapply(1:length(Cluster_markers), 
                                 function(i){
                                   markers_df <- Cluster_markers[[i]]
                                   markers_df$ensemble_ids = rownames(markers_df)
                                   ## merge the ensemble IDs in the dataframe with the HUGO terms 
                                   markers_df_merged <- merge(markers_df, seur_genes_df, 
                                                              by.x='ensemble_ids', 
                                                              by.y='V1', all.x=T, all.y=F,sort=F)
                                   #markers_df_merged2 <-  markers_df_merged[match(markers_df_merged$ensemble_ids, markers_df$ensemble_ids),]
                                   markers_df_merged2 <- markers_df_merged[order(markers_df_merged$logGER*(-log10(markers_df_merged$FDR)),decreasing = T),]
                                   return(markers_df_merged2)
                                 }, simplify = FALSE)

names(Cluster_markers_merged) = names(sCVdata_list$res.0.1@DEvsRest)         

lapply(Cluster_markers_merged, function(x)head(x, 20))               
