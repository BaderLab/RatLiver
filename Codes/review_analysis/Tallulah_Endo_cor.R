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
endo_data = readRDS('Objects/Endothelial/LSEC_harmony_cleaned.RDS')
endo_seur = GetAssayData(endo_data)

endo_data <- FindVariableFeatures(endo_data)
human_HVGs <- VariableFeatures(endo_data)
endo_data = endo_data[rownames(endo_data) %in% human_HVGs,] 
cluster_names_types = names(table(endo_data$sub_annotations))
cluster_names = endo_data$sub_annotations

### dividing the expression matrix based on the clusters
cluster_expression <- sapply(1:length(cluster_names_types), function(i){
  a_cluster_name = cluster_names_types[i]
  GetAssayData(endo_data, 'data')[,cluster_names == a_cluster_name] 
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
dim(human_data)

########################################################################
################ Rat - sinle nuc seq data ################
########################################################################
## importing the gene expression data
merged_samples = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
Resolution =  2.5 #0.6
resolutions = Resolution
merged_samples <- FindClusters(merged_samples, resolution = Resolution, verbose = FALSE)
table(merged_samples$SCT_snn_res.0.6)
table(merged_samples$SCT_snn_res.2.5)
merged_samples$cluster = as.character(merged_samples$SCT_snn_res.2.5)

table(merged_samples$cluster, merged_samples$annot_IM) ## cluster 8 seems to be the macrophage population
#mes_cluster_num = c(5, 7) # 0.6 res
endo_cluster_num = c(11,30) # 2.5 res
endo_cluster_num = c(11)
##################################################################
######## Subclustering the endoenchymal/endothelilal cluster

endo_data = merged_samples[,merged_samples$cluster %in% endo_cluster_num]
dim(endo_data)
endo_data_meta = endo_data@meta.data

##### removing the cells which will be subcluster-2 after subclustering (based on Sai's request)
endo_data = endo_data[,!colnames(endo_data) %in% endo_subcluster_c2_UMIs]


DefaultAssay(endo_data) <- 'RNA'
endo_data@assays$SCT <- NULL

endo_data = SCTransform(endo_data, vst.flavor = "v2", verbose = TRUE)

endo_data[["pca"]] <- NULL
endo_data[["harmony"]] <- NULL
endo_data[["umap"]] <- NULL
endo_data[["umap_h"]] <- NULL
endo_data@meta.data <- endo_data@meta.data[,!grepl('SCT_snn', colnames(endo_data@meta.data))]

endo_data <- RunPCA(endo_data,verbose=T)
plot(100 * endo_data@reductions$pca@stdev^2 / endo_data@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")
table(endo_data$sample_name)


endo_data <- RunHarmony(endo_data, group.by.vars = "sample_name", assay.use="RNA")

res = 1.0
endo_data <- RunUMAP(endo_data, reduction = "harmony", dims = 1:30, verbose = FALSE, reduction.key = "UMAP_") %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = res, verbose = FALSE)

table(endo_data$cluster)
table(endo_data$cluster[endo_data$SCT_snn_res.1==2])
table(endo_data$SCT_snn_res.1)

markers = c('Lyve1')
i = 1
gene_name = markers[i] 
df_umap <- data.frame(UMAP_1=getEmb(endo_data, 'umap')[,1], 
                      UMAP_2=getEmb(endo_data, 'umap')[,2], 
                      library_size= endo_data$nCount_RNA, 
                      mito_perc=endo_data$mito_perc, 
                      n_expressed=endo_data$nFeature_RNA,
                      cluster=endo_data$SCT_snn_res.1, 
                      orig_cluster=endo_data$cluster,
                      cell_status = endo_data$cell_status,
                      nuclear_fraction=endo_data$nuclear_fraction, 
                      Alb=GetAssayData(endo_data)['Alb',], 
                      a_gene = GetAssayData(endo_data)[gene_name,],
                      sample_name = endo_data$sample_name, 
                      strain = endo_data$strain, 
                      umi=colnames(endo_data))

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=a_gene))+geom_point(size=1.5, alpha=0.6)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(gene_name)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cluster))+geom_point(size=1.5, alpha=0.6)+theme_classic()
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=orig_cluster))+geom_point(size=1.5, alpha=0.6)+theme_classic()

saveRDS(endo_data, '~/rat_sham_sn_data/standardQC_results/sham_sn_endothelial_subclusters_merged_data_res2.5.rds')
saveRDS(endo_data, '~/rat_sham_sn_data/standardQC_results/sham_sn_endothelial_subclusters_merged_data_res1_onlyC11.rds')

###################################################################

### figuring out the cluster of origin for subcluster 2 in the res=1 and removing it from the subclusterg - cluster 30??
### this mainly includes cells from cluster-30 of the original map
endo_subcluster_c2_UMIs = rownames(df_umap)[df_umap$cluster==2]
saveRDS(endo_subcluster_c2_UMIs, '~/rat_sham_sn_data/standardQC_results/endo_subcluster_cluster11_30_C2_UMIs.rds') 

##########

rat_HVGs <- VariableFeatures(endo_data)
head(endo_data)
table(endo_data$annot_IM_g,endo_data$SCT_snn_res.1)

##### annotate the hepatocytes based on the mouse zonation layers
endo_data$cluster <- as.character(endo_data$SCT_snn_res.1)
cluster_names = endo_data$cluster 
cluster_names_types = names(table(cluster_names))

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
endo_cluster_average_exp <- get_scaled_by_gene(cluster_average_exp_df) ## scaling over the clusters of interest
endo_cluster_average_exp$rat_ID=rownames(endo_cluster_average_exp)
head(endo_cluster_average_exp)
dim(endo_cluster_average_exp)
endo_cluster_average_exp <- endo_cluster_average_exp[endo_cluster_average_exp$rat_ID %in% rat_HVGs,]

rat_cluster_average.df = endo_cluster_average_exp
head(rat_cluster_average.df)
head(human_data)

############################################################
merge_human_rat = merge(human_data, rat_cluster_average.df, by.x='rat_genes', by.y='rat_ID')
head(merge_human_rat)
dim(merge_human_rat)

merge_human_rat = merge_human_rat[-c(1, 2)] #10
num_human_clusters = ncol(human_data) - 2
pheatmap::pheatmap(cor(merge_human_rat)[1:num_human_clusters, (num_human_clusters+1):ncol(merge_human_rat)])
dim(merge_human_rat)
############################################################
merge_mouse_rat = merge(mouse_data, rat_cluster_average.df, by.x='rat_genes', by.y='rat_ID')
head(merge_mouse_rat)
dim(merge_mouse_rat)

merge_mouse_rat = merge_mouse_rat[-c(1, 2)] #10
pheatmap::pheatmap(cor(merge_mouse_rat)[1:3, 4:ncol(merge_mouse_rat)])



####################################################################
endo_data = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_endothelial_subclusters_merged_data_res2.5.rds')
endo_data <- readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_endothelial_subclusters_merged_data_res1_onlyC11.rds')

endo_data@meta.data <- endo_data@meta.data[,!grepl('SCT_snn', colnames(endo_data@meta.data))]
table(endo_data$cluster)

head(endo_data)

resolutions = seq(0.4,2, 0.2)
for (res in resolutions){
  endo_data <- FindClusters(endo_data, resolution = res, verbose = FALSE)
}


head(endo_data@meta.data)
your_cluster_results =data.frame(endo_data@meta.data[,colnames(endo_data@meta.data) %in% paste0('SCT_snn_res.', resolutions)])
head(your_cluster_results)


sCVdata_list <- CalcAllSCV(
  inD=endo_data,
  clusterDF=your_cluster_results,
  assayType='SCT', #specify assay slot of data
  DRforClust="harmony",#reduced dimensions for silhouette calc
  #exponent=exp(1), #log base of normalized data
  #pseudocount=1,
  #DRthresh=0.5, #gene filter - minimum detection rate
  testAll=T, #stop testing clusterings when no DE between clusters
  #FDRthresh=0.005,
  #calcSil=F, #use cluster::silhouette to calc silhouette widths
  calcDEvsRest=T,
  calcDEcombn= T #
)

sham_sn_merged_scCLustViz_object = '~/rat_sham_sn_data/standardQC_results/sham_sn_merged_endothelial_subclusters_sCVdata_res2.5.RData' ### find the results on run1
sham_sn_merged_scCLustViz_object = '~/rat_sham_sn_data/standardQC_results/sham_sn_endothelial_subclusters_sCVdata_res1_onlyC11.RData'

#saveRDS(sCVdata_list, '~/rat_sham_sn_data/standardQC_results/sham_sn_merged_sCVdata.rds') ### find the results on run1
save(endo_data, sCVdata_list, 
     file=sham_sn_merged_scCLustViz_object) ## new data scClustViz object

### DE for clustering resolution of 0.8 and 1.4

load(sham_sn_merged_scCLustViz_object)

runShiny(
  ## write the path to the file bellow:
  filePath= sham_sn_merged_scCLustViz_object,
  outPath="./",
  # Save any further analysis performed in the app to the
  # working directory rather than library directory.
  annotationDB="org.Rn.eg.db",
  # This is an optional argument, but will add annotations.
  imageFileType="png"
  #Set the file format of any saved figures from the app.
)
