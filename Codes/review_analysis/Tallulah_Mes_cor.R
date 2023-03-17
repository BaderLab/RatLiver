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
################ Rat - sinle nuc seq data ################
########################################################################
## importing the gene expression data
merged_samples = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
Resolution = 0.6
resolutions = Resolution
merged_samples <- FindClusters(merged_samples, resolution = Resolution, verbose = FALSE)
table(merged_samples$SCT_snn_res.0.6)



merged_samples <- readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
table(merged_samples$cluster, merged_samples$annot_IM) ## cluster 8 seems to be the macrophage population
mes_cluster_num = c(5,7)

##################################################################
######## Subclustering the mesenchymal/endothelilal cluster

mes_data = merged_samples[,merged_samples$cluster %in% mes_cluster_num]
mes_data_meta = mes_data@meta.data
DefaultAssay(mes_data) <- 'RNA'
mes_data@assays$SCT <- NULL

mes_data = SCTransform(mes_data, vst.flavor = "v2", verbose = TRUE)
mes_data <- RunPCA(mes_data,verbose=T)
plot(100 * mes_data@reductions$pca@stdev^2 / mes_data@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")
mes_data <- RunHarmony(mes_data, group.by.vars = "sample_name", assay.use="RNA")

res = 1.0
mes_data <- RunUMAP(mes_data, reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = res, verbose = FALSE)

dim(mes_data)
markers = c('Lyve1')
i = 1
gene_name = markers[i] 
df_umap <- data.frame(UMAP_1=getEmb(mes_data, 'umap')[,1], 
                      UMAP_2=getEmb(mes_data, 'umap')[,2], 
                      library_size= mes_data$nCount_RNA, 
                      mito_perc=mes_data$mito_perc, 
                      n_expressed=mes_data$nFeature_RNA,
                      cluster=mes_data$cluster, 
                      cell_status = mes_data$cell_status,
                      nuclear_fraction=mes_data$nuclear_fraction, 
                      Alb=GetAssayData(mes_data)['Alb',], 
                      a_gene = GetAssayData(mes_data)[gene_name,],
                      sample_name = mes_data$sample_name, 
                      strain = mes_data$strain, 
                      umi=colnames(mes_data))

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=a_gene))+geom_point(size=1.5, alpha=0.6)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(gene_name)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cluster))+geom_point(size=1.5, alpha=0.6)+theme_classic()

saveRDS(mes_data, '~/rat_sham_sn_data/standardQC_results/sham_sn_mesenchymal_subclusters_merged_data.rds')
###################################################################



rat_HVGs <- VariableFeatures(mes_data)
head(mes_data)
table(mes_data$annot_IM_g,mes_data$SCT_snn_res.1)

##### annotate the hepatocytes based on the mouse zonation layers
mes_data$cluster <- as.character(mes_data$SCT_snn_res.1)
cluster_names = mes_data$cluster 
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
mes_cluster_average_exp <- get_scaled_by_gene(cluster_average_exp_df) ## scaling over the clusters of interest
mes_cluster_average_exp$rat_ID=rownames(mes_cluster_average_exp)
head(mes_cluster_average_exp)
dim(mes_cluster_average_exp)
mes_cluster_average_exp <- mes_cluster_average_exp[mes_cluster_average_exp$rat_ID %in% rat_HVGs,]

rat_cluster_average.df = mes_cluster_average_exp
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



#################
mes_data <- readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_mesenchymal_subclusters_merged_data.rds')

resolutions = seq(0.6, 1.4, 0.2)
for (res in resolutions){
  mes_data <- FindClusters(mes_data, resolution = res, verbose = FALSE)
}


head(mes_data@meta.data)
your_cluster_results =data.frame(mes_data@meta.data[,colnames(mes_data@meta.data) %in% paste0('SCT_snn_res.', resolutions)])
head(your_cluster_results)


sCVdata_list <- CalcAllSCV(
  inD=mes_data,
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

sham_sn_merged_scCLustViz_object =  '~/rat_sham_sn_data/standardQC_results/sham_sn_merged_mesenchymal_subclusters_sCVdata.RData' ### find the results on run1
#saveRDS(sCVdata_list, '~/rat_sham_sn_data/standardQC_results/sham_sn_merged_sCVdata.rds') ### find the results on run1
save(mes_data, sCVdata_list, 
     file=sham_sn_merged_scCLustViz_object) ## new data scClustViz object





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





