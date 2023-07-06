source('Codes/Functions.R')
source('~/XSpecies/Codes/convert_human_to_ortholog_functions.R')
library(Hmisc)
library(RColorBrewer)
library(viridis)
library(reshape2)

get_scaled_by_gene <- function(x){
  y <- scale(t(x), center = TRUE, scale = TRUE)
  y[is.nan(y)] <- 0
  y2 <- t(y)
  y2 <- as.data.frame(y2)
  return(y2)
}

########################################################################
################ Rat - sinle nuc seq data ################
########################################################################
## importing the gene expression data
merged_samples = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')

########################################################################
############### Analysis based on resolution 2.5 ##############

Resolution = 2.5
merged_samples <- FindClusters(merged_samples, resolution = Resolution, verbose = FALSE)
table(merged_samples$SCT_snn_res.2.5)
rat_HVGs <- VariableFeatures(merged_samples)


#### defining the mega clusters based on 0.6 resolution
Hep0 = as.character(c(0, 4, 7, 10, 16, 25, 27)) # 27 is the little tail
Hep1 = as.character(c(1, 6, 15, 17, 12, 31)) # 31 is the tail
Hep2 = as.character(c(2, 5, 9, 21)) #
Hep3 = as.character(c(3, 8,13, 23, 32 ))


merged_samples$clusters = as.character(merged_samples$SCT_snn_res.2.5)
merged_samples$Hep_clusters = ifelse(merged_samples$clusters %in% as.character(Hep0), 'Hep0', merged_samples$clusters) 
merged_samples$Hep_clusters = ifelse(merged_samples$Hep_clusters %in% as.character(Hep1), 'Hep1', merged_samples$Hep_clusters) 
merged_samples$Hep_clusters = ifelse(merged_samples$Hep_clusters %in% as.character(Hep2), 'Hep2', merged_samples$Hep_clusters) 
merged_samples$Hep_clusters = ifelse(merged_samples$Hep_clusters %in% as.character(Hep3), 'Hep3', merged_samples$Hep_clusters) 
table(merged_samples$Hep_clusters)


meta_to_include =c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'mito_perc', 
                   'cell_status', 'nuclear_fraction', 'celltype_dQC', 'sample_name', 
                   'strain', 'Hep_clusters', 'SCT_snn_res.2.5' )

merged_samples_meta = merged_samples@meta.data
merged_samples@meta.data = merged_samples_meta[,colnames(merged_samples_meta) %in% meta_to_include]
head(merged_samples@meta.data)

hep_data = merged_samples[,merged_samples$SCT_snn_res.2.5 %in% c(Hep0, Hep1, Hep2, Hep3)]
hep_data_meta = hep_data@meta.data
head(hep_data_meta)
DefaultAssay(hep_data) <- 'RNA'
hep_data@assays$SCT <- NULL
GetAssayData(hep_data, 'counts')[100:200,]
GetAssayData(hep_data, 'data')[100:200,]

hep_data@reductions$pca <- NULL
hep_data@reductions$umap <- NULL
hep_data@reductions$umap_h <- NULL
hep_data@reductions$harmony <- NULL

hep_data = SCTransform(hep_data, vst.flavor = "v2", verbose = TRUE)
#hep_data = ScaleData(hep_data)

hep_data <- RunPCA(hep_data,verbose=T)
plot(100 * hep_data@reductions$pca@stdev^2 / hep_data@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")
hep_data <- RunHarmony(hep_data, group.by.vars = "sample_name")

res = 2.5
hep_data <- RunUMAP(hep_data, reduction = "harmony", dims = 1:30, verbose = FALSE) 
hep_data <-   FindNeighbors(hep_data, reduction = "harmony", dims = 1:30, verbose = FALSE)
hep_data <-  FindClusters(hep_data, resolution = res, verbose = FALSE)

dim(hep_data)
markers = c('Ptprc','Cd68', 'Cd163', 'Mrc1', 'Clec4f', 'Siglec5')

i = 5
gene_name = markers[i] 
gene_name = 'Clec4f'
df_umap <- data.frame(UMAP_1=getEmb(hep_data, 'umap')[,1], 
                      UMAP_2=getEmb(hep_data, 'umap')[,2], 
                      library_size= hep_data$nCount_RNA, 
                      mito_perc=hep_data$mito_perc, 
                      n_expressed=hep_data$nFeature_RNA,
                      cluster=hep_data$SCT_snn_res.2.5,
                      cell_status = hep_data$cell_status,
                      nuclear_fraction=hep_data$nuclear_fraction, 
                      Alb=GetAssayData(hep_data)['Alb',], 
                      a_gene = GetAssayData(hep_data)[gene_name,],
                      sample_name = hep_data$sample_name, 
                      Hep_clusters=hep_data$Hep_clusters,
                      strain = hep_data$strain, 
                      umi=colnames(hep_data))


ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cluster))+geom_point(size=1.7,alpha=0.8)+
  theme_classic()+#scale_color_manual('Cluster', values = c(colorPalatte))+
  theme(text = element_text(size=16.5),
        plot.title = element_text(hjust = 0.5),
        legend.title=element_text(size = 10))+
  ggtitle(paste0(''))#colorPalatte



ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=sample_name))+geom_point(size=1.3,alpha=0.7)+theme_classic()
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=Hep_clusters))+geom_point(size=1.3,alpha=0.7)+theme_classic()


ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=a_gene))+
  scale_color_viridis(direction = -1)+ggtitle(gene_name)+
  geom_point(alpha=0.8,size=1.6)+
  scale_color_viridis('Expresion\nValue', direction = -1)+theme_classic()+
  theme(text = element_text(size=16.5),
        plot.title = element_text(hjust = 0.5),
        legend.title=element_text(size = 10))+
  ggtitle(paste0(gene_name)) 


saveRDS(hep_data, '~/rat_sham_sn_data/standardQC_results/sham_sn_merged_Hep_subclusters_june30.rds')
########################################################################################
########### Make the scClustViz object of the Hep subcluster ##################
########################################################################################

resolutions = seq(0.6, 1.4, 0.2)
for (res in resolutions){
  hep_data <- FindClusters(hep_data, resolution = res, verbose = FALSE)
}


head(hep_data@meta.data)
your_cluster_results =data.frame(hep_data@meta.data[,colnames(hep_data@meta.data) %in% paste0('SCT_snn_res.', resolutions)])
head(your_cluster_results)


sCVdata_list <- CalcAllSCV(
  inD=hep_data,
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

#sham_sn_merged_scCLustViz_object =  '~/rat_sham_sn_data/standardQC_results/sham_sn_merged_macrophage_subclusters_sCVdata.RData' ### find the results on run1
save(hep_data, sCVdata_list, 
     file=sham_sn_merged_scCLustViz_object) ## new data scClustViz object

##################################################################



#############################################################################################
############### Making marker heatmap for the Hep subcluster data. ###########################
#############################################################################################
PP1_list= c("Sds", "Hal", "Cyp2f4", "Arg1", "Acly", "Ass1", "Gls2", "Agxt", "Uroc1", "Gldc", "Gls2")
PP2_list= c("Apoc3", "Serpina1", "Apoc1", "Apoe", "Itih4", "Apoa1", "Ttr", "Tf", "Alb", "Pigr", "Orm1", "Rpl3", "Fads1", "Aldh1b1", "Srd5a1", "Hsd17b13")
IZ_list= c("Saa4", "Hint1", "Cyp8b1", "Cyp2e1", "Cox7c", "Fabp1")
CV_list= c("Notum", "Cyp27a1", "Fabp7", "Akr1c1", "Gsta5", "Slc22a1", "Aox3",
           "Sult1e1", "Fmo1", "Oat", "Ahr", "Cyp7a1", "Glul", "Rhbg", "Cyp2e1", "Cyp1a2", "Por")
markers_list = c(CV_list, IZ_list, PP2_list, PP1_list)

################################
##### Generating average expression based heatmap

get_scaled_by_gene <- function(x){
  y <- scale(t(x), center = TRUE, scale = TRUE)
  y[is.nan(y)] <- 0
  y2 <- t(y)
  y2 <- as.data.frame(y2)
  return(y2)
}



merged_samples_subset = hep_data
#merged_samples_subset2 <- CreateAssayObject(GetAssayData(merged_samples_subset)[rownames(merged_samples_subset) %in% markers_list,])
merged_samples_subset <- merged_samples_subset[rownames(merged_samples_subset) %in% markers_list,]

dim(merged_samples_subset)
#merged_samples_subset2$manual_annotation = merged_samples_subset$manual_annotation
cluster_names = as.character(merged_samples_subset$SCT_snn_res.2.5 )
cluster_names_types = names(table(cluster_names))

### dividing the expression matrix based on the clusters
cluster_expression <- sapply(1:length(cluster_names_types), function(i){
  a_cluster_name = cluster_names_types[i]
  GetAssayData(merged_samples_subset, 'data')[,cluster_names == a_cluster_name] 
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
colnames(cluster_average_exp_df) = paste0(names(cluster_average_exp))
head(cluster_average_exp_df)


## scale and center all the genes in the matrix
cluster_average_exp_df_scaled <- get_scaled_by_gene(cluster_average_exp_df) ## scaling over the clusters of interest
colnames(cluster_average_exp_df_scaled)
rownames(cluster_average_exp_df_scaled)


anno_df_genes = c(rep(NA, nrow(cluster_average_exp_df_scaled)))
sapply(1:ncol(merged_samples_subset) , function(i){
  if (rownames(cluster_average_exp_df_scaled)[i] %in% CV_list) anno_df_genes[i] <<- 'CV'
  if (rownames(cluster_average_exp_df_scaled)[i] %in% IZ_list) anno_df_genes[i] <<- 'IZ'
  if (rownames(cluster_average_exp_df_scaled)[i] %in% PP2_list) anno_df_genes[i] <<- 'PP2'
  if (rownames(cluster_average_exp_df_scaled)[i] %in% PP1_list) anno_df_genes[i] <<- 'PP1'
})
length(anno_df_genes)

##### ordering clusters
#cluster_orders = as.character(c(Hep3, Hep2, Hep0, Hep1))
cluster_orders = as.character(colnames(cluster_average_exp_df_scaled))
#anno_df = data.frame(cluster=colnames(cluster_average_exp_df_scaled),annot=anno_df)

cluster_average_exp_df_scaled = cluster_average_exp_df_scaled[,cluster_orders]


#### ordering genes
gene_orders = c(CV_list, IZ_list, PP2_list, PP1_list) 

gene_orders = unique(gene_orders[gene_orders %in% rownames(cluster_average_exp_df_scaled)])
anno_df_genes = data.frame(genes=row.names(cluster_average_exp_df_scaled),annot=anno_df_genes)


cluster_average_exp_df_scaled = cluster_average_exp_df_scaled[gene_orders,]
dim(cluster_average_exp_df_scaled)
anno_df_genes = anno_df_genes[match(row.names(cluster_average_exp_df_scaled), anno_df_genes$genes),]



head(cluster_average_exp_df_scaled)
rownames(cluster_average_exp_df_scaled)
pheatmap(t(cluster_average_exp_df_scaled), cluster_rows = FALSE, color= inferno(20, direction = +1),  
         cluster_cols = FALSE)


