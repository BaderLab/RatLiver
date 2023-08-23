

## importing the gene expression data
merged_samples = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
Resolution = 2.5
merged_samples <- FindClusters(merged_samples, resolution = Resolution, verbose = FALSE)
table(merged_samples$SCT_snn_res.2.5)
merged_samples$clusters = merged_samples$SCT_snn_res.2.5


sample_types = names(table(merged_samples$sample_name))
median_exp_gene = sapply(1:length(sample_types), function(i) 
  median(merged_samples$nFeature_RNA[merged_samples$sample_name==sample_types[i]]), simplify = F)
names(sample_types) = sample_types

cluster_num = '4'
gene_name = 'Cps1'  #Sirpa
rownames(merged_samples)[grep(gene_name, rownames(merged_samples))] ## check if the gene is present in the dataset

df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap_h')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap_h')[,2], 
                      clusters=merged_samples$clusters,
                      expression_val=GetAssayData(merged_samples)[gene_name,],
                      a_cluster=merged_samples$cluster==cluster_num,
                      umi=colnames(merged_samples))


annot_info <- read.csv('figure_panel/Annotations_SingleNuc_Rat_July26.csv')
head(annot_info)
1:33 %in% annot_info$Cluster
annot_info <- annot_info[1:34, 1:4]
colnames(annot_info)[1] = 'clusters'

nrow(annot_info)
########## merging df_umap with annot_info data.frame
df_umap = merge(df_umap, annot_info, by.x='clusters', by.y='clusters', all.x=T, order=F)
#### re-ordering the rows based on UMIs of the initial data
df_umap <- df_umap[match(colnames(merged_samples),df_umap$umi),]
df_umap$umi == colnames(merged_samples)
sum(df_umap$umi != colnames(merged_samples))
df_umap$label_clust = paste0(df_umap$label, ' (',df_umap$clusters, ')')

colnames(merged_samples) == df_umap$umi


merged_samples$annotation = df_umap$label
merged_samples$cluster = df_umap$label_clust
merged_samples$cluster_number = as.character(merged_samples$SCT_snn_res.2.5)
table(as.character(merged_samples$SCT_snn_res.2.5), merged_samples$cluster )

merged_samples$strain
merged_samples$sample_name

head( merged_samples@meta.data[,c( "nCount_RNA","nFeature_RNA","mito_perc","sample_name",  
                                   "strain", "cluster",'cluster_number', 'annotation')])
merged_samples@meta.data = merged_samples@meta.data[,c( "nCount_RNA","nFeature_RNA","mito_perc","sample_name",  
                                                       "strain", "cluster",'cluster_number', 'annotation')]


merged_samples$sample_name=ifelse(merged_samples$sample_name=='DA_SHAM003_CST', 'DA-1', merged_samples$sample_name)
merged_samples$sample_name=ifelse(merged_samples$sample_name=='DA_SHAM004_CST', 'DA-2', merged_samples$sample_name)
merged_samples$sample_name=ifelse(merged_samples$sample_name=='LEW_SHAM1_CST', 'LEW-1', merged_samples$sample_name)
merged_samples$sample_name=ifelse(merged_samples$sample_name=='LEW_SHAM2_CST', 'LEW-2', merged_samples$sample_name)
table(merged_samples$sample_name)

saveRDS(merged_samples, 'cell_browser/snRNAseq_ratLiver_cellBrowser.rds')


########################
obj = readRDS('~/Spatial_Rat/Spatial_MacParland_Catia__Rat_6-1_A1_seurObj.rds')
obj@meta.data$PC1_score <- obj@reductions$pca@cell.embeddings[,1]
obj@meta.data$PC2_score <- obj@reductions$pca@cell.embeddings[,2]

obj@meta.data <- obj@meta.data[,c("nCount_Spatial","nFeature_Spatial", "percent.mt", "PC1_score", "PC2_score" )]
colnames(obj@meta.data)[c(4,5)] =  c("Zonation score (PC1)","Zonation score (PC2)")
head(obj@meta.data)
SpatialFeaturePlot(obj, features = c("Zonation score (PC1)","Zonation score (PC2)"))
saveRDS(obj, 'cell_browser/Spatial_Rat_6-1_A1_cellBrowser.rds')


########################
obj = readRDS('~/Spatial_Rat/Spatial_MacParland_Catia__Rat_6-1_B1_seurObj.rds')
obj@meta.data$PC1_score <- obj@reductions$pca@cell.embeddings[,1]
obj@meta.data$PC2_score <- obj@reductions$pca@cell.embeddings[,2]

obj@meta.data <- obj@meta.data[,c("nCount_Spatial","nFeature_Spatial", "percent.mt", "PC1_score", "PC2_score" )]
colnames(obj@meta.data)[c(4,5)] =  c("Zonation score (PC1)","Zonation score (PC2)")
head(obj@meta.data)
SpatialFeaturePlot(obj, features = c("Zonation score (PC1)","Zonation score (PC2)"))
saveRDS(obj, 'cell_browser/Spatial_Rat_6-1_B1_cellBrowser.rds')


obj=readRDS('cell_browser/Spatial_Rat_6-1_B1_cellBrowser.rds')
GetTissueCoordinates(obj)
obj$`Zonation score (PC1)`
########################

merged_samples <- readRDS('cell_browser/snRNAseq_ratLiver_cellBrowser.rds')
head(merged_samples@meta.data )

library(Seurat)
obj = readRDS('~/RatLiver/cell_browser/Spatial_Rat_6-1_A1_cellBrowser.rds')
SpatialFeaturePlot(obj, features = c("Zonation score (PC1)"), image.alpha=0)+scale_fill_gradient(low="white", high="tomato3")
SpatialFeaturePlot(obj, features = c("Alb","Cd68"))
head(obj@meta.data )

obj2 = readRDS('cell_browser/Spatial_Rat_6-1_B1_cellBrowser.rds')
head(obj2@meta.data )
SpatialFeaturePlot(obj2, features = c("Zonation score (PC1)"))



######################## MAKING FIGURE FOR GRAPHICAL ABSTRACT
merged_samples = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
merged_samples = readRDS('cell_browser/TLH_cellBrowser_2.rds')
merged_samples = readRDS('cell_browser/TLH_cellBrowser/TLH_cellBrowser.rds')
merged_samples = readRDS('cell_browser/Immune_enriched_cellBrowser_2.rds')
df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap_h')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap_h')[,2])

df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap')[,2])
a_color = 'moccasin'#gray50' #lavender' #'seashell3'
ggplot(df_umap, aes(UMAP_1, UMAP_2))+geom_point(color=a_color, size=0.8)+theme_classic()




rot_data <- readRDS('~/RatLiver/Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds') ## MT-removed

scores <- data.frame(rot_data$rotScores)
colnames(scores) = paste0('Varimax_', 1:ncol(scores))
embedd_df_rotated <- data.frame(scores)

pc_num = 5
rot_df <- data.frame(Varimax_1=embedd_df_rotated$Varimax_1,
                     emb_val=embedd_df_rotated[,pc_num],
                     strain=merged_samples$strain
)

ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=strain))+geom_point(alpha=0.7,size=4)+
  theme_classic()+ylab(paste0('Varimax ',pc_num))+xlab("Varimax 1")+
  scale_color_manual(name='Immune\nsubcluster',values = mycolors)+theme_classic()+
  theme(text = element_text(size=22))+#, legend.title = element_blank()
  ggtitle(paste0('Set-2 Immune-subclusters over Varimax 1 and ', pc_num))
