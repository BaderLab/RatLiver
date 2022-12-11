# can then write a Seurat object to a directory from which you can run cbBuild:
library(Seurat)
?ExportToCellbrowser
packageVersion('Seurat')
source("https://raw.githubusercontent.com/maximilianh/cellBrowser/master/src/cbPyLib/cellbrowser/R/ExportToCellbrowser-seurat2.R")

ExportToCellbrowser(pbmc_small, dir="pbmcSmall", dataset.name="pbmcSmall")
#Or, you can build a from this dataset into the htdocs directory, serve the result on port 8080 via http, and open a web browser from within R:
ExportToCellbrowser(pbmc_small, dir="pbmcSmall", cb.dir="htdocs", dataset.name="pbmcSmall", port=8080)



### link didn't open - 

source('Codes/Functions.R')
Initialize()
####################################################################
######## prepare the TLH data for conversion ########
old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)

merged_samples <- readRDS('Objects/merged_samples_oldSamples_mt40_lib1500_MTremoved.rds')
colnames(merged_samples) == colnames(your_scRNAseq_data_object)  ### cells have been named correctly in this data 

Idents(merged_samples) = as.character(merged_samples$res.0.6)
head(merged_samples)

merged_samples$cluster = as.character(merged_samples$res.0.6)
merged_samples@meta.data[,grep(colnames(merged_samples@meta.data), pattern = 'res.')] = NULL
head(merged_samples@meta.data)
merged_samples@meta.data$orig.ident = NULL
head(merged_samples@meta.data)

### adding the cell-type annotation ###
annotation_TLH = read.csv('TLH_annotation.csv')
meta_data = merged_samples@meta.data
meta_data$cell_id = rownames(meta_data)
meta_data2 = merge(meta_data, annotation_TLH, by.x='cluster', by.y='cluster', all.x=T, all.y=F)
meta_data2 = meta_data2[match(meta_data$cell_id, meta_data2$cell_id),]
sum(meta_data$cell_id != colnames(merged_samples))
sum(meta_data2$cell_id != colnames(merged_samples))
meta_data2$strain = gsub(meta_data2$strain, pattern = 'rat_', replacement = '')
rownames(meta_data2) = meta_data2$cell_id
meta_data2 <- meta_data2[,colnames(meta_data2) != 'cell_id']
head(meta_data2)
table(meta_data2$strain)
###########################

head(merged_samples@meta.data) 
merged_samples@meta.data <- meta_data2
head(merged_samples@meta.data)

merged_samples <- readRDS('cell_browser/TLH_cellBrowser_2.rds')
colnames(merged_samples@meta.data)[1] = 'cluster_number'
colnames(merged_samples@meta.data)[ncol(merged_samples@meta.data)] = 'cluster'
Idents(merged_samples) = as.character(merged_samples$cluster)

saveRDS(merged_samples, 'cell_browser/TLH_cellBrowser_2.rds')
####################################################################
######## prepare the immune-enriched data for conversion ########
#merged_samples <- readRDS('Objects/merged_samples_newSamples_MT-removed.rds')
new_data_scCLustViz_object <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_labelCor.RData"
load(new_data_scCLustViz_object)
merged_samples = your_scRNAseq_data_object
Idents(merged_samples) = as.character(sCVdata_list$res.0.6@Clusters)
head(merged_samples)
merged_samples$cluster = as.character(sCVdata_list$res.0.6@Clusters)
merged_samples@meta.data[,grep(colnames(merged_samples@meta.data), pattern = 'res.')] = NULL
head(merged_samples@meta.data)
merged_samples@meta.data$orig.ident = NULL
head(merged_samples@meta.data)

### adding the cell-type annotation ###
annotation_Immune = read.csv('immune_enriched_annotation.csv')
meta_data = merged_samples@meta.data
meta_data$cell_id = rownames(meta_data)
meta_data2 = merge(meta_data, annotation_Immune, by.x='cluster', by.y='cluster', all.x=T, all.y=F)
meta_data2 = meta_data2[match(meta_data$cell_id, meta_data2$cell_id),]
sum(meta_data$cell_id != colnames(merged_samples))
sum(meta_data2$cell_id != colnames(merged_samples))
meta_data2$strain = gsub(meta_data2$strain, pattern = 'rat_', replacement = '')
rownames(meta_data2) = meta_data2$cell_id
meta_data2 <- meta_data2[,colnames(meta_data2) != 'cell_id']
table(meta_data2$strain)
#############
head(merged_samples@meta.data) 
merged_samples@meta.data <- meta_data2
head(merged_samples@meta.data)


merged_samples <- readRDS('cell_browser/Immune_enriched_cellBrowser_2.rds')
colnames(merged_samples@meta.data)[1] = 'cluster_number'
colnames(merged_samples@meta.data)[ncol(merged_samples@meta.data)] = 'cluster'
Idents(merged_samples) = as.character(merged_samples$cluster)
saveRDS(merged_samples, 'cell_browser/Immune_enriched_cellBrowser_2.rds')

################################################################
######## prepare the immune-subclustering data for conversion ########
new_data_scCLustViz_object_Immune <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_ImmuneSub_c17Included_labelCor.RData"
load(new_data_scCLustViz_object_Immune)
merged_samples = your_scRNAseq_data_object

Idents(merged_samples) = as.character(merged_samples$cluster)
merged_samples$cluster = as.character(sCVdata_list$res.1@Clusters)
head(merged_samples)

### adding the cell-type annotation ###
annotation_subImmune = read.csv('Immune_subclusters_annotation.csv')
meta_data = merged_samples@meta.data
meta_data$cell_id = rownames(meta_data)
meta_data2 = merge(meta_data, annotation_subImmune, by.x='cluster', by.y='cluster', all.x=T, all.y=F)
meta_data2 = meta_data2[match(meta_data$cell_id, meta_data2$cell_id),]
sum(meta_data$cell_id != colnames(merged_samples))
sum(meta_data2$cell_id != colnames(merged_samples))
meta_data2$strain = gsub(meta_data2$strain, pattern = 'rat_', replacement = '')
rownames(meta_data2) = meta_data2$cell_id
meta_data2 <- meta_data2[,colnames(meta_data2) != 'cell_id']
table(meta_data2$strain)
#############
head(merged_samples@meta.data) 
merged_samples@meta.data <- meta_data2
head(merged_samples@meta.data)

# nCount_RNA nFeature_RNA init_clusters immune_clusters mito_perc sample_name
merged_samples$orig.ident <- NULL
head(merged_samples)

merged_samples = readRDS('cell_browser/Immune_subclusters_cellBrowser_2.rds')
colnames(merged_samples@meta.data)[1] = 'cluster_number'
colnames(merged_samples@meta.data)[ncol(merged_samples@meta.data)] = 'cluster'
Idents(merged_samples) = as.character(merged_samples$cluster)

saveRDS(merged_samples, 'cell_browser/Immune_subclusters_cellBrowser_2.rds')

####################################################################
############## Checking the samples using visualization ###########
###################################################################

rm(list=ls())
merged_samples <- readRDS('cell_browser/TLH_cellBrowser.rds')
merged_samples <- readRDS('cell_browser/Immune_enriched_cellBrowser.rds')
merged_samples = readRDS('cell_browser/Immune_subclusters_cellBrowser.rds')

cluster_num = '4'
gene_name = 'Ptprc'
df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap')[,2], 
                      expression_val=GetAssayData(merged_samples)[gene_name,],
                      clusters=merged_samples$cluster)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=expression_val))+
  geom_point(alpha=0.6,size=1.2)+
  scale_color_viridis('Expression\nValue', direction = -1)+theme_classic()+
  theme(text = element_text(size=16.5),
        plot.title = element_text(hjust = 0.5),
        legend.title=element_text(size = 10))+
  ggtitle(paste0(gene_name)) 




