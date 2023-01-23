plan("multiprocess", workers = 3)
options(error = function() {traceback(2, max.lines=100); if(!interactive()) quit(save="no", status=1, runLast=T)})
options(future.globals.maxSize = 200000 * 1024^2)

library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)
library(ggplot2)
library(DropletQC)
library(ggplot2)
library(patchwork)
library(dplyr)


sample_names = c('MacParland__SingleNuc_DA_SHAM003_CST_3pr_v3', 'MacParland__SingleNuc_DA_SHAM004_CST_3pr_v3' , 
                 'MacParland__SingleNuc_LEW_SHAM1_CST_3pr_v3', 'MacParland__SingleNuc_LEW_SHAM2_CST_3pr_v3')


###########################################################################
############ Calculating nuclear fraction using droplet QC #############
#devtools::install_github("powellgenomicslab/DropletQC", build_vignettes = FALSE)
main_dir = '~/SynologyDrive/211102_A00827_0441_BH27VGDMXY_MacParland/'

for(i in 2:length(sample_names)){
  sample_name = sample_names[i]
  outs = paste0(main_dir, sample_name, '/', 'outs/')
  
  nf <- nuclear_fraction_tags(outs=outs,verbose = TRUE)
  saveRDS(nf, paste0('~/rat_sham_sn_data/',sample_name , '_nucfrac.rds'))
}

###########################################################################
# ## import the data



### loading the required libraries
source('~/RatLiver/Codes/Functions.R')
Initialize()

i = 4
sample_name = sample_names[i]


input_from_10x <- paste0("~/rat_sham_sn_data/", sample_name,'/filtered_feature_bc_matrix/')
data_raw <- CreateSeuratObject(counts=Read10X(input_from_10x, gene.column = 2),
                               min.cells=0,min.features=1, 
                               project = sample_name)
dim(data_raw)
libSize <- colSums(GetAssayData(data_raw@assays$RNA))

data_raw = SCTransform(data_raw, vst.flavor = "v2", verbose = TRUE)
data_raw <- RunPCA(data_raw,verbose=T)
plot(100 * data_raw@reductions$pca@stdev^2 / data_raw@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")

data_raw <- RunUMAP(data_raw, reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

nf = readRDS(paste0('~/rat_sham_sn_data/',sample_name , '_nucfrac.rds'))
nf.umi.df <- data.frame(nf=nf, umi=data_raw$nCount_RNA)
empty_drops.df <- identify_empty_drops(nf_umi=nf.umi.df, include_plot = TRUE)
#cluster the data and then identify the damaged cells
empty_drops.df$cell_type = data_raw$seurat_clusters
empty_drops.dc <- identify_damaged_cells(empty_drops.df, verbose = TRUE, output_plots = TRUE)

data_raw$cell_status = empty_drops.dc[[1]]$cell_status
length(empty_drops.dc)
table(empty_drops.dc[[1]]$cell_status)
wrap_plots(empty_drops.dc[[2]], nrow = 5)

df_umap <- data.frame(UMAP_1=getEmb(data_raw, 'umap')[,1], 
                      UMAP_2=getEmb(data_raw, 'umap')[,2], 
                      library_size= data_raw$nCount_RNA, 
                      n_expressed=data_raw$nFeature_RNA,
                      cluster=data_raw$seurat_clusters, 
                      cell_status = empty_drops.dc[[1]]$cell_status, 
                      nuclear_fraction = empty_drops.dc[[1]]$nuclear_fraction, 
                      Alb=GetAssayData(data_raw)['Alb',])


ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=nuclear_fraction))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(sample_name)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=library_size))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(sample_name)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cell_status))+geom_point(alpha=0.3)+theme_classic()+ggtitle(sample_name)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=Alb))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(sample_name)

saveRDS(empty_drops.dc, paste0('~/rat_sham_sn_data/',sample_name , '_empty_drop_DF.rds'))
empty_drops.dc <- readRDS(paste0('~/rat_sham_sn_data/',sample_name , '_empty_drop_DF.rds'))



#####################################################
data <- CreateSeuratObject(counts=Read10X(input_from_10x, gene.column = 2),
                               min.cells=0,min.features=1, 
                               project = sample_name)
dim(data)
libSize <- colSums(GetAssayData(data@assays$RNA))


genes_df <- read.delim(paste0(input_from_10x,'features.tsv.gz'), header = F)
data[['RNA']] <- AddMetaData(data[['RNA']], genes_df$V2, col.name = 'symbol')
data[['RNA']] <- AddMetaData(data[['RNA']], genes_df$V1, col.name = 'ensembl')

MIT_PATTERN = '^Mt-'
mito_genes_index <- grep(pattern = MIT_PATTERN, data[['RNA']]@meta.features$symbol )
data[['RNA']]@meta.features$symbol[mito_genes_index]
data[["mito_perc"]] <- PercentageFeatureSet(data, features = mito_genes_index)
summary(data[["mito_perc"]]$mito_perc )

data$cell_status = empty_drops.dc[[1]]$cell_status
data$nuclear_fraction = empty_drops.dc[[1]]$nuclear_fraction
sum(rownames(empty_drops.dc[[1]]) != colnames(data))
table(data$cell_status)
hist(data$mito_perc)
sum(data$mito_perc >= 10)


data = data[,data$cell_status != 'empty_droplet']
data = data[,data$mito_perc < 10]
dim(data)

data = SCTransform(data, vst.flavor = "v2", verbose = TRUE)
data <- RunPCA(data,verbose=T)
plot(100 * data@reductions$pca@stdev^2 / data@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")

data <- RunUMAP(data, reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

df_umap <- data.frame(UMAP_1=getEmb(data, 'umap')[,1], 
                      UMAP_2=getEmb(data, 'umap')[,2], 
                      library_size= data$nCount_RNA, 
                      mito_perc=data$mito_perc, 
                      n_expressed=data$nFeature_RNA,
                      cluster=data$seurat_clusters, 
                      cell_status = data$cell_status,
                      nuclear_fraction=data$nuclear_fraction, 
                      Alb=GetAssayData(data)['Alb',])


ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=library_size))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(sample_name)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=mito_perc))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(sample_name)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cell_status))+geom_point(alpha=0.3)+theme_classic()+ggtitle(sample_name)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cluster))+geom_point(alpha=0.3)+theme_classic()+ggtitle(sample_name)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=Alb))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(sample_name)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=nuclear_fraction))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(sample_name)


saveRDS(data, paste0('~/rat_sham_sn_data/',sample_name , '_data_afterQC.rds'))

############################################################
############## Integrating the filtered files ##############
########################################################

list_files = list.files(path = '~/rat_sham_sn_data/', pattern ='_data_afterQC.rds', include.dirs = T, full.names = T)
files_rds <- lapply(list_files, readRDS)

sample_names = list.files(path = '~/rat_sham_sn_data/', pattern ='_data_afterQC.rds')
sample_names = gsub('_3pr_v3_data_afterQC.rds', '', sample_names)
sample_names = gsub('MacParland__SingleNuc_', '', sample_names)
strain = c('DA', 'DA', 'LEW', 'LEW')

files_rds = sapply(1:length(files_rds), 
                       function(i) {
                         files_rds[[i]]$sample_name = sample_names[i]
                         files_rds[[i]]$strain = strain[i]
                         return(files_rds[[i]])
                       },simplify = F)

names(files_rds) = sample_names
head(files_rds[[1]])
lapply(files_rds, dim)

merged_data <- merge(files_rds[[1]], c(files_rds[[2]], files_rds[[3]], files_rds[[4]] ), # 
                  add.cell.ids = names(files_rds), 
                  project = names(files_rds), 
                  merge.data = TRUE)
dim(merged_data)

Idents(merged_data) = merged_data$sample_name


merged_data = SCTransform(merged_data, vst.flavor = "v2", verbose = TRUE) ### run the pipeline without this line and see how it would differ - could not run pca and find variable genes
#merged_data <- FindVariableFeatures(merged_data) #SCT assay is comprised of multiple SCT models
merged_data <- RunPCA(merged_data,verbose=T)
plot(100 * merged_data@reductions$pca@stdev^2 / merged_data@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")

#### UMAP on PC components and clustering
merged_data <- RunUMAP(merged_data, reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

#### UMAP on corrected PC components and clustering
merged_data <- RunHarmony(merged_data, group.by.vars = "sample_name", assay.use="RNA")
merged_data <- RunUMAP(merged_data, reduction = "harmony", dims = 1:30)  %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resoution = 0.7, verbose = FALSE)


gene_name = 'Itgal'

df_umap <- data.frame(UMAP_1=getEmb(merged_data, 'umap')[,1], 
                      UMAP_2=getEmb(merged_data, 'umap')[,2], 
                      library_size= merged_data$nCount_RNA, 
                      mito_perc=merged_data$mito_perc, 
                      n_expressed=merged_data$nFeature_RNA,
                      cluster=merged_data$seurat_clusters, 
                      cell_status = merged_data$cell_status,
                      nuclear_fraction=merged_data$nuclear_fraction, 
                      Alb=GetAssayData(merged_data)['Alb',], 
                      sample_name = merged_data$sample_name, 
                      strain = merged_data$strain, 
                      gene=GetAssayData(merged_data)[gene_name,])


ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=library_size))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=mito_perc))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(sample_name)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=gene))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(gene_name)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cluster))+geom_point(alpha=0.3)+theme_classic()
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cell_status))+geom_point(alpha=0.3)+theme_classic()
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=strain))+geom_point(alpha=0.3)+theme_classic()
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=sample_name))+geom_point(alpha=0.3)+theme_classic()



saveRDS(merged_data, '~/rat_sham_sn_data/sham_singleNuc_merged.rds')



