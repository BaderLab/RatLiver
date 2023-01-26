### loading the required libraries
source('~/RatLiver/Codes/Functions.R')
Initialize()

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


LIB_SIZE_CUT_OFF_MAX = 20000
LIB_SIZE_CUT_OFF = 1000
MIT_CUT_OFF = 10
NUM_GENES_DETECTED = 100

i = 4
sample_name = sample_names[i]


input_from_10x <- paste0("~/rat_sham_sn_data/", sample_name,'/filtered_feature_bc_matrix/')
data_raw <- CreateSeuratObject(counts=Read10X(input_from_10x, gene.column = 2),
                               min.cells=1,min.features=0, 
                               project = sample_name)

MIT_PATTERN = '^Mt-'
mito_genes_index <- grep(pattern = MIT_PATTERN, rownames(data_raw) )
rownames(data_raw)[mito_genes_index]
data_raw[["mito_perc"]] <- PercentageFeatureSet(data_raw, features = mito_genes_index)
summary(data_raw[["mito_perc"]]$mito_perc )
libSize <- colSums(GetAssayData(data_raw@assays$RNA))
dim(data_raw)



######## fast evaluation of the results 
print(paste0('Total number of cells: ', ncol(data_raw)))
to_drop_mito <- data_raw$mito_perc > MIT_CUT_OFF
print(paste0('to_drop_mito: ',sum(to_drop_mito)))
print(paste0('to_drop_mito percentage: ', round(sum(to_drop_mito)*100/ncol(data_raw),2) ))


to_drop_lib_size <- data_raw$nCount_RNA < LIB_SIZE_CUT_OFF | data_raw$nCount_RNA > LIB_SIZE_CUT_OFF_MAX
print(paste0('to_drop_lib_size: ', sum(to_drop_lib_size)))
print(paste0('to_drop_lib_size percentage: ', round( sum(to_drop_lib_size)*100/ncol(data_raw),2) ))

print(paste0('remained after both filters: ', sum(!to_drop_lib_size & !to_drop_mito)))
print(paste0('Percentage of remained after both filters: ', 
             round(  sum(!to_drop_lib_size & !to_drop_mito)*100/ncol(data_raw),2) ))

to_drop_num_genes <- data_raw$nFeature_RNA < NUM_GENES_DETECTED
print(paste0('to_drop_num_genes: ', sum(to_drop_num_genes)))
print(paste0('to_drop_num_genes percentage: ', round( sum(to_drop_num_genes)*100/ncol(data_raw),2) ))

print(paste0('remained after both filters: ', sum(!to_drop_num_genes & !to_drop_mito)))
print(paste0('Percentage of remained after both filters: ', 
             round(  sum(!to_drop_num_genes & !to_drop_mito)*100/ncol(data_raw),2) ))



############ importing dropletQC results
df = data.frame(library_size= data_raw$nCount_RNA, 
                mito_perc=data_raw$mito_perc , 
                n_expressed=data_raw$nFeature_RNA)

empty_drops.dc <- readRDS(paste0('~/rat_sham_sn_data/DropletQC_results/cellstatus_DF/',sample_name , '_empty_drop_DF.rds'))
sum(rownames(empty_drops.dc$df) != rownames(df))

df = cbind(df, empty_drops.dc$df)
data_raw$cell_status  = df$cell_status
data_raw$nuclear_fraction  = df$nuclear_fraction
data_raw$celltype_dQC  = df$cell_type

##### Properties of empty droplets based on DropletQC results
summary(df$library_size[df$cell_status == 'empty_droplet'])
summary(df$mito_perc[df$cell_status == 'empty_droplet'])
summary(df$n_expressed[df$cell_status == 'empty_droplet'])


## Visualization of QC metrics
#pdf(paste0('Plots/QC/',sample_name,'/QC_',sample_name,'_',
#           'mito_',MIT_CUT_OFF,'_lib_',LIB_SIZE_CUT_OFF,'.pdf'))

ggplot(data.frame(data_raw$nCount_RNA), aes(data_raw.nCount_RNA))+
  geom_histogram(bins = 60,color='black',fill='pink',alpha=0.5)+
  theme_bw()+ggtitle('library size for all cells (before filter)')+xlab('Library sizes')+
  ylab('Number of cells')+labs(caption = sample_name)

ggplot(data.frame(data_raw$nFeature_RNA), aes(data_raw.nFeature_RNA))+
  geom_histogram(bins = 60,color='black',fill='blue',alpha=0.3)+
  theme_bw()+ggtitle('# expressed genes for all cells(before filtering)')+xlab('Number of expressed genes')+
  ylab('Number of cells')+labs(caption = sample_name)

ggplot(data.frame(data_raw$mito_perc), aes(data_raw.mito_perc))+
  geom_histogram(bins = 60,color='black',fill='green',alpha=0.3)+
  theme_bw()+ggtitle('proportion of reads mapped to Mt genes(before filtering)')+xlab('Mitochondrial proportion (%)')+
  ylab('Number of cells')+labs(caption = sample_name)

ggplot(df, aes(x=library_size, y=mito_perc, color=n_expressed))+geom_point()+
  theme_bw()+xlab('Library Size')+ylab('Mitochondrial transcript percent')+
  geom_hline(yintercept= MIT_CUT_OFF, linetype="dashed", color = "red")+
  geom_vline(xintercept = LIB_SIZE_CUT_OFF, linetype="dashed", color = "red3")+
  scale_color_viridis('#expressed\ngenes',direction = +1)+
  ggtitle(paste0(sample_name,'\n','Mitochondrial transcript threshold: ', MIT_CUT_OFF,'\nLibrary size threshold: ', LIB_SIZE_CUT_OFF))


ggplot(df, aes(x=library_size, y=mito_perc, color=cell_status))+geom_point(size=1.5,alpha=0.7)+
  theme_classic()+xlab('Library size')+ylab('Mitochondrial transcript percent')+
  geom_hline(yintercept= MIT_CUT_OFF, linetype="dashed", color = "red", size=1)+
  geom_vline(xintercept = LIB_SIZE_CUT_OFF, linetype="dashed", color = "red3", size=1)+
  #ggtitle('Library size and mitochondrial transcript cutoffs')+
  theme(axis.title = element_text(size = 17), axis.text = element_text(size = 16), plot.title=element_text(size=18))


ggplot(df, aes(x=n_expressed, y=mito_perc, color=library_size))+geom_point()+scale_color_viridis()+
  theme_bw()+xlab('num detected genes')+ylab('Mitochondrial transcript percent')+
  geom_hline(yintercept= MIT_CUT_OFF, linetype="dashed", color = "red")+labs(caption = sample_name)+
  geom_vline(xintercept = NUM_GENES_DETECTED, linetype="dashed", color = "red3", size=0.5)+
  ggtitle(paste0('mito threshold: ', MIT_CUT_OFF,' , num detected genes threshold: ', NUM_GENES_DETECTED, ' (before filter)'))


ggplot(df, aes(x=library_size, y=n_expressed, color=mito_perc))+geom_point()+labs(caption = sample_name)+
  theme_bw()+xlab('Library Size')+ylab('Number of expressed genes')+scale_color_viridis(option = 'magma')+ggtitle('before filter')





data <- data_raw[,!to_drop_mito & !to_drop_lib_size]
show(data)

df_filt = data.frame(library_size= data$nCount_RNA, 
                     mito_perc=data$mito_perc , 
                     cell_status=data$cell_status, 
                     nuclear_fraction=data$nuclear_fraction,
                     n_expressed= data$nFeature_RNA)

ggplot(df_filt, aes(x=library_size, y=mito_perc, color=library_size))+geom_point()+scale_color_viridis()+
  theme_bw()+xlab('Library Size')+ylab('Mitochondrial transcript percent')+labs(caption = sample_name)+
  ggtitle(paste0('mito threshold: ', MIT_CUT_OFF,' , library size threshold: ', LIB_SIZE_CUT_OFF, ' (after filter)'))

ggplot(df_filt, aes(x=library_size, y=mito_perc, color=cell_status))+geom_point()+
  theme_bw()+xlab('Library Size')+ylab('Mitochondrial transcript percent')+labs(caption = sample_name)+
  ggtitle(paste0('mito threshold: ', MIT_CUT_OFF,' , library size threshold: ', LIB_SIZE_CUT_OFF, ' (after filter)'))

ggplot(df_filt, aes(x=n_expressed, y=mito_perc, color=library_size))+geom_point()+scale_color_viridis()+
  theme_bw()+xlab('Number of detected genes')+ylab('Mitochondrial transcript percent')+labs(caption = sample_name)+
  ggtitle(paste0('mito threshold: ', MIT_CUT_OFF,' , num detected genes threshold: ', LIB_SIZE_CUT_OFF, ' (after filter)'))






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
                      mito_perc = data$mito_perc,
                      n_expressed=data$nFeature_RNA,
                      cluster=data$seurat_clusters, 
                      cell_status = data$cell_status, 
                      nuclear_fraction = data$nuclear_fraction, 
                      Alb=GetAssayData(data)['Alb',])


ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=nuclear_fraction))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(sample_name)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=library_size))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(sample_name)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cell_status))+geom_point(alpha=0.3)+theme_classic()+ggtitle(sample_name)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=Alb))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(sample_name)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=mito_perc))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(sample_name)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cluster))+geom_point(alpha=0.3)+theme_classic()+ggtitle(sample_name)
saveRDS(data, paste0('~/rat_sham_sn_data/standardQC_results/',sample_name , '_data_afterQC.rds'))



############################################################
############## Integrating the filtered files ##############
########################################################

list_files = list.files(path = '~/rat_sham_sn_data/standardQC_results', pattern ='_data_afterQC.rds', include.dirs = T, full.names = T)
files_rds <- lapply(list_files, readRDS)

sample_names = list.files(path = '~/rat_sham_sn_data/standardQC_results', pattern ='_data_afterQC.rds')
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

merged_data <- merged_data_backup


### run the pipeline without this line and see how it would differ - could not run pca and find variable genes
#### scaling the merged data instead of SCTransform does not solve the issue
### current pipeline: we are normlizing the individual samples first - merging them - normalizing the merged sample 
merged_data = SCTransform(merged_data, vst.flavor = "v2", verbose = TRUE) 
#merged_data <- ScaleData(merged_data)
#merged_data <- FindVariableFeatures(merged_data) #SCT assay is comprised of multiple SCT models

merged_data <- RunPCA(merged_data,verbose=T)
plot(100 * merged_data@reductions$pca@stdev^2 / merged_data@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")

#### UMAP on PC components and clustering
merged_data <- RunUMAP(merged_data, reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

#### UMAP on corrected PC components and clustering
merged_data <- RunHarmony(merged_data, group.by.vars = "sample_name")
merged_data <- RunUMAP(merged_data, reduction = "harmony", dims = 1:30, reduction.name = "umap_h")  %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

Resolution = 2.2
merged_data <- FindClusters(merged_data, resolution = Resolution, verbose = FALSE)
gene_name = 'Cd5l'

df_umap <- data.frame(UMAP_1=getEmb(merged_data, 'umap')[,1], 
                      UMAP_2=getEmb(merged_data, 'umap')[,2], 
                      UMAP_h_1=getEmb(merged_data, 'umap_h')[,1], 
                      UMAP_h_2=getEmb(merged_data, 'umap_h')[,2], 
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


ggplot(df_umap, aes(x=UMAP_h_1, y=UMAP_h_2, color=cluster))+geom_point(alpha=0.6)+theme_classic()+ggtitle(paste0('res: ', Resolution))

ggplot(df_umap, aes(x=UMAP_h_1, y=UMAP_h_2, color=gene))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(gene_name)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=gene))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(gene_name)

ggplot(df_umap, aes(x=UMAP_h_1, y=UMAP_h_2, color=library_size))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)
ggplot(df_umap, aes(x=UMAP_h_1, y=UMAP_h_2, color=mito_perc))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(sample_name)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=nuclear_fraction))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)


ggplot(df_umap, aes(x=UMAP_h_1, y=UMAP_h_2, color=cell_status))+geom_point(alpha=0.3)+theme_classic()
ggplot(df_umap, aes(x=UMAP_h_1, y=UMAP_h_2, color=strain))+geom_point(alpha=0.3)+theme_classic()
ggplot(df_umap, aes(x=UMAP_h_1, y=UMAP_h_2, color=sample_name))+geom_point(alpha=0.3)+theme_classic()

saveRDS(merged_data, '~/rat_sham_sn_data/standardQC_results/sham_sn_merged_standardQC.rds')

