source('Codes/Functions.R')
Initialize()
library(remotes)
library(ggplot2)
library(RColorBrewer)
library(stats)
library(ggpubr)

library(celda)
library(singleCellTK)
library(SingleCellExperiment)


######### importing the old samples raw data
## Define cut-off values
qc_info = list('MacParland_Sonya__RAT_SHAM_DA_M_10WK_003_strained'=list(LIB_SIZE_CUT_OFF=2000, MIT_CUT_OFF=20), 
               'McParland_Sonya__RAT_SHAM_DA_M_10WK_004_strained'=list(LIB_SIZE_CUT_OFF=1500, MIT_CUT_OFF=30),
               'McParland_Sonya__RAT_SHAM_LEW_M_13WK_001_strained'=list(LIB_SIZE_CUT_OFF = 2000,MIT_CUT_OFF = 40),
               'McParland_Sonya__RAT_SHAM_LEW_M_13WK_002_strained'=list(LIB_SIZE_CUT_OFF = 2000,MIT_CUT_OFF = 40))
data_norm_list = list(NA, NA, NA, NA)
names(data_norm_list) = sample_names
strain_info = c('DA', 'DA', 'LEW', 'LEW')

##############################
### making a list of paths to the filtered and raw matrices
# note: DA_M_10WK_004_strained is the 200515_A00827_0155_AHLV55DRXX_McParland_Sonya reseq sample
sample_names = list.dirs('~/RatLiver/Data/SoupX_data/SoupX_inputs/',recursive = FALSE, full.names = FALSE)
path_vector = rep(x =NA ,length(sample_names))
#path_vector_raw = rep(x =NA ,length(sample_names))


genes_df <- read.delim(paste0('Data/rat_Lew_01/features.tsv.gz'), header = F)
head(genes_df)

for(i in 1:length(sample_names)){
  sample_name = sample_names[i]
  path_vector[i] = paste0('~/RatLiver/Data/SoupX_data/SoupX_inputs/', sample_name)
  #path_vector_raw[i] =paste0('~/RatLiver/Data/SoupX_data/SoupX_inputs/', sample_name)
}

for(i in 1:length(sample_names)){
  sce <- importCellRanger(sampleDirs = path_vector[i], dataType='filtered')
  sce.raw <- importCellRanger(sampleDirs = path_vector[i], dataType='raw')
  sce <- decontX(sce, background = sce.raw)
  #saveRDS(sce, paste0('~/RatLiver/Data/DecontX_data/', sample_names[i], '_decontX_out.rds'))
}



### next to do
# check the umap of each sample colored based on contamination and Alb
# convert gene ids to gene names for query and mt-based QC
## merge the decont matrices and represent the total map similar to soupX results


i = 4
for(i in 1:length(sample_names)){
  sample_name = sample_names[i]
  sample_name
  sce <- readRDS(paste0('~/RatLiver/Data/DecontX_data/', sample_names[i], '_decontX_out.rds'))
  plotDecontXContamination(sce)
  
  #rownames(soupX_out$out) = gsub(rownames(soupX_out$out),pattern = '_', replacement = '-')
  data = CreateSeuratObject(counts = assay(sce, 'decontXcounts'))
  data$contamination_frac = sce$decontX_contamination
  data$sample_name = sample_name
  data$strain = strain_info[i]
  MIT_PATTERN = '^Mt-'
  mito_genes_index <- grep(pattern = MIT_PATTERN, genes_df$V2 )
  data[["mito_perc"]] <- PercentageFeatureSet(data, features = mito_genes_index)
  
  to_drop_mito = data$mito_perc > qc_info[[sample_name]]$MIT_CUT_OFF
  to_drop_lib_size = data$nCount_RNA < qc_info[[sample_name]]$LIB_SIZE_CUT_OFF
  to_drop_numgenes = data$nFeature_RNA < 100 
  
  summary(data$mito_perc )
  summary(data$nFeature_RNA )
  sum(!to_drop_mito & !to_drop_numgenes)
  sum(!to_drop_mito & !to_drop_lib_size)
  sum(to_drop_lib_size)
  sum(to_drop_numgenes)
  sum(to_drop_mito)
  
  data <- data[,!to_drop_mito & !to_drop_numgenes]
  data = SCTransform(data, vst.flavor = "v2", verbose = TRUE)
  
  data_norm_list[[i]] = data
}
rm(data)
gc()

merged_data <- merge(data_norm_list[[1]], c(data_norm_list[[2]], data_norm_list[[3]], data_norm_list[[4]] ), # 
                     add.cell.ids = names(data_norm_list), 
                     project = names(data_norm_list), 
                     merge.data = TRUE)
dim(merged_data)
Idents(merged_data) = merged_data$sample_name


rm(data_norm_list)
gc()
saveRDS(merged_data, '~/RatLiver/Data/DecontX_data//TLH_merged_decontX_decontaminated.rds')


merged_data <- readRDS( '~/RatLiver/Data/DecontX_data//TLH_merged_decontX_decontaminated.rds')
DefaultAssay(merged_data) <- 'SCT'
merged_data@assays$RNA <- NULL
gc()

merged_data = SCTransform(merged_data, vst.flavor = "v2", verbose = TRUE, assay = 'SCT', conserve.memory = TRUE) 

#saveRDS(merged_data, '~/RatLiver/Data/DecontX_data/TLH_merged_decontX_decontaminated_normed.rds')
merged_data = readRDS('~/RatLiver/Data/DecontX_data/TLH_merged_decontX_decontaminated_normed.rds')

merged_data <- RunPCA(merged_data,verbose=T)
plot(100 * merged_data@reductions$pca@stdev^2 / merged_data@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")

#### UMAP on corrected PC components and clustering
merged_data <- RunHarmony(merged_data, group.by.vars = "sample_name")
merged_data <- RunUMAP(merged_data, reduction = "harmony", dims = 1:30, reduction.name = "umap_h")  %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)


get_ens_id <- function(gene_symbol){
  genes_df$V1[genes_df$V2 == gene_symbol]
}

Resolution = 0.2
merged_data <- FindClusters(merged_data, resolution = Resolution, verbose = FALSE)
gene_name = 'Cd68'

is_cluster = ifelse(merged_data$seurat_clusters== 7, 'cluster', 'other')
df_umap <- data.frame(#UMAP_1=getEmb(merged_data, 'umap')[,1], 
  #UMAP_2=getEmb(merged_data, 'umap')[,2], 
  UMAP_h_1=getEmb(merged_data, 'umap_h')[,1], 
  UMAP_h_2=getEmb(merged_data, 'umap_h')[,2], 
  library_size= merged_data$nCount_RNA, 
  mito_perc=merged_data$mito_perc, 
  n_expressed=merged_data$nFeature_RNA,
  cluster=merged_data$seurat_clusters, 
  is_cluster=is_cluster,
  #cell_status = merged_data$cell_status,
  #nuclear_fraction=merged_data$nuclear_fraction, 
  Alb=GetAssayData(merged_data)[get_ens_id('Alb'),], 
  gene=GetAssayData(merged_data)[get_ens_id(gene_name),],
  sample_name = merged_data$sample_name, 
  strain = merged_data$strain)

ggplot(df_umap, aes(x=UMAP_h_1, y=UMAP_h_2, color=is_cluster))+geom_point(alpha=0.8)+theme_classic()#+ggtitle(paste0('res: ', Resolution))
ggplot(df_umap, aes(x=UMAP_h_1, y=UMAP_h_2, color=cluster))+geom_point(alpha=0.8)+theme_classic()#+ggtitle(paste0('res: ', Resolution))
ggplot(df_umap, aes(x=UMAP_h_1, y=UMAP_h_2, color=sample_name))+geom_point(alpha=0.3)+theme_classic()#+ggtitle(paste0('res: ', Resolution))

ggplot(df_umap, aes(x=UMAP_h_1, y=UMAP_h_2, color=gene))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(gene_name)+xlab('UMAP_1')+ylab('UMAP_2')
ggplot(df_umap, aes(x=UMAP_h_1, y=UMAP_h_2, color=Alb))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle('Alb')+xlab('UMAP_1')+ylab('UMAP_2')


genes = c('Alb', 'Tat', 'G6pc', 'Cps1', 'Tdo2')
genes = c('Lyve1', 'Id3') # Lsecs
genes = c('Igfbp7', 'Rbp1', 'Col3a1', 'Sparc') # stellate cells
genes = c('Nkg7', 'Cd7')

Idents(merged_data) = as.character(merged_data$seurat_clusters)
cluster7_markers = FindMarkers(merged_data, ident.1 = '7')
cluster7_markers <- cluster7_markers[order(cluster7_markers$avg_log2FC, decreasing = T),]
cluster7_markers <- cluster7_markers[order(cluster7_markers$p_val_adj, decreasing = F),]


