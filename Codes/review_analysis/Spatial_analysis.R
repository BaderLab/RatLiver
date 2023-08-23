library(Seurat)
library(Seuratobj)
library(ggplot2)
library(patchwork)
library(dplyr)
require(scales)

get_scaled_by_gene <- function(x){
  y <- scale(t(x), center = TRUE, scale = TRUE)
  y[is.nan(y)] <- 0
  y2 <- t(y)
  y2 <- as.data.frame(y2)
  return(y2)
}



get_varimax_rotated <- function(gene_exp_matrix, loading_matrix){
  
  ## gene_exp_matrix: gene expression matrix. rows named as genes and columns named as UMIs.
  ##                  attention: Data SHOULD be SCALED.
  ##                  Potentially gained from Seurat GetAssayData() function
  ## loading_matrix: the PCA loading matrix with rows named as gene names and 
  ##                 columns named as PC numbers. Potentially gained from Seurat Loadings() function
  
  ######### Varimax rotation
  initial_data <- t(gene_exp_matrix[rownames(gene_exp_matrix) %in% rownames(loading_matrix),])
  
  ## apply varimax rotation on the loadings
  varimax_res <- varimax(loading_matrix)
  rotatedLoadings <- varimax_res$loadings
  class(rotatedLoadings) <- "matrix";
  ## calculating the PC scores matrix
  invLoadings     <- t(pracma::pinv(rotatedLoadings))
  scores          <- scale(initial_data) %*% invLoadings ## this second scaling is not necessary
  
  ## compacting the rotated loading and score matrices in a list
  rotated_data <- list(rotLoadings=rotatedLoadings, rotScores = scores)
  return(rotated_data)
}


#### Tallula's spatial pipelinne:
#https://github.com/tallulandrews/SpatialAnalysis/blob/main/Spatial_Script.R
set.seed(101)

file1 = 'MacParland_Catia__Rat_6-1_A1_VIS'
slice <- file1

file2 = 'MacParland_Catia__Rat_6-1_B1_VIS'
slice <- file2

obj <- Load10X_Spatial(paste('~/Spatial_Rat/', slice, sep=""), 
                       filename="filtered_feature_bc_matrix.h5", 
                       slice=slice)

obj@meta.obj$orig.ident <- rep(slice, ncol(obj))
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^Mt-")


# QC
plot_count <- SpatialFeaturePlot(obj, features = "nCount_Spatial") + 
  theme(legend.position = "right")
plot_count
plot_feature <- SpatialFeaturePlot(obj, features = "nFeature_Spatial") + 
  theme(legend.position = "right")
plot_feature
plot_count + plot_feature


plot_mt <- SpatialFeaturePlot(obj, features = "percent.mt") + 
  theme(legend.position = "right")
plot_mt


# FeatureSelection
# Filter
detect_freq <- Matrix::rowMeans(obj@assays$Spatial@counts > 0)
highest <- apply(obj@assays$Spatial@counts, 1, max)
# genes with detection frequency of less than 0.05 and maximum read capture of less than 3 were removed during the quality control process. 
obj <- obj[detect_freq > 0.05 & highest >= 3,]
obj <- obj[!grepl("^Mt-", rownames(obj)),]

FindSpatiallyVariableFeatures(object = obj, slot = "scale.data", nfeatures=500)

# Normalization
obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE)

mean_expr <- Matrix::rowMeans(obj@assays$Spatial@data)


########## 
obj <- FindSpatiallyVariableFeatures(obj, assay = "SCT", features = VariableFeatures(obj)[1:1000],
                                       selection.method = "moransi")

top.features <- head(SpatiallyVariableFeatures(obj, selection.method = "moransi"), 500)


obj = readRDS('~/Spatial_Rat/Spatial_MacParland_Catia__Rat_6-1_A1_seurObj.rds')
dim(obj@assays$Spatial)


obj = readRDS('~/Spatial_Rat/Spatial_MacParland_Catia__Rat_6-1_B1_seurObj.rds')
colnames(obj)

################################


PP1_list= c("Sds", "Hal", "Cyp2f4", "Arg1", "Acly", "Ass1", "Gls2", "Agxt", 
            "Uroc1", "Gldc", "Gls2")
PP2_list= c("Apoc3", "Serpina1", "Apoc1", "Apoe", "Itih4", "Apoa1", "Ttr", 
            "Tf", "Alb", "Pigr", "Orm1", "Rpl3", "Fads1", "Aldh1b1", "Srd5a1", "Hsd17b13")
IZ_list= c("Saa4", "Hint1", "Cyp8b1", "Cyp2e1", "Cox7c", "Fabp1")
CV_list= c("Notum", "Cyp27a1", "Fabp7", "Akr1c1", "Gsta5", "Slc22a1", "Aox3",
           "Sult1e1", "Fmo1", "Oat", "Ahr", "Cyp7a1", "Glul", "Rhbg", "Cyp2e1", "Cyp1a2", "Por")
markers_list = c(CV_list, IZ_list, PP2_list, PP1_list)

################################
PC_markers = c("Glul", "Cyp27a1", "Akr1c2", "Nr1i3", "Csad", "Sult1e1", "Slc22a1", "Avpr1a", 
            "Fmo1", "Ahr", "Rgn", "Sord", "Oat", "Akr7a3", "Mup4", "Dhrs7l1", "Slc13a3", 
            "Cela1", "Cyp7a1", "Gstm1", "Slc27a5", "LOC100360095", "LOC100910877", "LOC259244", 
            "AABR07048474.1", "AABR07048463.1", "AABR07047899.1", "Cyp2e1",
            "C1r", "Cyp2c11", "Cyp1a2", "Cdh17", "Xpnpep2", "Cyp2c22", "Gstt3", "Fetub",
            "Gulo", "Ugt1a5", "Ces1d", "Ca3", "Scd")

PP_markers = c("Apoa4", "Orm1", "Hmgcs1", "Cdh1", "Apoc2", "Cyp4a2", "Atp1a1", "Wfdc2", "Srd5a1", 
              "Spink1l", "Itih4", "Wfdc21", "Gjb2", "Apoc3", "Gpx1", "Apoa1", "Alb", "Rarres2", 
              "Cfi", "Fabp1", "Ctsh", "Uroc1", "Etnppl", "Hsd17b13", "Slc38a4", "Sult2a6", "Ass1", 
              "Gldc", "Arg1", "Agxt", "Cyp2c7", "Crot", "Ftcd", "Sds", "Hal", "Mfsd2a", "Cyp4a2",
              "Cps1", "Cyp2c7", "Tdo2", "Slc25a47", "Gls2", "Bhmt")

markers_list = c(PC_markers,PP_markers)
################################

pdf('~/Spatial_Rat/Spatial_MacParland_Catia__Rat_6-1_A1_VIS_HepMarkers_2.pdf')
pdf('~/Spatial_Rat/Spatial_MacParland_Catia__Rat_6-1_B1_VIS_HepMarkers_2.pdf')
markers_list = markers_list[markers_list%in%rownames(obj)]

for(i in 1:length(markers_list)){
  p=SpatialFeaturePlot(obj, features = markers_list[i])+theme(legend.text = element_text(size = 8),
            legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"))
  print(p)
}
dev.off()
marker = "Glul"
marker = "Cd163" # Marco Cd68
SpatialFeaturePlot(obj, features = marker)+
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 20), 
        legend.key.size = unit(1, "cm"))

# Pipeline
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:20)
obj <- FindClusters(obj, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:20)

p1 <- DimPlot(obj, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(obj, label = TRUE, label.size = 3)

p1
p2


loading_mat <- read.csv(paste0('~/Spatial_Rat/', slice, '_pca_loadings.csv') )
head(loading_mat)
rownames(loading_mat) = (loading_mat$X)
loading_mat <- loading_mat[,-1]

npcs = 12
pca_scores = data.frame(obj@reductions$pca@cell.embeddings)[,1:npcs]
loading_mat = data.frame(obj@reductions$pca@feature.loadings)[,1:npcs]
head(loading_mat[order(loading_mat$PC_1, decreasing = F),],30)
#write.csv(loading_mat,paste0('~/Spatial_Rat/', slice, '_pca_loadings.csv'), row.names = T)
                         
                         
pericentral_spatial = c('Glul', 'Slc22a1', 'Fmo1', 'Oat', 'Cyp7a1', 'Cyp27a1', 'Sul1e1', 'Cyp1a2')
periportal_spatial = c('Aldh1b1', 'Srd5a1', 'Hsd17b13', 'Sds', 'Agxt', 'Gsl2', 'Ass1', 'Uroc1', 'Orm1')

# spatial plot of pca scores

obj@meta.data$PC1_score <- obj@reductions$pca@cell.embeddings[,1]
obj@meta.data$PC2_score <- obj@reductions$pca@cell.embeddings[,2]
obj@meta.data$PC3_score <- obj@reductions$pca@cell.embeddings[,3]
obj@meta.data$PC4_score <- obj@reductions$pca@cell.embeddings[,4]
obj@meta.data$PC5_score <- obj@reductions$pca@cell.embeddings[,5]
obj@meta.data$PC6_score <- obj@reductions$pca@cell.embeddings[,6]
pca_plot1 <- SpatialFeaturePlot(obj, features = c("PC1_score","PC2_score", "PC3_score", 
                                                  "PC4_score", "PC5_score", "PC6_score"))
pca_plot1

obj@meta.data$PC7_score <- obj@reductions$pca@cell.embeddings[,7]
obj@meta.data$PC8_score <- obj@reductions$pca@cell.embeddings[,8]
obj@meta.data$PC9_score <- obj@reductions$pca@cell.embeddings[,9]
obj@meta.data$PC10_score <- obj@reductions$pca@cell.embeddings[,10]
obj@meta.data$PC11_score <- obj@reductions$pca@cell.embeddings[,11]
obj@meta.data$PC12_score <- obj@reductions$pca@cell.embeddings[,12]
pca_plot2 <- SpatialFeaturePlot(obj, features = c("PC7_score","PC8_score", "PC9_score", 
                                                  "PC10_score", "PC11_score", "PC12_score"), )
pca_plot2

SpatialFeaturePlot(obj, features = 'PC1_score')+
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 20), 
        legend.key.size = unit(1, "cm"))

SpatialFeaturePlot(obj, features = 'PC2_score')+
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 20), 
        legend.key.size = unit(1, "cm"))

pc_num = 1
PC_df = data.frame(genes=rownames(loading_mat), loading=loading_mat[,pc_num])
PC_df = PC_df[order(PC_df$loading, decreasing =F),]
head(PC_df, 30)

pericentral_spatial_data = t(GetAssayData(obj)[rownames(obj) %in% pericentral_spatial,])
sum(rownames(pericentral_spatial_data) != rownames(pca_scores))
cor_mat = cor(cbind(pca_scores,pericentral_spatial_data))
colnames(cor_mat)
pheatmap(cor_mat[1:ncol(pca_scores),(ncol(pca_scores)+1):ncol(cor_mat) ])

head(pca_scores)
dim(pca_scores)
dim(loading_mat)



##############################################################
############### Varimax analysis ######################
##############################################################
pca_loading <- obj@reductions$pca@feature.loadings[,1:npcs]
pca_loading <- pca_loading/rowSums(pca_loading)
test_pca <- varimax(pca_loading)


gene_loadings_matrix <- c();
rotated_pca <- data.frame(obj@reductions$pca@cell.embeddings[,1:npcs] %*% test_pca$rotmat)
loading_mat <- data.frame(obj@reductions$pca@feature.loadings[,1:npcs] %*% test_pca$rotmat)
head(loading_mat)
colnames(loading_mat) <- paste("RotPC", 1:ncol(loading_mat), sep="_")

if (length(gene_loadings_matrix) == 0) {
  gene_loadings_matrix <- loading_mat
} else {
  all_genes <- sort(unique(c(rownames(gene_loadings_matrix), rownames(loading_mat))))
  gene_loadings_matrix <- gene_loadings_matrix[match(all_genes, rownames(gene_loadings_matrix)),]
  loading_mat <- loading_mat[match(all_genes, rownames(loading_mat)),]
  loading_mat[is.na(loading_mat)] <- 0
  gene_loadings_matrix[is.na(gene_loadings_matrix)] <- 0;
  rownames(gene_loadings_matrix) <- all_genes
  rownames(loading_mat) <- all_genes
  gene_loadings_matrix <- cbind(gene_loadings_matrix, loading_mat)
}


colnames(rotated_pca) <- paste("RotPC", 1:ncol(rotated_pca), sep="_")
obj@meta.data <- cbind(obj@meta.data, rotated_pca[,1:npcs])


#write.csv(loading_mat,paste0('~/Spatial_Rat/', slice, '_VarimaxPca_loadings.csv'), row.names = T)

#png(paste(tag, "varimax_plot.png", sep="_"), width=10, height=7, unit="in", res=300)
print(SpatialFeaturePlot(obj, features=c("RotPC_1", "RotPC_2", "RotPC_3", 
                                         "RotPC_4", "RotPC_5", "RotPC_6")))
#dev.off()

#png(paste(tag, "varimax_plot2.png", sep="_"), width=10, height=7, unit="in", res=300)
print(SpatialFeaturePlot(obj, features=c("RotPC_7", "RotPC_8", "RotPC_9", 
                                         "RotPC_10", "RotPC_11", "RotPC_12")))
#dev.off()

pc_num = 1
RotPC_df = data.frame(genes=rownames(loading_mat), loading=loading_mat[,pc_num])
RotPC_df = RotPC_df[order(RotPC_df$loading, decreasing =F),]
head(RotPC_df, 30)


pericentral_spatial_data = t(GetAssayData(obj)[rownames(obj) %in% pericentral_spatial,])
sum(rownames(pericentral_spatial_data) != rownames(rotated_pca))
cor_mat = cor(cbind(rotated_pca,pericentral_spatial_data))
colnames(cor_mat)
pheatmap(cor_mat[1:ncol(rotated_pca),(ncol(rotated_pca)+1):ncol(cor_mat) ])


######################################################################
############ Calculating the average expression of Hep clusters
######################################################################
## importing the gene expression data
merged_samples_sub = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
merged_samples_all = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC_allfeatures.rds')
dim(merged_samples_all)
dim(merged_samples_sub)

#### checking if the order of cells have been preserved in the data containing all the genes 
sum(colnames(merged_samples_all) != colnames(merged_samples_sub))
head(merged_samples_all)
head(merged_samples_sub)

Resolution = 2.5
resolutions = Resolution
merged_samples_sub <- FindClusters(merged_samples_sub, resolution = Resolution, verbose = FALSE)
table(merged_samples_sub$SCT_snn_res.2.5)

merged_samples_all$SCT_snn_res.2.5 = merged_samples_sub$SCT_snn_res.2.5
merged_samples = merged_samples_all

####################################################
#### defining the mega clusters based on 0.6 resolution
Hep0 = as.character(c(4, 7, 10, 16, 25, 27, 0)) # 27 is the little tail
Hep1 = as.character(c(31, 6, 12, 1, 17, 15)) # 31 is the tail
Hep2 = as.character(c(21, 9, 2, 5)) #
Hep3 = as.character(c(23, 3, 8, 13, 32, 26)) 

######################################################
#### calculating average expression to correlate between clusters and factors in two versions
######################################################
############# Version 1 ############### 
Hep0 = c(4, 7, 10, 16, 25, 27, 0)
Hep1 = c(31, 6, 12, 1, 17, 15, 26)
Hep2 = c(21, 9, 2, 5)
Hep3 = c(23, 3, 8, 13, 32)
##########################################
################ Version 2 ###############
Hep3 = c(23, 3, 8)
Hep2 = c(21, 9, 2, 5)
Hep0 = c(10, 25, 27, 0)
Hep1 = c(31, 6, 12, 1, 17, 15, 26)
HepX = c(13, 32)
HepY = c(4, 7, 16)
###########################
merged_samples$clusters = as.character(merged_samples$SCT_snn_res.2.5)
merged_samples$Hep_clusters = ifelse(merged_samples$clusters %in% as.character(Hep0), 'Hep0', merged_samples$clusters) 
merged_samples$Hep_clusters = ifelse(merged_samples$Hep_clusters %in% as.character(Hep1), 'Hep1', merged_samples$Hep_clusters) 
merged_samples$Hep_clusters = ifelse(merged_samples$Hep_clusters %in% as.character(Hep2), 'Hep2', merged_samples$Hep_clusters) 
merged_samples$Hep_clusters = ifelse(merged_samples$Hep_clusters %in% as.character(Hep3), 'Hep3', merged_samples$Hep_clusters) 

merged_samples$Hep_clusters = ifelse(merged_samples$Hep_clusters %in% as.character(HepX), 'HepX', merged_samples$Hep_clusters) 
merged_samples$Hep_clusters = ifelse(merged_samples$Hep_clusters %in% as.character(HepY), 'HepY', merged_samples$Hep_clusters) 


table(merged_samples$Hep_clusters)
####################################################

################################
##### Generating average expression based heatmap
merged_samples_subset = merged_samples[, merged_samples$Hep_clusters %in% c('Hep0', 'Hep1', 'Hep2', 'Hep3', 'HepX', 'HepY')]
merged_samples_subset = merged_samples[, merged_samples$Hep_clusters %in% c('Hep0', 'Hep1', 'Hep2', 'Hep3')]

dim(merged_samples_subset)

#merged_samples_subset2$manual_annotation = merged_samples_subset$manual_annotation
cluster_names = as.character(merged_samples_subset$SCT_snn_res.2.5 )
cluster_names = as.character(merged_samples_subset$Hep_clusters )
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

top.features <- head(SpatiallyVariableFeatures(obj, selection.method = "moransi"), 100)
cluster_average_exp_df_scaled <- get_scaled_by_gene(cluster_average_exp_df[rownames(cluster_average_exp_df) %in% top.features,]) ## scaling over the clusters of interes

#cluster_average_exp_df_scaled <- get_scaled_by_gene(cluster_average_exp_df) 
head(cluster_average_exp_df_scaled)
dim(cluster_average_exp_df_scaled)
#cluster_average_exp_df_scaled <- cluster_average_exp_df_scaled[,as.character(c(Hep3, Hep2, Hep0, Hep1))] 
#colnames(cluster_average_exp_df_scaled) = paste0('cluster_', colnames(cluster_average_exp_df_scaled))
#head(cluster_average_exp_df_scaled)
#write.csv(cluster_average_exp_df_scaled, file = '~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC_averageCluterExp.csv')
#write.table(cluster_average_exp_df_scaled, file = '~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC_averageCluterExp.txt', row.names = T)

cluster_average_exp_df_scaled$genes=rownames(cluster_average_exp_df_scaled)

########################################################################
################## merging the loading matrix and the cluster averages

loading_mat = loading_mat[,1:npcs]
loading_mat$genes=rownames(loading_mat)
loading_mat$genes=rownames(loading_mat)
head(loading_mat)

loading_mat = loading_mat[,colnames(loading_mat) %in% c(paste0('PC_', 1:4), 'genes')]

merge.df = merge(loading_mat, cluster_average_exp_df_scaled, by.x='genes', by.y='genes')
head(merge.df)
dim(merge.df)
merge.df = merge.df[,-1]

dim(merge.df)
cor_mat=cor(merge.df)
head(cor_mat)

cor_mat2 = cor_mat[1:(ncol(loading_mat)-1),(ncol(loading_mat)):ncol(cor_mat)]
cor_mat2 = cor_mat2[c(1,2),]
rownames(cor_mat2) = c('PC1', 'PC2')

pheatmap(cor_mat2, cluster_cols = F, cluster_rows = F, fontsize=16)
pheatmap(cor_mat2, cluster_cols = T, cluster_rows = F)



## pathway analysis to do
### sample-A: PC5 and PC8 + maybe rotpc18 
### sample-B: pc7, 10, 6, 8 + rotpc 10, 5


file1 = 'MacParland_Catia__Rat_6-1_A1_VIS'
file2 = 'MacParland_Catia__Rat_6-1_B1_VIS'
slice <- file1

slice
input_loading_file = paste0('~/Spatial_Rat/', slice, '_VarimaxPca_loadings.csv')
input_loading_file = paste0('~/Spatial_Rat/', slice, '_pca_loadings.csv')
loading_mat <- read.csv(input_loading_file)
head(loading_mat)


rnk_file_path = paste0('~/Spatial_Rat/pathway_analysis/', slice, '/pca/rank_files/')
rnk_file_path = paste0('~/Spatial_Rat/pathway_analysis/', slice, '/varimax_pca/rank_files/')

rnk_file_path
npc = 12
#### saving markers #####
for(i in 2:(npc+1)){
  
  df = data.frame(gene=loading_mat[,1],
                   score=loading_mat[,i])
  head(df)
  df_ord = df[order(df$score, decreasing = T),]
  head(df_ord)
  
  file_name <- names(loading_mat)[i]
  
  write.table(df_ord, 
              paste0(rnk_file_path, file_name,'.rnk'),
              quote = F, sep = '\t',row.names = FALSE, col.names = FALSE)
}


##### Run GSEA on the DE list



####### Running GSEA on the markers list ########
gsea_jar <- '~/GSEA_Linux_4.1.0/gsea-cli.sh'
dest_gmt_file = '~/Rat_GOBP_AllPathways_no_GO_iea_May_01_2020_symbol.gmt'


#### chose between PCA or varimax loadings
rnk_file_path = paste0('~/Spatial_Rat/pathway_analysis/', slice, '/pca/rank_files/')
rnk_file_path = paste0('~/Spatial_Rat/pathway_analysis/', slice, '/varimax_pca/rank_files/')
rnk_file_path

working_dir = paste0('~/Spatial_Rat/pathway_analysis/', slice, '/pca/gsea_results/')
working_dir = paste0('~/Spatial_Rat/pathway_analysis/', slice, '/varimax_pca/gsea_results/')
working_dir

#### factors to check for the 
### sample-A: PC5 and PC8 + maybe rotpc18 
### sample-B: pc6, 7, 8, 10 + rotpc 10, 5

### zonation signatures:
# sample A: PC1 , PC2 - RotPC2 RotPC8
# sample B: PC1, PC2 - RotPC1, RotPC7

pattern_file =  'PC_2.rnk'
pattern_file = 'RotPC_8.rnk'

list.files(rnk_file_path)
rnk_files <- list.files(rnk_file_path,pattern = pattern_file, full.names = T)

#for(i in 2:length(rnk_files)){
i = 1
rnk_file = rnk_files[i]
print(rnk_file)
  
analysis_name = gsub('.rnk', '',pattern_file)
analysis_name
working_dir = working_dir
  
GSEA_command = paste("",gsea_jar,  "GSEAPreRanked -gmx", dest_gmt_file, "-rnk" ,
                       rnk_file , "-collapse false -nperm 1000 -scoring_scheme weighted -rpt_label ",
                       analysis_name,"  -plot_top_x 20 -rnd_seed 12345  -set_max 200 -set_min 15 -zip_report false -out" , #set min=15
                       working_dir, " > gsea_output.txt")
system(GSEA_command)
#}





