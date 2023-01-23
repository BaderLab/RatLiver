library(scmap)
library(celldex)
library(plyr)
library(stats)
library(ggpubr)

library(RColorBrewer)
library(viridis)
library(scales)

source('~/RatLiver/Codes/Functions.R')
Initialize()

get_scaled_by_gene <- function(x){
  y <- scale(t(x), center = TRUE, scale = TRUE)
  y[is.nan(y)] <- 0
  y2 <- t(y)
  y2 <- as.data.frame(y2)
  return(y2)
}

########################################################
############ Importing the new data to be annotated #########
merged_samples = readRDS('~/rat_sham_sn_data/sham_singleNuc_merged.rds')
### converting the seurat object to singleCellExperiment
merged_samples.sce = as.SingleCellExperiment(merged_samples)
rowData(merged_samples.sce)$feature_symbol = rownames(merged_samples.sce)


############ generating the reference data ############
new_data_scCLustViz_object <- "~/RatLiver/Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_labelCor.RData"
old_data_scClustViz_object <- "~/RatLiver/Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"

#########
load(new_data_scCLustViz_object)
load(old_data_scClustViz_object)

ref_data = your_scRNAseq_data_object
ref_data$sample_name = sapply(str_split(colnames(ref_data), '_'), '[[', 2)
ref_data$sample_name = ifelse(ref_data$orig.ident=='rat_DA_01_reseq', 'DA-1', 
                              ifelse(ref_data$orig.ident=='rat_DA_M_10WK_003', 'DA-2',
                                     ifelse(ref_data$orig.ident=='rat_Lew_01', 'Lew-1', 'Lew-2')))

df_umap <- data.frame(UMAP_1=getEmb(ref_data, 'umap')[,1], 
                      UMAP_2=getEmb(ref_data, 'umap')[,2], 
                      clusters=ref_data$cluster,
                      sample=ref_data$sample_name,
                      umi=colnames(ref_data))

##### adding the final annotations to the dataframe
set1_info <- read.csv('figure_panel/set-1-final-info.csv')
set1_info <- read.csv('~/RatLiver/figure_panel/set-2-final-info-updated.csv')
colnames(set1_info)[1] = 'clusters'
#set1_info <- set1_info[1:18,] # 14
set1_info$clusters = as.character(set1_info$clusters)
df_umap = merge(df_umap, set1_info[,1:3], by.x='clusters', by.y='clusters', all.x=T, order=F)
#### re-ordering the rows based on UMIs of the initial data
df_umap <- df_umap[match(colnames(ref_data),df_umap$umi),]
df_umap$umi == colnames(ref_data)

ref_data.sce <- as.SingleCellExperiment(ref_data)
df_umap$label = gsub(pattern = 'Inflammatory',replacement = 'Inf', df_umap$label)
ref_data.sce$cell_type = df_umap$label
table(ref_data.sce$cell_type)


# Create scmap-cluster reference by first selecting the most informative features
rowData(ref_data.sce)$feature_symbol = rownames(ref_data.sce)
ref_data.sce <- scmap::selectFeatures(ref_data.sce, suppress_plot=FALSE)
# Inspect the first 50 genes selected by scmap
rownames(ref_data.sce)[which(rowData(ref_data.sce)$scmap_features)][1:200]

colData(ref_data.sce)$cell_type1 =  ref_data.sce$cell_type
ref_data.sce <- indexCluster(ref_data.sce)
head(metadata(ref_data.sce)$scmap_cluster_index)
heatmap(as.matrix(metadata(ref_data.sce)$scmap_cluster_index))


##################################################################
######### Projecting the new data on the reference map ##########
scmapCluster_results <- scmapCluster(
  projection = merged_samples.sce, 
  index_list = list(
    yan = metadata(ref_data.sce)$scmap_cluster_index))

plot(getSankey(
  colData(ref_data.sce)$cell_type1, 
  scmapCluster_results$scmap_cluster_labs[,'yan'],
  plot_height = 400))


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
                      umi=colnames(merged_samples))



df_umap$scmap_cluster = scmapCluster_results$scmap_cluster_labs
head(df_umap)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=scmap_cluster))+geom_point(size=1,alpha=0.6)+theme_bw()
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cluster))+geom_point(size=1,alpha=0.6)+theme_bw()

##################################################################





##################################################################
############## calculating the cluster average expression of the new data ############## 

nb.cols <- length(names(table(merged_samples$seurat_clusters)))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols) #Pastel1

merged_samples <- FindNeighbors(merged_samples, reduction = "harmony", dims = 1:30)
merged_samples <- FindClusters(merged_samples, resolution = 0.4)
df_umap$cluster = merged_samples$SCT_snn_res.0.4

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cluster))+geom_point(size=1,alpha=0.6)+
  theme_bw()+scale_color_manual(values = c(mycolors)) #colorPalatte
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cluster))+geom_point(size=1,alpha=0.6)+theme_classic()

merged_samples$cluster = df_umap$cluster
cluster_names_types = names(table(df_umap$cluster))

### dividing the expression matrix based on the clusters
cluster_expression <- sapply(1:length(cluster_names_types), function(i){
  a_cluster_name = cluster_names_types[i]
  GetAssayData(merged_samples, 'data')[,merged_samples$cluster == a_cluster_name] 
}, simplify = F)

names(cluster_expression) = cluster_names_types
lapply(cluster_expression, head)

## calculate the average expression of each gene in each cluster
cluster_average_exp <- lapply(cluster_expression, function(x){
  df = data.frame(average=rowSums(x)/ncol(x))
  return(df)
})
lapply(cluster_average_exp, dim)
lapply(cluster_average_exp, head)
## Concatenate all the clusters together to make a matrix
cluster_average_exp_df = do.call(cbind,cluster_average_exp)
#colnames(cluster_average_exp_df) = paste0('cluster_',names(cluster_average_exp))
colnames(cluster_average_exp_df) = paste0('c_',names(cluster_average_exp))
head(cluster_average_exp_df)

## scale and center all the genes in the matrix
cluster_average_exp_df <- get_scaled_by_gene(cluster_average_exp_df) ## scaling over the clusters of interest
cluster_average_exp_df$rat_ID=rownames(cluster_average_exp_df)
head(cluster_average_exp_df)

cluster_average_exp_df = cluster_average_exp_df[rownames(cluster_average_exp_df)%in%VariableFeatures(merged_samples),]
head(cluster_average_exp_df)
dim(cluster_average_exp_df)


##################################################################

ref_cluster_average_exp_df <- readRDS('~/RatLiver/Results/rat_old_cluster_average_exp_all.rds')
ref_cluster_average_exp_df <- readRDS('~/RatLiver/Results/rat_new_cluster_average_exp_all.rds')

colnames_set1 = c("Hep (0)", "Hep (1)", "Hep (12)", "Hep (15)", "Hep (16)", "Hep (2)", "Hep (4)", 
  "Hep (6)", "Hep (8)", "Lyz2/Cd74 Mo/Mac (9)", "Endothelial (11)", "Endothelial (3)", "Lymphocyte (13)",  
  "Marco/Cd5l Mac (10)",  "Marco/Cd5l Mac (5)", "Mesenchymal (14)", "Mesenchymal (7)", "rat_ID" )

colnames_set2 = c("pDC (17)", "Cd3+ (10)", "Erythroid (5)", "Hep (0)" , "Hep (1)", "Hep (14)",
                  "Hep (15)", "Hep (2)", "Hep (3)", "Hep (9)", "Mo/Mac/cDC (11)", "Endo (4)",
                  "Endo (6)", "Marco/Cd5l Mac (8)", "B cell (12)", "gd T cell (7)", 
                  "Marco/Cd5l Mac (13)", "Mesenchymal (16)", "rat_ID" )

colnames(ref_cluster_average_exp_df) = colnames_set1
colnames(ref_cluster_average_exp_df) = colnames_set2
nFeatures = 2000

old_data_scClustViz_object <- "~/RatLiver/Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
new_data_scCLustViz_object <- "~/RatLiver/Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_labelCor.RData"

load(old_data_scClustViz_object)
load(new_data_scCLustViz_object)

ref_data = your_scRNAseq_data_object
DefaultAssay(ref_data) <- 'RNA'
ref_data <- FindVariableFeatures(ref_data, nfeatures=nFeatures)
HVGs <- VariableFeatures(ref_data)

ref_cluster_average_exp_df = ref_cluster_average_exp_df[rownames(ref_cluster_average_exp_df)%in%HVGs,]


########################################################################
######## merging the Healthy Rat ref data with the new samples ######## 
average_exp_merged = merge(cluster_average_exp_df, ref_cluster_average_exp_df, by.x='rat_ID', by.y='rat_ID')

dim(average_exp_merged) ## additional layers to ne to be added to the filter -> 10126 one-2-one orthos
head(average_exp_merged)
rat_cor_mat = cor(average_exp_merged[,-1],method = 'pearson')
colnames(rat_cor_mat)

colnames(rat_cor_mat) <- gsub('Inflammatory', 'Inf', colnames(rat_cor_mat))
rownames(rat_cor_mat) = colnames(rat_cor_mat)
number_of_clusters = length(cluster_names_types)
rat_cor_mat = rat_cor_mat[1:number_of_clusters,(number_of_clusters+1):ncol(rat_cor_mat)] 
pheatmap(rat_cor_mat,color = inferno(20),  clustering_method='ward.D2')


get_matched_label <- function(index, cor.mat, thr){
  label = 'unknown'
  a_clust_value = cor.mat[index,]
  tmp_label =  names(which.max(cor.mat[index,]))
  if(a_clust_value[tmp_label]>thr) label = tmp_label
  return(label)
}

annotation.df = data.frame(cluster=rownames(rat_cor_mat),  
                           annotation=sapply(1:nrow(rat_cor_mat), 
                                             function(i) get_matched_label(index=i, rat_cor_mat, thr=0.3) ))

annotation.df$annotation_g = sub(" \\(.*", "", annotation.df$annotation)
#saveRDS(annotation.df, 'Objects/Integrated_Syngeneic_res0.5_set1_based_annotation.rds')
#saveRDS(annotation.df, 'Objects/Integrated_Rejection_res0.3_set2_based_annotation.rds')
#saveRDS(annotation.df, 'Objects/Integrated_Tolerance_res0.2_set2_based_annotation.rds')
dev.off()
gridExtra::grid.table(annotation.df[,-3])
dev.off()


##################################################################################
###########  Annotation based on similarity with the mouse samples ###############

mouse_cluster_average_df <- readRDS('~/RatLiver/Results/mouse_cluster_average_exp_HVGs.rds')
colnames_to_change = colnames(mouse_cluster_average_df)[1:(ncol(mouse_cluster_average_df)-1)] 
colnames(mouse_cluster_average_df)[1:(ncol(mouse_cluster_average_df)-1)] = substr(colnames_to_change,4,nchar(colnames(mouse_cluster_average_df)))
### converting the IDs
rat_to_mouse_genes <- readRDS('~/XSpecies/rat_to_mouse_genes.rds')
rat_to_mouse_genes = rat_to_mouse_genes[rat_to_mouse_genes$mmusculus_homolog_orthology_type=='ortholog_one2one',]
rat_to_mouse_genes <- rat_to_mouse_genes[,c('mmusculus_homolog_associated_gene_name', 'symbol')]
dim(rat_to_mouse_genes)
head(rat_to_mouse_genes)

#### adding the rat ortholog gene symbols
mouse_cluster_average.df.homolog <- merge(mouse_cluster_average_df, rat_to_mouse_genes, by.x='mouse_genes',
                                          by.y='mmusculus_homolog_associated_gene_name', sort=F)
head(mouse_cluster_average.df.homolog)
head(cluster_average_exp_df) 

merged_rat_mouse = merge(cluster_average_exp_df, mouse_cluster_average.df.homolog, by.y='symbol', by.x='rat_ID')
head(merged_rat_mouse)
dim(merged_rat_mouse)

number_of_clusters = length(cluster_names_types)

cor.mat <- cor(merged_rat_mouse[,!colnames(merged_rat_mouse) %in% c('mouse_genes', 'rat_ID')],method = 'pearson')
rownames(cor.mat) <- gsub('Inflammatory', 'Inf', rownames(cor.mat))
cor.mat = cor.mat[1:number_of_clusters,(number_of_clusters+1):ncol(cor.mat)]
pheatmap(cor.mat,color = inferno(20),main='', clustering_method = 'ward.D2') 

annotation.df = data.frame(cluster=rownames(cor.mat),  
                           annotation=sapply(1:nrow(cor.mat), 
                                             function(i) get_matched_label(index=i, cor.mat, thr=0.3) ))

annotation.df$annotation_g = sub(" \\(.*", "", annotation.df$annotation)
annotation.df = annotation.df[match(paste0('c_', 0:(number_of_clusters-1)),annotation.df$cluster),]
gridExtra::grid.table(annotation.df[,-3])
dev.off()


############### Adding the annotation results to the umap and seurat data object ###############

df_umap$cluster = as.character(df_umap$cluster)
annotation.df$clusters_num = as.character(0:(nrow(annotation.df)-1))

df_umap = merge(df_umap, annotation.df[,2:4], by.x='cluster', by.y='clusters_num', all.x=T, order=F)
#### re-ordering the rows based on UMIs of the initial data
df_umap <- df_umap[match(colnames(merged_samples),df_umap$umi),]
df_umap$umi == colnames(merged_samples)

head(df_umap)

### saving the umap containing the annotation based on all maps
#colnames(df_umap)[15:ncol(df_umap)] = c('annot_mm', 'annot_mm_g', 'annot_TLH', 'annot_TLH_g', 'annot_IM', 'annot_IM_g')
#saveRDS(df_umap, '~/rat_sham_sn_data/merged_df_umap_annot.rds')


#saveRDS(df_umap, 'Objects/Integrated_Syngeneic_mt15_lib1500_df_UMAP.rds')
df_umap$label = paste0(df_umap$annotation_g, '(', df_umap$cluster, ')')
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=label))+geom_point(size=1,alpha=0.5)+
  theme_classic()+scale_color_manual(values = c(colorPalatte, 'black')) #mycolors
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=label))+geom_point(size=2,alpha=1)+
  theme_classic()+scale_color_manual(values = c(colorPalatte, 'black')) #mycolors

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=annot_mm))+geom_point(size=1.5,alpha=0.7)+
  theme_classic()+scale_color_manual(values = c(colorPalatte, 'black')) #mycolors
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=annot_TLH))+geom_point(size=1.5,alpha=0.7)+
  theme_classic()+scale_color_manual(values = c(colorPalatte, 'black')) #mycolors
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=annot_IM))+geom_point(size=1.5,alpha=0.7)+
  theme_classic()+scale_color_manual(values = c(colorPalatte, 'black')) #mycolors

df_umap$label_set2 = df_umap$label  
###### generating count barplots

df <- data.frame(sample_type = merged_samples$sample_name, 
                 cluster = as.character(df_umap$cluster))
rownames(df) = NULL
counts <- ddply(df, .(df$sample_type, df$cluster), nrow)
names(counts) <- c("sample_type", "cluster", "Freq")

###### ordering the bars in decreasing order
freq.df = data.frame(table(df$cluster))
freq.df = freq.df[order(freq.df$Freq, decreasing = T),]
cluster_orders = as.character(freq.df$Var1)
counts$cluster= factor(counts$cluster, levels = cluster_orders ) 


ggplot(data=counts, aes(x=cluster, y=Freq, fill=sample_type)) +
  geom_bar(stat="identity",color='black')+theme_classic()+#+scale_fill_brewer(palette = "Blues")+
  ylab('Counts')+xlab('Clusters')+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=10,angle=90,color='black'),
        legend.title = element_blank())+
  xlab('')

counts_split <- split( counts , f = counts$cluster )
counts_split_norm <- lapply(counts_split, function(x) {x$Freq=x$Freq/sum(x$Freq);x})
counts_norm <- do.call(rbind,counts_split_norm )

ggplot(data=counts_norm, aes(x=cluster, y=Freq, fill=sample_type)) +
  geom_bar(stat="identity",color='black',alpha=0.9)+theme_classic()+
  ylab('Fraction of sample per cell type (%)')+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=12.5,angle=90,color='black'),
        legend.title = element_blank()) +  
  xlab('')






