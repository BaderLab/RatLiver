library(scmap)
library(celldex)
library(plyr)
library(stats)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(scales)
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
#merged_samples = readRDS('~/rat_sham_sn_data/DropletQC_results/sham_sn_merged_data_dropletQC.rds')
merged_samples = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_standardQC.rds')
merged_samples <- readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
gene_names = c('Ptprc', 'Calcrl', 'Nkg7', 'Cd3e', 'Marco', 'Lyz2', 'Cd19', 'Ms4a1', 'Stab2')
gene_names = c('Ptprc', 'Sox9', 'Acta2', 'Cyp2e1', 'Pck1', 'Cyp1a2')

merged_samples <- FindNeighbors(merged_samples,reduction="harmony",verbose=T)

Resolution = 0.6
merged_samples <- FindClusters(merged_samples, resolution = Resolution, verbose = FALSE)


i = 7
gene_name = gene_names[i]
df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap_h')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap_h')[,2], 
                      library_size= merged_samples$nCount_RNA, 
                      mito_perc=merged_samples$mito_perc, 
                      n_expressed=merged_samples$nFeature_RNA,
                      cluster=merged_samples$SCT_snn_res.0.6, 
                      cell_status = merged_samples$cell_status,
                      nuclear_fraction=merged_samples$nuclear_fraction, 
                      Alb=GetAssayData(merged_samples)['Alb',], 
                      sample_name = merged_samples$sample_name, 
                      gene=GetAssayData(merged_samples)[gene_name,],
                      strain = merged_samples$strain, 
                      umi=colnames(merged_samples))

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=gene))+geom_point(alpha=0.6,size=1.4)+
  scale_color_viridis('Expression\nValue', direction = -1)+theme_classic()+
  theme(text = element_text(size=16.5),
        plot.title = element_text(hjust = 0.5),
        legend.title=element_text(size = 10))+
  ggtitle(paste0(gene_name)) 



ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cluster))+geom_point(size=1,alpha=0.6)+theme_classic()

##################################################################
############## calculating the cluster average expression of the new data ############## 

nb.cols <- length(names(table(merged_samples$seurat_clusters)))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols) #Pastel1

merged_samples <- FindNeighbors(merged_samples, reduction = "harmony", dims = 1:30)
merged_samples <- FindClusters(merged_samples, resolution = 0.7)
df_umap$cluster = merged_samples$SCT_snn_res.0.7

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cluster))+geom_point(size=1,alpha=0.6)+
  theme_classic()+scale_color_manual(values = c(colorPalatte)) #mycolors
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
###### importing cluster average means of reference maps 
ref_cluster_average_exp_df <- readRDS('~/RatLiver/Results/rat_old_cluster_average_exp_all.rds')
ref_cluster_average_exp_df <- readRDS('~/RatLiver/Results/rat_new_cluster_average_exp_all.rds')

### refine the annotation names for each of the clusters
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

dev.off()
gridExtra::grid.table(annotation.df[,-3])
dev.off()




#annotation.df_set1 = annotation.df
#annotation.df_set2 = annotation.df

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

#annotation.df_mm = annotation.df
############### Adding the annotation results to the umap and seurat data object ###############
## use each of the annotation.df individually and add their meta data to the df_umap
#df_umap_backup = df_umap
#annotation.df = annotation.df_mm
df_umap$cluster = as.character(df_umap$cluster)
annotation.df$clusters_num = as.character(0:(nrow(annotation.df)-1))

df_umap = merge(df_umap, annotation.df[,2:4], by.x='cluster', by.y='clusters_num', all.x=T, order=F)
#### re-ordering the rows based on UMIs of the initial data
df_umap <- df_umap[match(colnames(merged_samples),df_umap$umi),]
df_umap$umi == colnames(merged_samples)

head(df_umap)
colnames(df_umap)

### saving the umap containing the annotation based on all maps
colnames(df_umap)[(ncol(df_umap)-5):ncol(df_umap)] = c('annot_TLH', 'annot_TLH_g', 'annot_IM', 'annot_IM_g', 'annot_mm', 'annot_mm_g')
#saveRDS(df_umap, '~/rat_sham_sn_data/standardQC_results/merged_df_umap_annot.rds')


df_umap$label = paste0(df_umap$annotation_g, '(', df_umap$cluster, ')')
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=label))+geom_point(size=1,alpha=0.5)+
  theme_classic()+scale_color_manual(values = c(colorPalatte, 'black')) #mycolors
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=label))+geom_point(size=1.6,alpha=0.8)+
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




#######################################################
#### adding annotation info to the seurat object######
merged_metadata = cbind(merged_samples@meta.data, df_umap[,c(13, 15:ncol(df_umap))])
head(merged_metadata)
merged_samples@meta.data = merged_metadata
head(merged_samples)
#saveRDS(merged_samples, '~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
merged_samples <- readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')


