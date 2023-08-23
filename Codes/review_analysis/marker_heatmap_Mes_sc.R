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

########################################################################################
################ Rat - single cell RNAseq - TLH map -  mesenchymal clusters ################
########################################################################################

rat_cluster_average.df <- readRDS('Results/rat_old_cluster_average_exp_all.rds')
mes_sub_seur = readRDS('Results/old_samples/Mesenchymal_subclusters.rds')
sc_HVGs = VariableFeatures(FindVariableFeatures(mes_sub_seur))

umap_df = getUmapDF(mes_sub_seur)
umap_df$gcluster = mes_sub_seur$final_cluster
umap_df$res.0.2 = mes_sub_seur$SCT_snn_res.0.2
umap_df$res.0.4 = mes_sub_seur$SCT_snn_res.0.4

ggplot(umap_df,aes(UMAP_1, UMAP_2, color=gcluster))+geom_point()+theme_bw()

ggplot(umap_df,aes(UMAP_1, UMAP_2, color=res.0.4))+
  geom_point(alpha=0.6,size=1.8)+theme_classic()+
  scale_color_brewer(name='subcluster',palette = 'Set2')+
  theme(text = element_text(size=13.5),
        plot.title = element_text(hjust = 0.5),
        legend.title=element_text(size = 10))+
  ggtitle('scRNA-seq TLH map mesenchymal subclusters') 


mes_sub_seur$mes_subclust = paste0('rat_mes_', mes_sub_seur$SCT_snn_res.0.4)

cluster_names_types = names(table(mes_sub_seur$mes_subclust ))
cluster_names = mes_sub_seur$mes_subclust

### dividing the expression matrix based on the clusters
cluster_expression <- sapply(1:length(cluster_names_types), function(i){
  a_cluster_name = cluster_names_types[i]
  GetAssayData(mes_sub_seur, 'data')[,cluster_names == a_cluster_name] 
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


cluster_average_exp_df <- cluster_average_exp_df[,!colnames(cluster_average_exp_df) %in% c("rat_mes_2","rat_mes_3")]
## scale and center all the genes in the matrix
cluster_average_exp_df <- get_scaled_by_gene(cluster_average_exp_df) ## scaling over the clusters of interest
head(cluster_average_exp_df)
cluster_average_exp_df$rat_ID = rownames(cluster_average_exp_df)
rat_cluster_average.df = cluster_average_exp_df
head(rat_cluster_average.df)
colnames(rat_cluster_average.df)[1:4] = paste0('Mes_',c(0,1,4,5))
head(rat_cluster_average.df)


##### single cell TLH map highly variables genes
sc_mes_subclusters_average = rat_cluster_average.df
sc_HVGs


########## selecting the number of related genes
FB_list = c("Dpt", "Entpd2", "Gsn", "Col1a1", "Flt1")
HSC_list = c("Reln", "Pth1", "Lrat", "Hgf", "Pth1r", "Ecm1")
VSMC_list = c("Acta2", "Tagln", "Tpm2")

markers_list = c(FB_list, HSC_list, VSMC_list)
dim(rat_cluster_average.df)


#### ordering the genes based on the markers list
markers_list[!markers_list%in% rownames(rat_cluster_average.df)]
markers_list = markers_list[markers_list%in% rownames(rat_cluster_average.df)]
rat_cluster_average_subset = rat_cluster_average.df[markers_list, ]
rat_cluster_average_subset = rat_cluster_average_subset[,colnames(rat_cluster_average_subset)!='rat_ID']

anno_df_genes = c(rep(NA, nrow(rat_cluster_average_subset)))
sapply(1: nrow(rat_cluster_average_subset) , function(i){
  if (rownames(rat_cluster_average_subset)[i] %in% FB_list) anno_df_genes[i] <<- 'FB'
  if (rownames(rat_cluster_average_subset)[i] %in% HSC_list) anno_df_genes[i] <<- 'HSC'
  if (rownames(rat_cluster_average_subset)[i] %in% VSMC_list) anno_df_genes[i] <<- 'VSMC'
})
anno_df_genes
anno_df_genes = data.frame(Annotation=anno_df_genes)
rownames(anno_df_genes) = markers_list


head(rat_cluster_average_subset)
rownames(rat_cluster_average_subset)
head(anno_df_genes)
colnames(anno_df_genes) = 'Group'
pheatmap(rat_cluster_average_subset, cluster_rows = FALSE, color= inferno(20, direction = +1),  
         cluster_cols = FALSE, annotation_row = anno_df_genes)





########################################################################################
################ Rat - single nuc RNAseq - mesenchymal clusters ################
########################################################################################


mes_data <- readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_mesenchymal_subclusters_merged_data_res2.5_cluster24Only.rds')
mes_data$mes_subclust = paste0('Mes_', mes_data$SCT_snn_res.1)
mes_data$res.1 = paste0(mes_data$SCT_snn_res.1)



umap_df = getUmapDF(mes_data)
umap_df$res.1 = mes_data$res.1

ggplot(umap_df,aes(UMAP_1, UMAP_2, color=res.1))+
  geom_point(alpha=0.6,size=2.2)+theme_classic()+
  scale_color_brewer(name='subcluster',palette = 'Set1')+
  theme(text = element_text(size=13.5),
        plot.title = element_text(hjust = 0.5),
        legend.title=element_text(size = 10))+
  ggtitle('snRNA-seq map mesenchymal subclusters') 




sn_HVGs = VariableFeatures(FindVariableFeatures(mes_data))

cluster_names_types = names(table(mes_data$mes_subclust ))
cluster_names = mes_data$mes_subclust

### dividing the expression matrix based on the clusters
cluster_expression <- sapply(1:length(cluster_names_types), function(i){
  a_cluster_name = cluster_names_types[i]
  GetAssayData(mes_data, 'data')[,cluster_names == a_cluster_name] 
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
head(cluster_average_exp_df)
cluster_average_exp_df$rat_ID = rownames(cluster_average_exp_df)
rat_cluster_average.df = cluster_average_exp_df
head(rat_cluster_average.df)

sn_mes_subclusters_average = rat_cluster_average.df
sn_HVGs

########## selecting the number of related genes

FB_list =  c("Dpt", "Serpinf1", "Mgst1", "Clu", "Flt1")
HSC_list = c("Reln", "Pth1", "Lrat", "Hgf", "Pth1r")
VSMC_list = c("Acta2", "Tagln")
markers_list = c(FB_list, HSC_list, VSMC_list)


markers_list = c(FB_list, HSC_list, VSMC_list)
dim(rat_cluster_average.df)


#### ordering the genes based on the markers list
markers_list[!markers_list%in% rownames(rat_cluster_average.df)]
markers_list = markers_list[markers_list%in% rownames(rat_cluster_average.df)]
rat_cluster_average_subset = rat_cluster_average.df[markers_list, ]
rat_cluster_average_subset = rat_cluster_average_subset[,colnames(rat_cluster_average_subset)!='rat_ID']

anno_df_genes = c(rep(NA, nrow(rat_cluster_average_subset)))
sapply(1: nrow(rat_cluster_average_subset) , function(i){
  if (rownames(rat_cluster_average_subset)[i] %in% FB_list) anno_df_genes[i] <<- 'FB'
  if (rownames(rat_cluster_average_subset)[i] %in% HSC_list) anno_df_genes[i] <<- 'HSC'
  if (rownames(rat_cluster_average_subset)[i] %in% VSMC_list) anno_df_genes[i] <<- 'VSMC'
})
anno_df_genes
anno_df_genes = data.frame(Annotation=anno_df_genes)
rownames(anno_df_genes) = markers_list


#rat_cluster_average_subset <- rat_cluster_average_subset[,paste0('Mes_',c(5, 4, 1, 0))]

head(rat_cluster_average_subset)
rownames(rat_cluster_average_subset)
colnames(anno_df_genes) = 'Group'
rat_cluster_average_subset = rat_cluster_average_subset[,c('Mes_1', 'Mes_0', 'Mes_2')]
pheatmap(rat_cluster_average_subset, cluster_rows = FALSE, color= inferno(20, direction = +1),  
         cluster_cols = FALSE, annotation_row = anno_df_genes)



################################################################################

head(sn_mes_subclusters_average)
head(sc_mes_subclusters_average)

sn_mes_subclusters_average2 = sn_mes_subclusters_average
sc_mes_subclusters_average2 = sc_mes_subclusters_average

colnames(sn_mes_subclusters_average2)[1:3] = paste0(colnames(sn_mes_subclusters_average2)[1:3], '.sn')
colnames(sc_mes_subclusters_average2)[1:4] = paste0(colnames(sc_mes_subclusters_average2)[1:4], '.sc')

head(sc_mes_subclusters_average2)
head(sn_mes_subclusters_average2)

merged_df = merge(sn_mes_subclusters_average2, sc_mes_subclusters_average2, by.x='rat_ID', by.y='rat_ID')
head(merged_df)

dim(merged_df)
dim(merged_df[c(sc_HVGs,sn_HVGs),])
merged_df = merged_df[merged_df$rat_ID %in% c(sc_HVGs,sn_HVGs),]

merged_df = merged_df[merged_df$rat_ID %in% sc_HVGs[sc_HVGs%in%sn_HVGs],]
cor_mat = cor(merged_df[,-1])
pheatmap(cor_mat[1:(ncol(sn_mes_subclusters_average2)-1), (ncol(sn_mes_subclusters_average2)):ncol(cor_mat)])
sc_HVGs











