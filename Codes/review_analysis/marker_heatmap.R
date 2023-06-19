source('Codes/Functions.R')
source('~/XSpecies/Codes/convert_human_to_ortholog_functions.R')
library(Hmisc)
library(RColorBrewer)
library(viridis)
library(reshape2)
########################################################
############ Importing the nuc-seq data  #########

## importing the gene expression data
merged_samples_sub = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
merged_samples_all = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC_allfeatures.rds')


'Cd68'%in% rownames(GetAssayData(merged_samples_all, 'scale.data'))
'Cd68'%in% rownames(GetAssayData(merged_samples_sub, 'scale.data'))

markers_test <- c('Alb', 'Apoa1', 'Apoc1', 'Apoc3', 'Apoe', 'Fabp1', 'Itih4', 'Orm1', 'Pigr', 'Serpina1', 'Tf', 'Ttr')

markers_test %in% rownames(GetAssayData(merged_samples_all, 'scale.data'))
markers_test[!markers_test %in% rownames(GetAssayData(merged_samples_all, 'scale.data'))]

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

CV_Hep = c(23, 8, 3, 13, 32, 26, 12, 17, 1, 15) 
CV_Hep_Markers = c('Akr1c1', 'Ahr', 'Glul', 
                   'Notum', 'Rcan1', 'Cyp2f4', 'Cyp27a1', 'Cyp7a1', 'Cyp2e1', 'Cyp8b1')

other_Hep = c(21, 6, 20, 31, 9)
other_Hep_Markers = c('Mt2A', 'Hamp')

Periportal_Hep = c(27, 0, 4, 10, 16, 7, 2, 5) #removed 15
Periportal_Hep_Markers=c('Alb', 'Apoa1', 'Apoc1', 'Apoc3', 'Apoe', 'Fabp1', 
                         'Itih4', 'Orm1', 'Pigr', 'Serpina1', 'Tf', 'Ttr') 


Cholagiocytes = 29
Cholagiocytes_Markers = c('Epcam', 'Sox9', 'Anxa4')


General_Macrophage_Markers = c('Cd68', 'Clec4f')

Non_inf_mac = 19 
Non_inf_mac_Markers = c('Cd5l', 'Marco', 'Cd163', 'C1qa', 'C1qc', 'Aif1', 
                        'Slc11a', 'Hmox1', 'Vsig4', 'Ptas1', 'Ccl6')

Inf_mac = 33 
Inf_mac_Markers = c('Cd74', 'RT1-Ba', 'RT1-Bb', 'RT1-Da', 'RT1-Db1', 'Lyz2')

Mesenchymal = c(24) 
Mesenchymal_Markers = c('Col3a1', 'Colec10', 'Colec11', 'Ecm1') #'Stab2', 'Vwf', 


Endothelial = c(11, 30) 
Endothelial_Markers = c('Aqp1', 'Fcgr2b', 'Gpr182', 'Lyve1', 'Stab2')


clusters_to_include = c(CV_Hep, other_Hep, Periportal_Hep, Cholagiocytes, 
                        Non_inf_mac, Inf_mac, Mesenchymal, Endothelial)

markers_list = c(CV_Hep_Markers, other_Hep_Markers, Periportal_Hep_Markers, Cholagiocytes_Markers, 
                 General_Macrophage_Markers, Non_inf_mac_Markers, Inf_mac_Markers, 
                 Mesenchymal_Markers, Endothelial_Markers) 


merged_samples@active.ident <- merged_samples$SCT_snn_res.2.5 
merged_samples_subset = merged_samples[, merged_samples$SCT_snn_res.2.5 %in% clusters_to_include]
dim(merged_samples_subset)
dim(merged_samples)

annot = c(rep(NA, ncol(merged_samples_subset)))
sapply(1:ncol(merged_samples_subset) , function(i){
  if (merged_samples_subset$SCT_snn_res.2.5[i] %in% CV_Hep) annot[i] <<- 'Hep 1'
  if (merged_samples_subset$SCT_snn_res.2.5[i] %in% other_Hep) annot[i] <<- 'Hep'
  if (merged_samples_subset$SCT_snn_res.2.5[i] %in% Periportal_Hep) annot[i] <<- 'Hep 3'
  if (merged_samples_subset$SCT_snn_res.2.5[i] %in% Cholagiocytes) annot[i] <<- 'cholagiocyte'
  if (merged_samples_subset$SCT_snn_res.2.5[i] %in% Non_inf_mac) annot[i] <<- 'nonInfMac'
  if (merged_samples_subset$SCT_snn_res.2.5[i] %in% Inf_mac) annot[i] <<- 'InfMac'
  if (merged_samples_subset$SCT_snn_res.2.5[i] %in% Mesenchymal) annot[i] <<- 'mesenchymal'
  if (merged_samples_subset$SCT_snn_res.2.5[i] %in% Endothelial) annot[i] <<- 'endothelial'
  if(i %% 100 == 0 ) print(i)
})

merged_samples_subset$manual_annotation = annot



merged_samples_subset@active.ident <- factor(x = merged_samples_subset@active.ident, levels = unique(clusters_to_include))
merged_samples_subset@active.ident <- factor(x = merged_samples_subset$manual_annotation, levels = unique(merged_samples_subset$manual_annotation))
#saveRDS(merged_samples_subset, '~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC_allfeatures_heatmap_v.rds')
#merged_samples_subset <- readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC_allfeatures_heatmap_v.rds')


DoHeatmap(
  merged_samples_subset,
  features = markers_list,
  cells = NULL,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data", #"scale.data", data, counts
  assay = 'RNA',
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)

c('Ptas1', 'Slc11a', 'Cd68', 'Anxa4', 'Sox9', 'Cyp7a1') %in% rownames(merged_samples_subset)
counts(merged_samples_sub)
dim(merged_samples_sub)





################################
##### Generating average expression based heatmap

get_scaled_by_gene <- function(x){
  y <- scale(t(x), center = TRUE, scale = TRUE)
  y[is.nan(y)] <- 0
  y2 <- t(y)
  y2 <- as.data.frame(y2)
  return(y2)
}

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

anno_df = c(rep(NA, ncol(cluster_average_exp_df_scaled)))
sapply(1:ncol(merged_samples_subset) , function(i){
  if (colnames(cluster_average_exp_df_scaled)[i] %in% CV_Hep) anno_df[i] <<- 'Hep 1'
  if (colnames(cluster_average_exp_df_scaled)[i] %in% other_Hep) anno_df[i] <<- 'Hep 2'
  if (colnames(cluster_average_exp_df_scaled)[i] %in% Periportal_Hep) anno_df[i] <<- 'Hep 3'
  if (colnames(cluster_average_exp_df_scaled)[i] %in% Cholagiocytes) anno_df[i] <<- 'Cholagiocyte'
  if (colnames(cluster_average_exp_df_scaled)[i] %in% Non_inf_mac) anno_df[i] <<- 'Non-inf Mac'
  if (colnames(cluster_average_exp_df_scaled)[i] %in% Inf_mac) anno_df[i] <<- 'Inf Mac'
  if (colnames(cluster_average_exp_df_scaled)[i] %in% Mesenchymal) anno_df[i] <<- 'Mesenchymal'
  if (colnames(cluster_average_exp_df_scaled)[i] %in% Endothelial) anno_df[i] <<- 'Endothelial'
})


anno_df_genes = c(rep(NA, nrow(cluster_average_exp_df_scaled)))
sapply(1:ncol(merged_samples_subset) , function(i){
  if (rownames(cluster_average_exp_df_scaled)[i] %in% CV_Hep_Markers) anno_df_genes[i] <<- 'Hep 1'
  if (rownames(cluster_average_exp_df_scaled)[i] %in% other_Hep_Markers) anno_df_genes[i] <<- 'Hep 2'
  if (rownames(cluster_average_exp_df_scaled)[i] %in% Periportal_Hep_Markers) anno_df_genes[i] <<- 'Hep 3'
  if (rownames(cluster_average_exp_df_scaled)[i] %in% Cholagiocytes_Markers) anno_df_genes[i] <<- 'Cholagiocyte'
  if (rownames(cluster_average_exp_df_scaled)[i] %in% c(General_Macrophage_Markers, Non_inf_mac_Markers)) anno_df_genes[i] <<- 'Non-inf Mac'
  if (rownames(cluster_average_exp_df_scaled)[i] %in% Inf_mac_Markers) anno_df_genes[i] <<- 'Inf Mac'
  if (rownames(cluster_average_exp_df_scaled)[i] %in% Mesenchymal_Markers) anno_df_genes[i] <<- 'Mesenchymal'
  if (rownames(cluster_average_exp_df_scaled)[i] %in% Endothelial_Markers) anno_df_genes[i] <<- 'Endothelial'
})
length(anno_df_genes)

##### ordering clusters
cluster_orders = as.character(c(CV_Hep, other_Hep, Periportal_Hep, Cholagiocytes, Non_inf_mac,
                   Inf_mac, Mesenchymal, Endothelial))
anno_df = data.frame(cluster=colnames(cluster_average_exp_df_scaled),annot=anno_df)

cluster_average_exp_df_scaled = cluster_average_exp_df_scaled[,cluster_orders]
anno_df = anno_df[match(colnames(cluster_average_exp_df_scaled), anno_df$cluster),]
head(anno_df)
dim(cluster_average_exp_df_scaled)



#### ordering genes
gene_orders = c(CV_Hep_Markers, other_Hep_Markers, Periportal_Hep_Markers, Cholagiocytes_Markers, 
                General_Macrophage_Markers, Non_inf_mac_Markers, Inf_mac_Markers, 
                Mesenchymal_Markers, Endothelial_Markers) 

gene_orders = unique(gene_orders[gene_orders %in% rownames(cluster_average_exp_df_scaled)])
anno_df_genes = data.frame(genes=row.names(cluster_average_exp_df_scaled),annot=anno_df_genes)


cluster_average_exp_df_scaled = cluster_average_exp_df_scaled[gene_orders,]
dim(cluster_average_exp_df_scaled)
anno_df_genes = anno_df_genes[match(row.names(cluster_average_exp_df_scaled), anno_df_genes$genes),]



dim(anno_df)
dim(anno_df_genes)

row.names(anno_df) = anno_df$cluster
row.names(anno_df_genes) = anno_df_genes$genes

head(anno_df)
head(anno_df_genes)

anno_df = anno_df[,-1,drop=F]
anno_df_genes = anno_df_genes[,-1,drop=F ]

dim(cluster_average_exp_df_scaled)

#mat_breaks <- seq(min(cluster_average_exp_df_scaled), max(cluster_average_exp_df_scaled), length.out = 10)
annot_info <- read.csv('figure_panel/Annotations_SingleNuc_Rat_June_8_2023.csv')
annot_info <- annot_info[1:34, 1:4]
colnames(annot_info)[1] = 'clusters'


show_col(color_df$colors)
c(CV_Hep_Markers, other_Hep_Markers, Periportal_Hep_Markers, Cholagiocytes_Markers, 
  General_Macrophage_Markers, Non_inf_mac_Markers, Inf_mac_Markers, 
  Mesenchymal_Markers, Endothelial_Markers) 


fc <- colorRampPalette(c("green", "darkgreen"))
annot_info.df = data.frame(table(annot_info$label))

fc <- colorRampPalette(c('palevioletred1'))
cvHep_c = fc(annot_info.df$Freq[annot_info.df$Var1=='Hep 1']) #cvHep

fc <- colorRampPalette(c('orchid3'))
Hep_c = fc(annot_info.df$Freq[annot_info.df$Var1=='Hep 2'])

fc <- colorRampPalette(c('maroon1'))
ppHep_c = fc(annot_info.df$Freq[annot_info.df$Var1=='Hep 3']) #ppHep

fc <- colorRampPalette(c('mistyrose'))
unknown_c = fc(annot_info.df$Freq[annot_info.df$Var1=='Unknown/High Mito'])

fc <- colorRampPalette(c('cyan')) #coral
chol_c = fc(annot_info.df$Freq[annot_info.df$Var1=='Cholangiocytes'])

infMac_c = '#117733'
nonInfMac_c = '#999933' #'#47A265'
LSEC_c = c('#FFE729', '#FFE729')
Stellate_c = '#6699CC'

anno_df_col = anno_df
anno_df_col$colors=c(cvHep_c, Hep_c, ppHep_c, chol_c, nonInfMac_c, infMac_c, Stellate_c, LSEC_c)

##### makeing the input list for the pheatmap annotation colors
show_col(unique(anno_df_col$colors))
anno_df_col_list = unique(anno_df_col$colors)
names(anno_df_col_list) = unique(anno_df_col$CellType)

colnames(anno_df) = 'CellType'
colnames(anno_df_genes) = 'CellType'
pheatmap(t(cluster_average_exp_df_scaled), cluster_rows = FALSE,
         cluster_cols = FALSE, 
         annotation_row = data.frame(anno_df),
         annotation_col = data.frame(anno_df_genes),
         color= inferno(20, direction = +1),  
         border_color= FALSE, 
         annotation_names_row=FALSE,
         annotation_names_col=FALSE,
         annotation_colors = list(CellType=anno_df_col_list),
         annotation_legend = TRUE
         )
?pheatmap
