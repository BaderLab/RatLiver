source('Codes/Functions.R')
source('~/XSpecies/Codes/convert_human_to_ortholog_functions.R')
library(Hmisc)
library(RColorBrewer)
library(viridis)
library(reshape2)
########################################################
############ Importing the nuc-seq data  #########

## importing the gene expression data
merged_samples = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')

Resolution = 2.5
merged_samples <- FindClusters(merged_samples, resolution = Resolution, verbose = FALSE)
merged_samples$clusters = merged_samples$SCT_snn_res.2.5

df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap_h')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap_h')[,2], 
                      clusters=merged_samples$clusters,
                      umi=colnames(merged_samples))

##### adding the final annotations to the dataframe
annot_info <- read.csv('figure_panel/Annotations_SingleNuc_Rat_June_8_2023.csv')
colnames(annot_info)[1] = 'clusters'

########## merging df_umap with annot_info data.frame
df_umap = merge(df_umap, annot_info, by.x='clusters', by.y='clusters', all.x=T, order=F)
df_umap <- df_umap[match(colnames(merged_samples),df_umap$umi),]
data.frame(table(df_umap$label))
df_umap$label2 <- paste0(df_umap$label, ' (', df_umap$clusters, ')')
table(df_umap$label2)
write.csv(data.frame(table(df_umap$label2)), 'figure_panel/Map_CellType_Freq/snRNAseq_map_cellFreq.csv')




####################################################################
old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)
merged_samples <- your_scRNAseq_data_object
merged_samples$cluster = as.character(sCVdata_list$res.0.6@Clusters)
##### adding the final annotations to the dataframe
set1_info <- read.csv('figure_panel/set-1-final-info.csv')
colnames(set1_info)[1] = 'clusters'
set1_info$Annotation = gsub(pattern = 'Naive', replacement ='ab' ,set1_info$Annotation)

df_umap <- data.frame(umi=colnames(merged_samples),clusters=merged_samples$cluster )

set1_info$clusters = as.character(set1_info$clusters)
df_umap = merge(df_umap, set1_info[,1:3], by.x='clusters', by.y='clusters', all.x=T, order=F)
#### re-ordering the rows based on UMIs of the initial data
df_umap <- df_umap[match(colnames(merged_samples),df_umap$umi),]
df_umap$umi == colnames(merged_samples)
sum(df_umap$umi != colnames(merged_samples))
data.frame(table(df_umap$label))
write.csv(data.frame(table(df_umap$label)), 'figure_panel/Map_CellType_Freq/TLH_map_cellFreq.csv')

rm(list=ls())
###### set-2 updated labels
set_2_cl_ord = c( "Hep (1)", "Hep (3)", "Hep (9)", "Hep (14)", "Hep (15)", "Hep (0)" , "Hep (2)",  "LSEC (6)" , 
                  "LSEC (4)", "Stellate (16)", "Non-Inflammatory Mac (8)",   "Non-Inflammatory Mac (13)", 
                  "Inflammatory Mac (11)",  "gd T cell (7)" , "Naive T cell (10)", "Mature B cell (12)", 
                  "pDC (17)", "Erythroid (5)")

new_data_scCLustViz_object <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_labelCor.RData"
set1_info <- read.csv('figure_panel/set-2-final-info-updated.csv')
load(new_data_scCLustViz_object)
merged_samples <- your_scRNAseq_data_object
merged_samples$cluster = as.character(sCVdata_list$res.0.6@Clusters)

colnames(set1_info)[1] = 'clusters'
set1_info$Annotation = gsub(pattern = 'Naive', replacement ='ab' ,set1_info$Annotation)

df_umap <- data.frame(umi=colnames(merged_samples),clusters=merged_samples$cluster )

set1_info$clusters = as.character(set1_info$clusters)
df_umap = merge(df_umap, set1_info[,1:3], by.x='clusters', by.y='clusters', all.x=T, order=F)
#### re-ordering the rows based on UMIs of the initial data
df_umap <- df_umap[match(colnames(merged_samples),df_umap$umi),]
df_umap$umi == colnames(merged_samples)
sum(df_umap$umi != colnames(merged_samples))
data.frame(table(df_umap$label))
write.csv(data.frame(table(df_umap$label)), 'figure_panel/Map_CellType_Freq/ImmuneEnriched_map_cellFreq.csv')




