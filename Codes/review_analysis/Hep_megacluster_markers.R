source('Codes/Functions.R')
source('~/XSpecies/Codes/convert_human_to_ortholog_functions.R')
library(Hmisc)
library(RColorBrewer)
library(viridis)
library(reshape2)
library(scales)
########################################################
############ Importing the nuc-seq data  #########

## importing the gene expression data
merged_samples = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
# merged_samples_all = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC_allfeatures.rds')

Resolution = 2.5
merged_samples <- FindClusters(merged_samples, resolution = Resolution, verbose = FALSE)

Resolution = 0.6
merged_samples <- FindClusters(merged_samples, resolution = Resolution, verbose = FALSE)

mapping_table = table(merged_samples$SCT_snn_res.2.5, merged_samples$SCT_snn_res.0.6)

#### defining the mega clusters based on 0.6 resolution
Hep0 = c(0, 4, 7, 10, 16, 25, 27) # 27 is the little tail
Hep1 = c(1, 6, 15, 17, 12, 31) # 31 is the tail
Hep2 = c(2, 5, 9, 21) #
Hep3 = c(3, 8,13, 23, 32 )


merged_samples$clusters = as.character(merged_samples$SCT_snn_res.2.5)
merged_samples$Hep_clusters = ifelse(merged_samples$clusters %in% as.character(Hep0), 'Hep0', merged_samples$clusters) 
merged_samples$Hep_clusters = ifelse(merged_samples$Hep_clusters %in% as.character(Hep1), 'Hep1', merged_samples$Hep_clusters) 
merged_samples$Hep_clusters = ifelse(merged_samples$Hep_clusters %in% as.character(Hep2), 'Hep2', merged_samples$Hep_clusters) 
merged_samples$Hep_clusters = ifelse(merged_samples$Hep_clusters %in% as.character(Hep3), 'Hep3', merged_samples$Hep_clusters) 
table(merged_samples$Hep_clusters)


cluster_num = '4'
gene_name = 'Cd5l'  #Sirpa
rownames(merged_samples)[grep(gene_name, rownames(merged_samples))] ## check if the gene is present in the dataset

df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap_h')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap_h')[,2], 
                      clusters=merged_samples$clusters,
                      expression_val=GetAssayData(merged_samples)[gene_name,],
                      a_cluster=merged_samples$cluster==cluster_num,
                      umi=colnames(merged_samples),
                      Hep_clusters=merged_samples$Hep_clusters)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=Hep_clusters))+geom_point(alpha=0.7, size=2)+
  theme_classic()+scale_color_manual(name='Hep_clusters',values = colorPalatte)+
  theme(text = element_text(size=16),legend.title = element_blank())#


ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=expression_val))+
  geom_point(alpha=0.6,size=1.2)+
  scale_color_viridis('Expression\nValue', direction = -1)+theme_classic()+
  theme(text = element_text(size=16.5),
        plot.title = element_text(hjust = 0.5),
        legend.title=element_text(size = 10))+
  ggtitle(paste0(gene_name)) 

##### adding the final annotations to the dataframe
#annot_info <- read.csv('figure_panel/Annotations_SingleNuc_Rat_May_31_2023.csv')
#annot_info <- read.csv('figure_panel/Annotations_SingleNuc_Rat_June_8_2023.csv')
annot_info <- read.csv('figure_panel/Annotations_SingleNuc_Rat_June_21_2023.csv')

#annot_info <- read.csv('figure_panel/Annotations_SingleNuc_Rat_June_20_2023.csv')
annot_info <- annot_info[1:35, 1:4]
colnames(annot_info)[1] = 'clusters'

# clusters_not_included 
names(table(merged_samples$SCT_snn_res.2.5))[!names(table(merged_samples$SCT_snn_res.2.5)) %in% as.character(annot_info$clusters)]


########## merging df_umap with annot_info data.frame
df_umap = merge(df_umap, annot_info, by.x='clusters', by.y='clusters', all.x=T, order=F)
#### re-ordering the rows based on UMIs of the initial data
df_umap <- df_umap[match(colnames(merged_samples),df_umap$umi),]
df_umap$umi == colnames(merged_samples)
sum(df_umap$umi != colnames(merged_samples))
df_umap$label_clust = paste0(df_umap$label, ' (',df_umap$clusters, ')')

annot_info$label_clust = paste0(annot_info$label, ' (',annot_info$clusters, ')')
library("RColorBrewer")
title = '' # 'Mt-genes removed'
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=label))+geom_point(alpha=0.8, size=1)+
  theme_classic()+scale_color_manual(values = colorPalatte)+
  #scale_color_brewer(palette = 'Dark2')+ #name='clusters',
  theme(text = element_text(size=15),legend.title = element_blank())#

table(df_umap$clusters, df_umap$label)
##### checking if the order of cells in the merged samples and df object are the same
sum(colnames(merged_samples) != df_umap$umi)
merged_samples$final_label = df_umap$label


table(merged_samples$final_label )
dim(merged_samples)
DE_groups = c('Hep 3', 'Hep 1')
merged_samples_sub = merged_samples[,merged_samples$final_label %in% DE_groups]
dim(merged_samples_sub)
Idents(merged_samples_sub) = merged_samples_sub$final_label
Hep_levels= levels(Idents(merged_samples_sub))
Hep_levels
DE_groups
Hep_markers = FindMarkers(merged_samples_sub, 
            ident.1=DE_groups[1],
            ident.2=DE_groups[2],
            logfc.threshold = 0,
            min.pct=0,
            min.cells.feature = 1,
            min.cells.group = 1)


#names(Hep_markers) <- DE_groups


#file_name = gsub(' ', '',paste0('~/rat_sham_sn_data/standardQC_results/Hep_DEs/',DE_groups[1],
#                                '_vs_' ,DE_groups[2],'_megacluster_markers_updatedClusters.csv'))
#file_name = gsub(' ', '',paste0('~/rat_sham_sn_data/standardQC_results/Hep_DEs/',DE_groups[1],
#                                '_vs_' ,DE_groups[2],'_megacluster_markers_oldClusters.csv'))
file_name = gsub(' ', '',paste0('~/rat_sham_sn_data/standardQC_results/Hep_DEs/',DE_groups[1],
                                '_vs_' ,DE_groups[2],'_megacluster_markers_updatedClusters_hep1remove12and26.csv'))

file_name
write.csv(Hep_markers, file_name, col.names = T, row.names = T, quote = F)

#saveRDS(Hep_markers, paste0('~/rat_sham_sn_data/standardQC_results/Hep1_Hep3_megacluster_markers.rds'))
#write.csv(Hep_markers, '~/rat_sham_sn_data/standardQC_results/Hep1_Hep3_megacluster_markers.csv', col.names = T, row.names = T, quote = F)
#write.csv(Hep_markers, '~/rat_sham_sn_data/standardQC_results/Hep3_Hep1_megacluster_markers.csv', col.names = T, row.names = T, quote = F)


