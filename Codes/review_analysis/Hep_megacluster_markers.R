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
table(merged_samples$SCT_snn_res.2.5)
merged_samples$clusters = merged_samples$SCT_snn_res.2.5


cluster_num = '4'
gene_name = 'Cd5l'  #Sirpa
rownames(merged_samples)[grep(gene_name, rownames(merged_samples))] ## check if the gene is present in the dataset

df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap_h')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap_h')[,2], 
                      clusters=merged_samples$clusters,
                      expression_val=GetAssayData(merged_samples)[gene_name,],
                      a_cluster=merged_samples$cluster==cluster_num,
                      umi=colnames(merged_samples))



ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=expression_val))+
  geom_point(alpha=0.6,size=1.2)+
  scale_color_viridis('Expression\nValue', direction = -1)+theme_classic()+
  theme(text = element_text(size=16.5),
        plot.title = element_text(hjust = 0.5),
        legend.title=element_text(size = 10))+
  ggtitle(paste0(gene_name)) 

##### adding the final annotations to the dataframe
#annot_info <- read.csv('figure_panel/Annotations_SingleNuc_Rat_May_31_2023.csv')
annot_info <- read.csv('figure_panel/Annotations_SingleNuc_Rat_June_8_2023.csv')
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

title = '' # 'Mt-genes removed'
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=label))+geom_point(alpha=0.8, size=1)+
  theme_classic()+#scale_color_manual(values = Colors)+ #name='clusters',
  theme(text = element_text(size=15),legend.title = element_blank())#

##### checking if the order of cells in the merged samples and df object are the same
sum(colnames(merged_samples) != df_umap$umi)
merged_samples$final_label = df_umap$label


dim(merged_samples)
merged_samples_sub = merged_samples[,merged_samples$final_label %in% c('Hep 1' ,'Hep 3')]
dim(merged_samples_sub)
Idents(merged_samples_sub) = merged_samples_sub$final_label
levels= levels(Idents(merged_samples_sub))

Hep_markers = FindMarkers(merged_samples_sub, 
            ident.1=levels[1],
            ident.2=levels[2],
            logfc.threshold = 0,
            min.pct=0,
            min.cells.feature = 1,
            min.cells.group = 1)


names(Hep_markers) <- levels
Cluster_markers = Hep_markers


saveRDS(Hep_markers, paste0('~/rat_sham_sn_data/standardQC_results/Hep1_Hep3_megacluster_markers.rds'))
write.csv(Hep_markers, '~/rat_sham_sn_data/standardQC_results/Hep1_Hep3_megacluster_markers.csv', 
          col.names = T, row.names = T, quote = F)


