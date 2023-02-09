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
library(nebula)

## check their tutorial
## https://github.com/lhe17/nebula

#merged_data <- readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_standardQC.rds')
#a_cluster = 9 ## cluster 9 is the macrophage cluster based on merged_data$SCT_snn_res.1 0f sham_sn_merged_standardQC.rds
#merged_data <- merged_data[,merged_data$SCT_snn_res.1==a_cluster]

merged_data <- readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
macrophage_ids = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_macrophage_IDs.rds')

table(merged_samples$cluster, merged_samples$annot_IM) ## cluster 8 seems to be the macrophage population
mac_cluster_num = 8

merged_data <- merged_data[,colnames(merged_data)%in%macrophage_ids]
dim(merged_data)
sample_data = list(count=GetAssayData(merged_data@assays$RNA))
sample_data$sid = merged_data$sample_name
sample_data$offset = merged_data$nCount_RNA

#sample_data$pred = data.frame(strain=as.factor(merged_data$strain),
#                cluster=as.factor(merged_data$SCT_snn_res.1))

sample_data$pred = data.frame(strain=as.factor(merged_data$strain))

sample_data$count[1:5,1:5]
head(sample_data$sid)
table(sample_data$sid)
head(sample_data$pred)

df = model.matrix(~strain, data=sample_data$pred)
data_g = group_cell(count=sample_data$count,id=sample_data$sid,pred=df)


### use library size as an offset if you're using the count data. If not, each cell will be treated equally
re = nebula(sample_data$count,sample_data$sid,pred=df,offset=sample_data$offset) 
re.df = data.frame(re$summary)
head(re.df[order(re.df$logFC_strainLEW, decreasing = T),],30)
head(re.df[order(re.df$p_strainLEW, decreasing = F),],30)

head(re$summary)
re.df.ord = re.df[order(re.df$logFC_strainLEW, decreasing = T),]
re.df.ord = re.df.ord[re.df.ord$p_strainLEW < 0.05,]
head(re.df.ord, 30)
df_umap <- data.frame(UMAP_1=getEmb(merged_data, 'umap')[,1], 
                      UMAP_2=getEmb(merged_data, 'umap')[,2], 
                      UMAP_h_1=getEmb(merged_data, 'umap_h')[,1], 
                      UMAP_h_2=getEmb(merged_data, 'umap_h')[,2], 
                      library_size= merged_data$nCount_RNA, 
                      mito_perc=merged_data$mito_perc, 
                      n_expressed=merged_data$nFeature_RNA,
                      cluster=merged_data$SCT_snn_res.1, 
                      cell_status = merged_data$cell_status,
                      nuclear_fraction=merged_data$nuclear_fraction, 
                      Alb=GetAssayData(merged_data)['Alb',], 
                      #gene=GetAssayData(merged_data)[gene_name,],
                      is_cluster=ifelse(merged_data$SCT_snn_res.1 == a_cluster,'1','0'),
                      sample_name = merged_data$sample_name, 
                      strain = merged_data$strain)

ggplot(df_umap, aes(x=UMAP_h_1, y=UMAP_h_2, color=cluster))+geom_point(alpha=0.8)+theme_classic()+scale_color_manual(values = c(colorPalatte, 'black') )#+ggtitle(paste0('res: ', Resolution))
ggplot(df_umap, aes(x=UMAP_h_1, y=UMAP_h_2, color=is_cluster))+geom_point(alpha=0.8)+theme_classic()

ggplot(df_umap, aes(x=UMAP_h_1, y=UMAP_h_2, color=gene))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(gene_name)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=gene))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(gene_name)
