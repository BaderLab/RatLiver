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
re <- readRDS('~/rat_sham_sn_data/standardQC_results/nebula_nonInfMac_subcluster_results.rds')

re.df = data.frame(re$summary)
re.df = re.df[,c(2,6,8)]
head(re.df[order(re.df$logFC_strainLEW, decreasing = T),],30)
head(re.df[order(re.df$p_strainLEW, decreasing = F),],20)



write.csv(re.df[order(re.df$p_strainLEW, decreasing = F),], 
          '~/rat_sham_sn_data/standardQC_results/nebula_nonInfMac_subcluster_results.csv')


head(re$summary)
re.df.ord = re.df[order(re.df$logFC_strainLEW, decreasing = T),]
re.df.ord = re.df.ord[re.df.ord$p_strainLEW < 0.05,]
table_vis = head(re.df.ord, 50)

table_vis = head(re.df[order(re.df$p_strainLEW, decreasing = F),],20)
table_vis = data.frame(gene=table_vis$gene, pval_LEW=table_vis$p_strainLEW, logFC_LEW=table_vis$logFC_strainLEW)
gridExtra::grid.table(table_vis)
dev.off()

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



######## correlation between the nebula results and varimax factors

### loading the geneset
rot_data <- readRDS('Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds') ## MT-removed
rotatedLoadings <- rot_data$rotLoadings
head(rotatedLoadings)

var_num = 15
### defining the signature based on varimax results
rotatedLoadings_ord = rotatedLoadings[order(rotatedLoadings[,var_num], decreasing = T),]
varimax_df = data.frame(gene=rownames(rotatedLoadings_ord) ,loading=rotatedLoadings_ord[,var_num])
head(varimax_df)
nebula_df = data.frame(gene=re.df$gene,
                       nebula_p=re.df$p_strainLEW, 
                       nebula_logFC=re.df$logFC_strainLEW) #

nebula_df$score = -log(re.df$p_strainLEW)*re.df$logFC_strainLEW

re.df.ord = re.df[order(re.df$logFC_strainLEW, decreasing = T),] ## lewis is positive
re.df.ord = re.df[order(re.df$p_strainLEW, decreasing = F),]

head(re.df.ord,30)

head(nebula_df)
merged_df = merge(nebula_df, varimax_df, by.x='gene', by.y='gene')

head(merged_df)
plot(merged_df$nebula_p, merged_df$loading)
plot(merged_df$nebula_logFC, merged_df$loading)
plot(merged_df$score, merged_df$loading)


ggplot(merged_df, aes(nebula_p, loading, color=nebula_logFC))+geom_point(size=1.5, alpha=0.6)+
  theme_classic()+scale_color_viridis(direction = -1)

ggplot(merged_df, aes(nebula_logFC, loading, color=nebula_p))+geom_point(size=1.5, alpha=0.6)+
  theme_classic()+scale_color_viridis(direction = -1)

ggplot(merged_df, aes(score, loading, color=nebula_p))+geom_point(size=1.5, alpha=0.6)+
  theme_classic()+scale_color_viridis(direction = 1)




############ Visualizing the DE gene list using a heatmap with strain columns
re <- readRDS('~/rat_sham_sn_data/standardQC_results/nebula_nonInfMac_subcluster_results.rds')
re.df = data.frame(re$summary)
re.df = re.df[,c(2,6,8)]
head(re.df[order(re.df$p_strainLEW, decreasing = F),],20)
re.df_ord = re.df[order(re.df$p_strainLEW, decreasing = F),]

merged_samples = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
macrophage_ids = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_macrophage_IDs.rds')
Resolution = 0.6
resolutions = Resolution
merged_samples <- FindClusters(merged_samples, resolution = Resolution, verbose = FALSE)
table(merged_samples$SCT_snn_res.0.6)

merged_samples <- merged_samples[,colnames(merged_samples)%in%macrophage_ids]
dim(merged_samples)
sample_data = list(count=GetAssayData(merged_samples@assays$RNA))

strain_names_types = names(table(merged_samples$strain))
### dividing the expression matrix based on the clusters
strain_expression <- sapply(1:length(strain_names_types), function(i){
  a_strain_name = strain_names_types[i]
  GetAssayData(merged_samples, 'data')[,merged_samples$strain == a_strain_name] 
}, simplify = F)

names(strain_expression) = strain_names_types
lapply(strain_expression, head)


## calculate the average expression of each gene in each cluster
strain_average_exp <- lapply(strain_expression, function(x){
  df = data.frame(average=rowSums(x)/ncol(x))
  return(df)
})
lapply(strain_average_exp, dim)

## Concatenate all the clusters together to make a matrix
strain_average_exp_df = do.call(cbind,strain_average_exp)
colnames(strain_average_exp_df) = names(strain_average_exp)
head(strain_average_exp_df)
strain_average_exp_df$gene = rownames(strain_average_exp_df)


###################################################
merged_info = merge(re.df, strain_average_exp_df, by.x='gene', by.y='gene', sort=F)
merged_info_ord = merged_info[match(re.df_ord$gene, merged_info$gene),]
head(merged_info_ord)
merged_info_ord_sub = merged_info_ord[1:25,c("gene", 'DA', 'LEW')]
row.names(merged_info_ord_sub) <- merged_info_ord_sub$gene
merged_info_ord_sub <- merged_info_ord_sub[,-1]

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(as.matrix(merged_info_ord_sub), n = 50)


pheatmap(merged_info_ord_sub, cluster_rows = F, cluster_cols = F)

pheatmap(merged_info_ord_sub, cluster_rows = F, cluster_cols = F,  
         color = inferno(length(mat_breaks) - 1),
         breaks= mat_breaks)


pheatmap(merged_info_ord_sub, cluster_rows = F, cluster_cols = F,
         color=colorRampPalette(c("navy", "white", "red"))(50))
