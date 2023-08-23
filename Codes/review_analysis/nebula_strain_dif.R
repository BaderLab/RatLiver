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


library(gridExtra)
library(grid)
grid.ftable <- function(d, padding = unit(4, "mm"), ...) {
  
  nc <- ncol(d)
  nr <- nrow(d)
  
  ## character table with added row and column names
  extended_matrix <- cbind(c("", rownames(d)),
                           rbind(colnames(d),
                                 as.matrix(d)))
  
  ## string width and height
  w <- apply(extended_matrix, 2, strwidth, "inch")
  h <- apply(extended_matrix, 2, strheight, "inch")
  
  widths <- apply(w, 2, max)
  heights <- apply(h, 1, max)
  
  padding <- convertUnit(padding, unitTo = "in", valueOnly = TRUE)
  
  x <- cumsum(widths + padding) - 0.5 * padding
  y <- cumsum(heights + padding) - padding
  
  rg <- rectGrob(x = unit(x - widths/2, "in"),
                 y = unit(1, "npc") - unit(rep(y, each = nc + 1), "in"),
                 width = unit(widths + padding, "in"),
                 height = unit(heights + padding, "in"))
  
  tg <- textGrob(c(t(extended_matrix)), x = unit(x - widths/2, "in"),
                 y = unit(1, "npc") - unit(rep(y, each = nc + 1), "in"),
                 just = "center")
  
  g <- gTree(children = gList(rg, tg), ...,
             x = x, y = y, widths = widths, heights = heights)
  
  grid.draw(g)
  invisible(g)
}

## check their tutorial
## https://github.com/lhe17/nebula

#merged_data <- readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_standardQC.rds')
#a_cluster = 9 ## cluster 9 is the macrophage cluster based on merged_data$SCT_snn_res.1 0f sham_sn_merged_standardQC.rds
#merged_data <- merged_data[,merged_data$SCT_snn_res.1==a_cluster]

merged_data <- readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
macrophage_ids = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_macrophage_IDs.rds')



Resolution = 2.5
merged_data <- FindClusters(merged_data, resolution = Resolution, verbose = FALSE)
table(merged_data$SCT_snn_res.2.5)
merged_data$clusters = merged_samples$SCT_snn_res.2.5

table(merged_data$cluster, merged_data$annot_IM) ## cluster 8 seems to be the macrophage population
table(merged_data$clusters, merged_data$annot_IM) 
#mac_cluster_num = 8
mac_cluster_num = 19

macrophage_ids = colnames(merged_data)[merged_data$clusters == mac_cluster_num]

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
re.df$score = -log10(re.df$p_strainLEW) * re.df$logFC_strainLEW
re.df = re.df[order(re.df$score, decreasing = T),]
head(re.df,20)


length(varimax_df$gene)
re.df$gene[1:20][re.df$gene[1:20] %in% varimax_df$gene[1:10]]


write.csv(re.df[order(re.df$p_strainLEW, decreasing = F),], 
          '~/rat_sham_sn_data/standardQC_results/nebula_nonInfMac_subcluster_results.csv')

table_vis = data.frame(gene=re.df$gene, pval_LEW=re.df$p_strainLEW, logFC_LEW=re.df$logFC_strainLEW)
table_vis$score = -log10(table_vis$pval_LEW) * table_vis$logFC_LEW
table_vis = table_vis[order(table_vis$score, decreasing = T),]


colnames(table_vis)=c('Gene', 'p value', 'logFC(LEW/DA)', 'score')
head(table_vis, 20)
gridExtra::grid.table(head(table_vis, 15))
dev.off()


table_vis <- table_vis[1:11,-4]
grid.newpage()
table_vis$`logFC(LEW/DA)` = round(table_vis$`logFC(LEW/DA)`, 2) 
table_vis$`p value` = formatC(table_vis$`p value`, format = "e", digits = 2)
row.names(table_vis) = table_vis$Gene
table_vis = table_vis[,-1]
grid.ftable(table_vis, gp = gpar(fill = rep(c("white", "white"), each = 6)),)






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
                      #is_cluster=ifelse(merged_data$SCT_snn_res.1 == a_cluster,'1','0'),
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
