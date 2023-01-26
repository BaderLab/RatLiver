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


merged_samples <- readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
##################################################################
######## Subclustering the macrophage cluster

mac_data = merged_samples[, merged_samples$cluster=='7']
mac_data_meta = mac_data@meta.data
DefaultAssay(mac_data) <- 'RNA'
mac_data@assays$SCT <- NULL


mac_data = SCTransform(mac_data, vst.flavor = "v2", verbose = TRUE)
mac_data <- RunPCA(mac_data,verbose=T)
plot(100 * mac_data@reductions$pca@stdev^2 / mac_data@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")
mac_data <- RunHarmony(mac_data, group.by.vars = "sample_name", assay.use="RNA")

res = 1.2
mac_data <- RunUMAP(mac_data, reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = res, verbose = FALSE)

markers1 = c('Ptprc','Cd68', 'Cd163', 'Mrc1', 'Clec4f', 'Il18', 'Ccl3', 'Timp2', 'Ly6al', 'Hmox1', 'Siglec5')
markers2 = c('Lyz2', 'Xcr1', 'Clec10a')
markers3 = c('Clec9a', 'Xcr1', 'Clec10a')

markers = c(markers1, markers2, markers3)
i = 15
gene_name = markers[i] 
gene_name
df_umap <- data.frame(UMAP_1=getEmb(mac_data, 'umap')[,1], 
                      UMAP_2=getEmb(mac_data, 'umap')[,2], 
                      library_size= mac_data$nCount_RNA, 
                      mito_perc=mac_data$mito_perc, 
                      n_expressed=mac_data$nFeature_RNA,
                      cluster=mac_data$seurat_clusters, 
                      cell_status = mac_data$cell_status,
                      nuclear_fraction=mac_data$nuclear_fraction, 
                      Alb=GetAssayData(mac_data)['Alb',], 
                      a_gene = GetAssayData(mac_data)[gene_name,],
                      sample_name = mac_data$sample_name, 
                      strain = mac_data$strain, 
                      umi=colnames(mac_data))

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=a_gene))+geom_point(size=1.5, alpha=0.6)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(gene_name)



ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cluster))+geom_point(size=1,alpha=0.6)+
  theme_classic()+scale_color_manual(values = c(colorPalatte))+ggtitle(paste0('resolution: ', res))#colorPalatte
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=sample_name))+geom_point(size=1.3,alpha=0.7)+theme_classic()


data.frame(table(mac_data$SCT_snn_res.1.2))


df <- data.frame(sample_type = mac_data$sample_name, 
                 cluster = as.character(mac_data$SCT_snn_res.1.2))
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


Idents(mac_data) = mac_data$strain
strain_mac_makers = FindMarkers(mac_data, 
                                ident.1 = 'LEW', 
                                ident.2 = 'DA')

strain_mac_makers_ord = strain_mac_makers[order(strain_mac_makers$avg_log2FC, decreasing = T),]
strain_mac_makers_ord = strain_mac_makers[order(strain_mac_makers$p_val_adj, decreasing = T),]
head(strain_mac_makers_ord, 30)




library(devtools)
install_github("lhe17/nebula")

library(nebula)
re = nebula(count = input_data2[,1:num_cells_to_include],
            id = designMatrix2$strainLEW,
            pred=designMatrix2[1:num_cells_to_include,][,c(1,2)])



