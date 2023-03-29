##### merging the TLH and snRNA-seq maps
### importing the normalized samples
#### scale them individually 

old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)
TLH_data <- your_scRNAseq_data_object
TLH_data@meta.data$cluster = as.character(sCVdata_list$res.0.6@Clusters)
TLH_data$strain = unlist(lapply(str_split(TLH_data$orig.ident, '_'), '[[', 2))
dim(TLH_data)
head(TLH_data)
### adding the cell-type annotation ###
annotation_TLH = read.csv('TLH_annotation.csv')
meta_data = TLH_data@meta.data
meta_data$cell_id = rownames(meta_data)
meta_data2 = merge(meta_data, annotation_TLH, by.x='cluster', by.y='cluster', all.x=T, all.y=F)
meta_data2 = meta_data2[match(meta_data$cell_id, meta_data2$cell_id),]
sum(meta_data$cell_id != colnames(TLH_data))
sum(meta_data2$cell_id != colnames(TLH_data))
meta_data2$strain = gsub(meta_data2$strain, pattern = 'rat_', replacement = '')
rownames(meta_data2) = meta_data2$cell_id
meta_data2 <- meta_data2[,colnames(meta_data2) != 'cell_id']
head(meta_data2)
dim(meta_data2)
table(meta_data2$strain)
TLH_data$sample_name = TLH_data$orig.ident 

TLH_data$annotation = meta_data2$annotation
TLH_data = ScaleData(TLH_data)

table(TLH_data$annotation)
######################################################

sn_sham_data = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
Resolution = 0.6
sn_sham_data <- FindClusters(sn_sham_data, resolution = Resolution, verbose = FALSE)
table(sn_sham_data$SCT_snn_res.0.6)
sn_sham_data$annotation = sn_sham_data$annot_IM_g

sn_sham_data = ScaleData(sn_sham_data)

######################################################

merged_data <- merge(sn_sham_data, TLH_data, # 
                     add.cell.ids = c('snMap','scTLH'), 
                     merge.data = TRUE)
dim(merged_data)

Idents(merged_data) = merged_data$sample_name

table(merged_data$sample_name)


### run the pipeline without this line and see how it would differ - could not run pca and find variable genes
#### scaling the merged data instead of SCTransform does not solve the issue
### current pipeline: we are normlizing the individual samples first - merging them - normalizing the merged sample 
merged_data = SCTransform(merged_data, vst.flavor = "v2", verbose = TRUE) 
merged_data <- ScaleData(merged_data)
#merged_data <- FindVariableFeatures(merged_data) #SCT assay is comprised of multiple SCT models

merged_data <- RunPCA(merged_data,verbose=T)
plot(100 * merged_data@reductions$pca@stdev^2 / merged_data@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")

#### UMAP on PC components and clustering
merged_data <- RunUMAP(merged_data, reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.6, verbose = FALSE)

#### UMAP on corrected PC components and clustering
merged_data <- RunHarmony(merged_data, group.by.vars = "sample_name")
merged_data <- RunUMAP(merged_data, reduction = "harmony", dims = 1:30, reduction.name = "umap_h")  %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.6, verbose = FALSE)

merged_data <- readRDS('~/rat_sham_sn_data/standardQC_results/TLH_snMap_merged.rds')
merged_data$annotation2 = ifelse(is.na(merged_data$annotation), merged_data$annot_IM_g, merged_data$annotation)
table(merged_data$annotation2)
merged_data$annotation2_g = gsub("\\s*\\([^\\)]+\\)","",as.character(merged_data$annotation2))

umap_df = data.frame(annotation=merged_data$annotation2,
                     annotation_g=merged_data$annotation2_g,
                     sample_name =merged_data$sample_name,
                     strain=merged_data$strain,
                     cluster=merged_data$SCT_snn_res.0.6,
                     data=ifelse(grepl(pattern = '*_CST', merged_data$sample_name), 'snRNA-seq', 'scRNA-seq'),
                     Embeddings(merged_data, reduction = 'umap'),
                     Embeddings(merged_data, reduction = 'umap_h'))


##### making UMAP plots ####
ggplot(umap_df, aes(x=UMAP_1, y=UMAP_2, color=sample_name))+geom_point(size=1,alpha=0.6)+
  theme_classic()+ scale_color_manual(values = colorPalatte)
ggplot(umap_df, aes(x=UMAP_1, y=UMAP_2, color=data))+geom_point(size=1.32,alpha=0.3)+theme_classic()
ggplot(umap_df, aes(x=UMAP_1, y=UMAP_2, color=annotation))+geom_point(size=1,alpha=0.6)+
  theme_classic()+ scale_color_manual(values = colorPalatte)
ggplot(umap_df, aes(x=UMAP_1, y=UMAP_2, color=annotation_g))+geom_point(size=1,alpha=0.6)+
  theme_classic()+ scale_color_manual(values = colorPalatte)

ggplot(umap_df, aes(x=umap_h_1, y=umap_h_2, color=sample_name))+geom_point(size=1,alpha=0.3)+
  theme_classic()+ scale_color_manual(values = colorPalatte)
ggplot(umap_df, aes(x=umap_h_1, y=umap_h_2, color=data))+geom_point(size=1.32,alpha=0.3)+theme_classic()
ggplot(umap_df, aes(x=umap_h_1, y=umap_h_2, color=annotation))+geom_point(size=1,alpha=0.6)+
  theme_classic()+ scale_color_manual(values = colorPalatte)
ggplot(umap_df, aes(x=umap_h_1, y=umap_h_2, color=annotation_g))+geom_point(size=1,alpha=0.6)+
  theme_classic()+ scale_color_manual(values = colorPalatte)
ggplot(umap_df, aes(x=umap_h_1, y=umap_h_2, color=cluster))+geom_point(size=1,alpha=0.6)+
  theme_classic()+ scale_color_manual(values = colorPalatte)

ggplot(umap_df, aes(x=UMAP_1, y=UMAP_2, color=cluster))+geom_point(size=3,alpha=0.9)+
  theme_classic()+ scale_color_manual(values = mycolors)




#### barplots of the contribution of each map to general clusters  ##### 
counts <- ddply(umap_df, .(umap_df$data, umap_df$cluster), nrow)
names(counts) <- c("data", "cluster", "Freq")
counts$cluster= factor(counts$cluster, levels = as.character(0:(length(table(counts$cluster))-1)) ) 


ggplot(data=counts, aes(x=cluster, y=Freq, fill=data)) +
  geom_bar(stat="identity",color='black')+theme_classic()+
  ylab('Counts')+xlab('Clusters')+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=10,angle=90,color='black'),
        legend.title = element_blank())+xlab('')+scale_fill_manual(values=colorPalatte)

counts_split <- split( counts , f = counts$cluster )
counts_split_norm <- lapply(counts_split, function(x) {x$Freq=x$Freq/sum(x$Freq);x})
counts_norm <- do.call(rbind,counts_split_norm )
ggplot(data=counts_norm, aes(x=cluster, y=Freq, fill=data)) +
geom_bar(stat="identity",color='black')+theme_classic()+
  ylab('Contribution of map per cluster (%)')+xlab('Cluster')+
  theme(text = element_text(size=15),
    axis.text.x = element_text(size=11,angle=90,color='black'),
    legend.title = element_blank()) +  xlab('')+scale_fill_manual(values=colorPalatte)




