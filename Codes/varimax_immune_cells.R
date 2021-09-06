####################################
### applying Varimax to the immune populations of the new samples ####
### immune cell clusters: cluster 4, 7, 8, 9, 11
### each of the samples are normalized first (all cells included), and then merged. 
### Afterwards immune cells are selected and then the immune subset will be scaled

## loading packages
source('Codes/Functions.R')
source('Codes/FactorAnalysisUtils.R')
Initialize()


##### Functions #####
clean_name <- function(x){
  x=x[-length(x)] 
  if(x[1]=='merged') x=x[-1]
  return(paste(x, collapse = '_'))
}


qc_variables <- c( "Library_size", "num_expressed_genes","mito_perc" )

check_qc_cor <- function(merged_samples, pca_embedding_df, main){
  df_cor <- data.frame(pca_embedding_df,
                       Library_size=merged_samples$nCount_RNA,
                       num_expressed_genes=merged_samples$nFeature_RNA, 
                       mito_perc=merged_samples$mito_perc)
  df_cor <- cor(df_cor[,!colnames(df_cor)%in% c('clusters')])
  df_cor_sub <- df_cor[colnames(df_cor)%in%qc_variables, !colnames(df_cor)%in%qc_variables]
  pheatmap(df_cor_sub, fontsize_row=10, main=main, color = magma(20,direction = 1))
}


## selecting the immune clusters indices
Immune_cells_UMI <- readRDS('Immune_cells_UMI.rds')
merged_samples <- readRDS('Objects/merged_samples_newSamples_MT-removed.rds')
merged_samples$clusters <- as.character(merged_samples$res.0.6)
UMIs_to_include <- colnames(merged_samples) %in% Immune_cells_UMI
colnames(merged_samples)[UMIs_to_include]
immune_cells_cluster_info <- merged_samples$clusters[UMIs_to_include]
  
### importing the raw files and applying normalization to the immune cells only
data_file_Lew <- 'Objects/rat_DA_M09_WK_008_3pr_v3/seur_QC_rat_DA_M09_WK_008_3pr_v3_mito_50_lib_1000.rds'
data_file_DA <- 'Objects/rat_LEW_M09_WK_009_3pr_v3/seur_QC_rat_LEW_M09_WK_009_3pr_v3_mito_50_lib_1000.rds'
sample_Lew <- readRDS(data_file_Lew)
sample_DA <- readRDS(data_file_DA)
samples_list <- c(sample_DA, sample_Lew)
names(samples_list) <- c('rat_DA_M09_WK_008', 'rat_LEW_M09_WK_009')

##### normalize and merge the samples ####

samples_raw <- lapply(samples_list, 
                      function(x) CreateSeuratObject(GetAssayData(x,assay='RNA')))

### (1) checking if separate normalization of immune cells in each sample changes the results -> YES!
samples_raw <- sapply(1:length(samples_list), 
                      function(i) CreateSeuratObject(GetAssayData(samples_list[[i]][,paste0(names(samples_list[i]),'_',colnames(samples_list[[i]])) %in% Immune_cells_UMI],
                                                                  assay='RNA')), simplify = F)
samples_scTransform <- lapply(samples_raw, function(x) {
  x = SCTransform(x, conserve.memory=F, verbose=T, return.only.var.genes=F,
                 variable.features.n = nrow(x[['RNA']]), 
                 ## according to the paper scaling is not recommended prior to PCA:
                 ## https://www.biorxiv.org/content/10.1101/576827v2
                 do.scale = FALSE, ### default value 
                 do.center = TRUE) ### default value ))
  x = CreateSeuratObject(GetAssayData(x, assay='SCT'))
  return(x)
}) 



all_merged_samples <- merge(samples_scTransform[[1]], samples_scTransform[[2]], 
                            add.cell.ids = names(samples_scTransform), 
                            project = "rat_data", 
                            merge.data = TRUE)

merged_samples <- CreateSeuratObject(GetAssayData(all_merged_samples, assay='RNA'))

merged_samples_sub <- merged_samples[,UMIs_to_include]
merged_samples_sub <- ScaleData(merged_samples_sub)




#### Removing all the Mt-genes for later steps
MIT_PATTERN = '^Mt-'
mito_genes_index <- grep(pattern = MIT_PATTERN, row.names(merged_samples_sub) )
merged_samples_sub[["mito_perc"]] <- PercentageFeatureSet(merged_samples_sub, features = mito_genes_index)

rownames(merged_samples_sub)[mito_genes_index]
merged_samples_sub <- merged_samples_sub[-mito_genes_index,]
grep(pattern = MIT_PATTERN, rownames(merged_samples_sub))



### (2) checking if separate normalization of immune cells in each sample changes the results 
### YES! but differences does not seem to be too much

#saveRDS(merged_samples_sub, 'merged_new_samples_SeperatlyNormed_immune.rds')
x = your_scRNAseq_data_object[rownames(your_scRNAseq_data_object) %in% rownames(merged_samples_sub),]
y = merged_samples_sub[rownames(merged_samples_sub) %in% rownames(your_scRNAseq_data_object),]
summary(as.vector(abs(GetAssayData(x) - GetAssayData(y))))


## add the initial clustering info
merged_samples_sub$init_clusters <- immune_cells_cluster_info

merged_samples <- readRDS('Results/new_samples/Immune_subclusters_updated.rds') 

## add the clustering resolution of the subclustering results
load('Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_ImmuneSub_c17Included.RData')
merged_samples_sub$immune_clusters = sCVdata_list$RNA_snn_res.1@Clusters
merged_samples_sub$immune_clusters = sCVdata_list$res.1@Clusters

###### Hep removed dataset - needs additional scaling for the PCA step
merged_samples_sub = readRDS('Results/new_samples/Immune_subclusters_labelCorrect_HepRemoved.rds')
merged_samples_sub = ScaleData(merged_samples_sub)


temp.df <- data.frame(str_split_fixed(colnames(merged_samples_sub), pattern = '_', n=6))
merged_samples_sub$sample_name <- temp.df$X2

  
### needed data for the varimax PCA 
merged_samples_sub <- RunPCA(merged_samples_sub, features = rownames(merged_samples_sub))  

loading_matrix = Loadings(merged_samples_sub, 'pca')
gene_exp_matrix = GetAssayData(merged_samples_sub, assay = 'RNA') #SCT
dim(gene_exp_matrix)
dim(loading_matrix)


rot_data <- get_varimax_rotated(gene_exp_matrix, loading_matrix)
rotatedLoadings <- rot_data$rotLoadings
scores <- data.frame(rot_data$rotScores)
colnames(scores) = paste0('PC_', 1:ncol(scores))
dim(scores)
dim(rotatedLoadings)
dim(merged_samples_sub) ### running PCA removes around 20 genes which have low variance
#saveRDS(rot_data, 'Results/new_samples/immune_varimax_results_cl17Inc.rds')
#saveRDS(rot_data, 'Results/new_samples/immune_varimax_results_cl17Inc_HepRemoved.rds')

merged_samples <- merged_samples_sub
##### Visualize the rotated PCs and evaluate mapping with clusters ######

### which PCs to evaluate in the rotated data

### PCA
pca_embedding_df <- data.frame(Embeddings(merged_samples,reduction = 'pca'))
PC_standardDev <- apply(pca_embedding_df, 2, sd)
elbow_plot(PC_standardDev, title = 'PCA on last 2 rat samples - immune cells')
ElbowPlot(merged_samples)
top_pc = 12

### varimax PCA
rot_standardDev <- apply(rot_data$rotScores, 2, sd)
elbow_plot(rot_standardDev, title = 'Varimax on last 2 rat samples - immune cells')
ndims = ncol(rot_data$rotScores)
rot_percVar = (rot_standardDev^2 * 100)/sum(rot_standardDev^2)
names(rot_percVar) <- paste0(1:(length(rot_percVar)))
perc_variance_threshold = 0.5
PCs_to_check <- as.numeric(names(rot_percVar[rot_percVar>perc_variance_threshold  ]))
PCs_to_check <- c(1:20,31:50)

#### PC and varimax correlation with technical covariates ####

check_qc_cor(merged_samples, pca_embedding_df[,1:top_pc], 
             main='PC embeddings correlation with technical covariates')

embedd_df_rotated <- data.frame(scores)
embedd_df_rotated_2 <- embedd_df_rotated[,PCs_to_check]
colnames(embedd_df_rotated_2) <- paste0('Var.PC_', PCs_to_check)
check_qc_cor(merged_samples, embedd_df_rotated_2, 
             main='Varimax-PC embeddings correlation with technical covariates')

pheatmap(cor(embedd_df_rotated[,c(PCs_to_check)]))
lapply(embedd_df_rotated, head)


pca_embedding_df$cluster <- merged_samples$init_clusters
##### Visualizing the PCA plots before the rotation #####

plot_dir = 'Plots/'
pdf(paste0(plot_dir, 'PCA_plots_mergedNewSamples_immuneCells_c17Inc.pdf'),width = 14,height=12)

for(i in 1:top_pc){
  df = data.frame(PC_1=pca_embedding_df$PC_1,
                  emb_val=pca_embedding_df[,i],
                  cluster=pca_embedding_df$cluster,
                  Library_size=merged_samples$nCount_RNA,
                  num_expressed_genes=merged_samples$nFeature_RNA,
                  sample_name=merged_samples$sample_name)
  
  p1=ggplot(df, aes(x=PC_1, y=emb_val, color=cluster))+geom_point(alpha=0.8)+
    theme_classic()+ylab(paste0('PC_',i))+scale_color_manual(values=colorPalatte)
  p2=ggplot(df, aes(x=PC_1, y=emb_val, color=Library_size))+geom_point(alpha=0.8)+
    theme_classic()+ylab(paste0('PC_',i))+scale_color_viridis(direction =-1 )
  p3=ggplot(df, aes(x=PC_1, y=emb_val, color=num_expressed_genes))+geom_point(alpha=0.8)+
    theme_classic()+ylab(paste0('PC_',i))+scale_color_viridis(direction =-1 )
  p4=ggplot(df, aes(x=PC_1, y=emb_val, color=sample_name))+geom_point(alpha=0.8)+
    theme_classic()+ylab(paste0('PC_',i))
  
  gridExtra::grid.arrange(p1,p2,p3,p4,ncol=2,nrow=2)
}
dev.off()


embedd_df_rotated = embedd_df_rotated_2
#colnames(embedd_df_rotated) <- paste0('Varimax_', 1:ncol(embedd_df_rotated))
##### Visualizing the PCA plots after the varimax rotation #####
pdf(paste0(plot_dir, 'VarimaxPCA_scTrans_mergedNewSamples_immuneCells_c17Inc.pdf'),width = 14,height=13) 

for(i in PCs_to_check){ 
  pc_num = i
  rot_df <- data.frame(Varimax_1=embedd_df_rotated$Varimax_1,
                       emb_val=embedd_df_rotated[,pc_num],
                       #cluster=as.character(merged_samples$final_cluster),
                       subcluster=as.character(merged_samples$cluster),
                       label = as.character(merged_samples$label),
                       Library_size=merged_samples$nCount_RNA,
                       num_expressed_genes=merged_samples$nFeature_RNA,
                       sample_name = merged_samples$sample_name)
  
  p1=ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=subcluster))+geom_point()+
    theme_classic()+ylab(paste0('Varimax_',i))+scale_color_manual(values=colorPalatte)
  p2=ggplot(rot_df, aes(x=subcluster, y=emb_val, fill=subcluster))+geom_boxplot()+
    theme_classic()+scale_fill_manual(values=colorPalatte)
  p3=ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=sample_name))+geom_point()+
    theme_classic()+ylab(paste0('Varimax_',i))
  p4=ggplot(rot_df, aes(x=sample_name, y=emb_val, fill=sample_name))+geom_boxplot()+theme_classic()
  
  gridExtra::grid.arrange(p1,p2,p3,p4,ncol=2,nrow=2)
}
dev.off()                      



##### checking the sample composition of each cluster ##### 

clusters_pred = as.character(merged_samples$immune_clusters)
df <- data.frame(sample_type = merged_samples$sample_name, 
                 cluster = as.character(clusters_pred))
rownames(df) = NULL
head(df)

#### based on sample-type ##### 
counts <- ddply(df, .(df$sample_type, df$cluster), nrow)
names(counts) <- c("sample_type", "cluster", "Freq")
counts$cluster= factor(counts$cluster, levels = as.character(0:(length(unique(df$cluster))-1)) ) 
ggplot(data=counts, aes(x=cluster, y=Freq, fill=sample_type)) +
  geom_bar(stat="identity",color='black')+theme_classic()+scale_fill_brewer(palette = "Set3")+
  ylab('Counts')

counts_split <- split( counts , f = counts$cluster )
counts_split_norm <- lapply(counts_split, function(x) {x$Freq=x$Freq/sum(x$Freq);x})
counts_norm <- do.call(rbind,counts_split_norm )
ggplot(data=counts_norm, aes(x=cluster, y=Freq, fill=sample_type)) +
  geom_bar(stat="identity",color='black')+theme_classic()+scale_fill_brewer(palette = "Set3")



##### checking res=1.0 clustering resolution #####
merged_samples <- FindNeighbors(merged_samples,
                                reduction="harmony",
                                dims=1:PC_NUMBER,
                                verbose=T)
merged_samples <- FindClusters(merged_samples,resolution=2.2,verbose=F)
table(as.character(merged_samples$RNA_snn_res.1), merged_samples$immune_clusters)


merged_samples <- RunHarmony(merged_samples, "sample_name",assay.use="RNA")

###### UMAP ###### 
a_gene_name = 'Alb'  #'Alb' 'Marco' Lyz2
a_gene_expression <- GetAssayData(merged_samples, assay='RNA')[a_gene_name,]
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=gene_expression))+geom_point()+
  theme_classic()+scale_color_viridis(direction = -1)+ggtitle(a_gene_name)

PC_NUMBER = top_pc
merged_samples <- RunUMAP(merged_samples, dims=1:PC_NUMBER, reduction="harmony")
df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap')[,2], 
                      library_size= merged_samples$nCount_RNA, 
                      n_expressed=merged_samples$nFeature_RNA,
                      mito_perc=merged_samples$mito_perc,
                      clusters= as.character(merged_samples$immune_clusters), 
                      sample_name=merged_samples$sample_name, 
                      strain=merged_samples$sample_name,
                      #gene_expression=a_gene_expression,
                      test_cluster=as.character(merged_samples$RNA_snn_res.1),
                      detectedGenes_More250=merged_samples$nFeature_RNA>250
)

pdf(paste0(plot_dir, 'umap_plots_newSamples_Mt_removed_immuneSub.pdf'))
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=library_size))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=n_expressed))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=mito_perc))+geom_point()+theme_classic()+scale_color_viridis(direction=1)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=gene_expression))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)+ggtitle(a_gene_name)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=sample_name))+geom_point(alpha=0.7)+theme_classic()+scale_color_manual(values = colorPalatte)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=strain))+geom_point(alpha=0.6)+theme_classic()+scale_color_manual(values = colorPalatte)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=clusters))+geom_point()+theme_classic()+scale_color_manual(values = colorPalatte)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=detectedGenes_More250))+geom_point(alpha=0.6)+theme_classic()+scale_color_manual(values = colorPalatte)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=test_cluster))+geom_point(alpha=0.7)+theme_classic()+scale_color_manual(values = colorPalatte)

dev.off()


###### tSNE ###### 
# Oat > pericentral
# Hpx > periportal
a_gene_name = 'Alb'  #'Alb' 'Marco' Lyz2
sum(a_gene_name %in% rownames(merged_samples))
a_gene_expression <- GetAssayData(merged_samples, assay='RNA')[a_gene_name,]

merged_samples <- RunTSNE(merged_samples,dims=1:PC_NUMBER, reduction="harmony")
df_tsne <- data.frame(tSNE_1=getEmb(merged_samples, 'tsne')[,1], 
                      tSNE_2=getEmb(merged_samples, 'tsne')[,2], 
                      library_size= merged_samples$nCount_RNA, 
                      n_expressed=merged_samples$nFeature_RNA,
                      mito_perc=merged_samples$mito_perc,
                      clusters=as.character(merged_samples$immune_clusters), 
                      sample_name=merged_samples$sample_name,
                      strain=merged_samples$sample_name,
                      gene_expression=a_gene_expression)

pdf(paste0(plot_dir, 'tsne_plots_newSamples_Mt_removed_immuneSub.pdf'))
ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=library_size))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)
ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=n_expressed))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)
ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=mito_perc))+geom_point()+theme_classic()+scale_color_viridis(direction=1)
ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=gene_expression))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)+ggtitle(a_gene_name)

ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=sample_name))+geom_point()+theme_classic()+scale_color_manual(values = colorPalatte)
ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=strain))+geom_point()+theme_classic()+scale_color_manual(values = colorPalatte)
ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2))+geom_point(aes(color=clusters),alpha=0.7,size=2)+theme_classic()+
  scale_color_manual(values = colorPalatte)
dev.off()




#### saving loadings #####
immune_loading_dir <- '~/RatLiver/Results/new_samples/rotated_loadings_immune/'

for(pc_index in PCs_to_check){ 
  print(pc_index)
  loading_df <- data.frame(genes=rownames(rotatedLoadings),values=rotatedLoadings[,pc_index])
  loading_df_ord <- loading_df[order(loading_df$values, decreasing = T),]
  rownames(loading_df_ord) <- NULL
  print(head(loading_df_ord, 25))
  print(tail(loading_df_ord, 25))
  write.table(loading_df_ord, 
              paste0(immune_loading_dir, 'rot_PC',pc_index,'_loadings.rnk'), 
              col.names = F, row.names = F, quote = F, sep = '\t')
}





##### annotate based on sub-cluster
cluster_of_interest = 'cluster_9'
subcluster_info <- readRDS(paste0('Results/subclusters/new_samples/', 
                                  cluster_of_interest ,'_subclust.rds'))
subclust.df = data.frame(umi=names(subcluster_info$res.1.48@Clusters), 
              subcluster=as.character(subcluster_info$res.1.48@Clusters))
pc_num = 19
rot_df <- data.frame(umi=rownames(embedd_df_rotated),
                     PC_1=embedd_df_rotated$PC_1,
                     emb_val=embedd_df_rotated[,pc_num],
                     cluster=merged_samples_sub$final_clusters,
                     Library_size=merged_samples_sub$nCount_RNA,
                     num_expressed_genes=merged_samples_sub$nFeature_RNA,
                     strain=merged_samples_sub$sample_name,
                     sample_name = merged_samples_sub$sample_name)
rot_df.2 <- merge(rot_df, subclust.df, by.x='umi', by.y='umi', all.x=T)
rot_df.2$subcluster_label <- ifelse(is.na(rot_df.2$subcluster), rot_df.2$cluster, 
                                    paste0(cluster_of_interest,'.sub_',rot_df.2$subcluster))
head(rot_df.2)

ggplot(rot_df.2, aes(x=PC_1, y=emb_val, color=subcluster_label))+geom_point()+
  theme_classic()+ylab(paste0('PC_',pc_num))+scale_color_manual(values=colorPalatte)

ggplot(rot_df, aes(x=PC_1, y=emb_val, color=cluster))+geom_point()+
  theme_classic()+ylab(paste0('PC_',pc_num))+scale_color_manual(values=colorPalatte)
ggplot(rot_df, aes(x=cluster, y=emb_val, fill=cluster))+geom_violin()+
  theme_classic()+scale_fill_manual(values=colorPalatte)
ggplot(rot_df, aes(x=PC_1, y=emb_val, color=strain))+geom_point()+
  theme_classic()+ylab(paste0('PC_',pc_num))
ggplot(rot_df, aes(x=strain, y=emb_val, fill=strain))+geom_violin()+theme_classic()




