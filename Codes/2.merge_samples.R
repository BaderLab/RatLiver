source('Codes/Functions.R')
source('Codes/FactorAnalysisUtils.R')
Initialize()
library(plyr)

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


##### Load the new samples' data #####
sample_name_Lew <- 'rat_LEW_M09_WK_009_3pr_v3'
sample_name_DA <- 'rat_DA_M09_WK_008_3pr_v3'
data_file_Lew <- 'Objects/rat_DA_M09_WK_008_3pr_v3/seur_QC_rat_DA_M09_WK_008_3pr_v3_mito_50_lib_1000.rds'
data_file_DA <- 'Objects/rat_LEW_M09_WK_009_3pr_v3/seur_QC_rat_LEW_M09_WK_009_3pr_v3_mito_50_lib_1000.rds'

sample_Lew <- readRDS(data_file_Lew)
sample_DA <- readRDS(data_file_DA)
samples_list <- c(sample_DA, sample_Lew)
names(samples_list) <- c('rat_DA_M09_WK_008', 'rat_LEW_M09_WK_009')


##### Load the old samples' data #####
data_file_Lew_01 <- 'Objects/rat_Lew_01/seur_QC_rat_Lew_01_mito_40_lib_2000.rds'
data_file_Lew_02 <- 'Objects/rat_Lew_02/seur_QC_rat_Lew_02_mito_40_lib_2000.rds'
data_file_DA_01 <- 'Objects/rat_DA_01_reseq/seur_QC_rat_DA_01_reseq_mito_30_lib_1500.rds'
data_file_DA_M_10WK <- 'Objects/rat_DA_M_10WK_003/seur_QC_rat_DA_M_10WK_003_mito_20_lib_2000.rds'

data_file_Lew_01 <- 'Objects/rat_Lew_01/seur_QC_rat_Lew_01_mito_40_lib_1500.rds'
data_file_Lew_02 <- 'Objects/rat_Lew_02/seur_QC_rat_Lew_02_mito_40_lib_1500.rds'
data_file_DA_01 <- 'Objects/rat_DA_01_reseq/seur_QC_rat_DA_01_reseq_mito_40_lib_1500.rds'
data_file_DA_M_10WK <- 'Objects/rat_DA_M_10WK_003/seur_QC_rat_DA_M_10WK_003_mito_40_lib_1500.rds'



samples_list <- c(data_file_Lew_01, data_file_Lew_02, data_file_DA_01, data_file_DA_M_10WK)
samples_list <- lapply(samples_list, readRDS)
names(samples_list) <- c('rat_Lew_01', 'rat_Lew_02', 'rat_DA_01_reseq', 'rat_DA_M_10WK_003')



##### normalize and merge the samples ####
samples_raw <- lapply(samples_list, function(x) CreateSeuratObject(GetAssayData(x, assay='RNA')))

samples_scTransform <- lapply(samples_raw, function(x) {
  x= SCTransform(x, conserve.memory=F, verbose=T, return.only.var.genes=F,
                 variable.features.n = nrow(x[['RNA']]), 
                 ## according to the paper scaling is not recommended prior to PCA:
                 ## https://www.biorxiv.org/content/10.1101/576827v2
                 do.scale = FALSE, ### default value 
                 do.center = TRUE) ### default value ))
  x = CreateSeuratObject(GetAssayData(x, assay='SCT'))
  return(x)
}) 

all_merged_samples <- merge(samples_scTransform[[1]], c(samples_scTransform[[2]], samples_scTransform[[3]], samples_scTransform[[4]] ), # 
                            add.cell.ids = names(samples_scTransform), 
                            project = "rat_data", 
                            merge.data = TRUE)

###  selecting the count data for cNMF
#samples_scTransform <- lapply(samples_list, function(x) CreateSeuratObject(GetAssayData(x, 'counts')))
# merged_samples <- CreateSeuratObject(GetAssayData(all_merged_samples, 'counts'))

merged_samples <- CreateSeuratObject(GetAssayData(all_merged_samples, assay='RNA'))
dim(merged_samples)
rowSums(merged_samples)

merged_samples <- ScaleData(merged_samples)
### merge the count files
### normalize and scale the count files
### run PCA on the merged data set

## number of feature in the raw rat data: 32883
## number of features in the normalized LEW data: 12410
## number of features in the normalized DA data: 13565

## 1351 gene in the Dark Agouti sample are either not expressed or have uniform expression in Lewis rat
## 196 gene in the Lewis sample are either not expressed or have uniform expression in Dark Agouti rat

### how exactly is scTransform working in this case???

### list of genes which are not included in the other sample ###
### all these genes are included in the merged sample, and although their values are non-zero(although low) in their original datasets,
### in the final merged dataset, their value will be set as zero 
DA_genes_not_in_Lew <- rownames(samples_scTransform$rat_DA_M09_WK_008)[!rownames(samples_scTransform$rat_DA_M09_WK_008) %in% 
                                                                         rownames(samples_scTransform$rat_LEW_M09_WK_009)]
Lew_genes_not_in_DA <- rownames(samples_scTransform$rat_LEW_M09_WK_009)[!rownames(samples_scTransform$rat_LEW_M09_WK_009) %in% 
                                                                          rownames(samples_scTransform$rat_DA_M09_WK_008)]


#### add the sample names ##### 
sample_names <- str_split(string = colnames(merged_samples), pattern = '_')
sample_names <- unlist(lapply(sample_names, clean_name))
merged_samples$sample_name <- sample_names
### ratio of the cells from each sample (DA cells are twice the Lewis rat)
table(merged_samples@meta.data$sample_name) 

merged_samples$strain = ifelse(merged_samples$sample_name %in% c('rat_DA_01_reseq', 'rat_DA_M_10WK_003', 'rat_DA_M09_WK_008'), 'rat_DA', 'rat_LEW')

MIT_PATTERN = '^Mt-'
mito_genes_index <- grep(pattern = MIT_PATTERN, row.names(merged_samples) )
merged_samples[["mito_perc"]] <- PercentageFeatureSet(merged_samples, features = mito_genes_index)


#### checking one of the Lewis genes as an example #### 
test_gene <- Lew_genes_not_in_DA[1]
df <- data.frame(gene_exp_val=GetAssayData(merged_samples, assay='RNA')[test_gene,],
                 sample_name=merged_samples$sample_name)
ggplot(df, aes(y=gene_exp_val, x=sample_name))+geom_violin()  
summary(df$gene_exp_val[df$sample_name=='rat_DA_M09_WK_008']) #### normalized expression is zero in DA
summary(GetAssayData(sample_DA, assay='RNA')[test_gene,]) #### raw expression is low but non-zero in DA
var(GetAssayData(sample_DA, assay='RNA')[test_gene,])  #### variance of the raw expression is low -> removed by the sctransform step


#### Removing all the Mt-genes for later steps
rownames(merged_samples)[mito_genes_index]
merged_samples <- merged_samples[-mito_genes_index,]
grep(pattern = MIT_PATTERN, rownames(merged_samples))

##### Run PCA and Harmony on the data  #####
merged_samples <- RunPCA(merged_samples, features = rownames(merged_samples))  


## the genes that are not expressed in one of the samples, are removed throughout the PCA 
length(DA_genes_not_in_Lew) + length(Lew_genes_not_in_DA)
length(rownames(merged_samples)) - nrow(Loadings(merged_samples, 'pca'))
loading_genes <- rownames(Loadings(merged_samples, 'pca'))
sum(DA_genes_not_in_Lew %in% loading_genes)
sum(Lew_genes_not_in_DA %in% loading_genes)


merged_samples <- RunHarmony(merged_samples, "sample_name",assay.use="RNA")



########   Cluster the merged samples   #########
dir = 'Results/'
merged_samples <- FindNeighbors(merged_samples,
                                reduction="harmony",
                                dims=1:PC_NUMBER,
                                verbose=T)

sCVdata_list <- get_cluster_object(seur = merged_samples, 
                                   max_seurat_resolution = 2.4, ## change this to higher values
                                   FDRthresh = 0.05, # FDR threshold for statistical tests
                                   min_num_DE = 1,
                                   seurat_resolution = 0, # Starting resolution is this plus the jump value below.
                                   seurat_resolution_jump = 0.2,
                                   DE_bw_clust = TRUE)


merged_samples@meta.data <- cbind(merged_samples@meta.data,
                                 data.frame(lapply(sCVdata_list, function(x) x@Clusters)))
head(merged_samples@meta.data)

######################################
######### Run Varimax on the merged samples #########

### needed data for the varimax PCA 
loading_matrix = Loadings(merged_samples, 'pca')
gene_exp_matrix = GetAssayData(merged_samples, assay = 'RNA')
dim(gene_exp_matrix)
dim(loading_matrix)


### Apply the varimax rotation

# rot_data <- readRDS('Results/rotatedPCA_oldSamples.rds')
rot_data <- readRDS('Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds') ## MT-removed
rot_data <- readRDS('Results/new_samples/varimax_rotated_object_new.rds') ## MT-removed

rot_data <- readRDS('Results/new_samples/immune_varimax_results.rds') ## MT-removed

rot_data <- get_varimax_rotated(gene_exp_matrix, loading_matrix)
rotatedLoadings <- rot_data$rotLoadings
scores <- data.frame(rot_data$rotScores)
colnames(scores) = paste0('Varimax_', 1:ncol(scores))
dim(scores)
dim(rotatedLoadings)
dim(merged_samples) ### running PCA removes around 20 genes which have low variance
head(rotatedLoadings)


##### Visualize the rotated PCs and evaluate mapping with clusters ######

### which PCs to evaluate in the rotated data

### PCA
pca_embedding_df <- data.frame(Embeddings(merged_samples,reduction = 'pca'))
PC_standardDev <- apply(pca_embedding_df, 2, sd)
elbow_plot(PC_standardDev, title = 'PCA on last 2 rat samples')
ElbowPlot(merged_samples)
top_pc = 15

### varimax PCA
rot_standardDev <- apply(rot_data$rotScores, 2, sd)
elbow_plot(rot_standardDev, title = 'Varimax on last 2 rat samples')
ndims = 50
rot_percVar = (rot_standardDev^2 * 100)/sum(rot_standardDev^2)
names(rot_percVar) <- paste0(1:(length(rot_percVar)))
perc_variance_threshold = 0.5
PCs_to_check <- as.numeric(names(rot_percVar[rot_percVar>perc_variance_threshold  ]))
PCs_to_check <- 1:31

#### PC and varimax correlation with technical covariates ####
pca_embedding_df$clusters = as.character(merged_samples$res.0.6)

check_qc_cor(merged_samples, pca_embedding_df[,1:top_pc], 
             main='PC embeddings correlation with technical covariates')

embedd_df_rotated <- data.frame(scores)
embedd_df_rotated_2 <- embedd_df_rotated[,PCs_to_check]
colnames(embedd_df_rotated_2) <- paste0('Var.PC_', PCs_to_check)
check_qc_cor(merged_samples, embedd_df_rotated_2, 
             main='Varimax-PC embeddings correlation with technical covariates')

pheatmap(cor(embedd_df_rotated[,c(PCs_to_check)]))
lapply(embedd_df_rotated, head)

##### Visualizing the PCA plots before the rotation #####

plot_dir = 'Plots/'
pdf(paste0(plot_dir, 'PCA_plots_mergedNewSamples_MT-removed.pdf'),width = 14,height=12)
#pdf(paste0(plot_dir, 'PCA_plots_mergedOldSamples_mt40_lib1500_MTremoved.pdf'),width = 14,height=12)


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




##### Visualizing the PCA plots after the varimax rotation #####

pdf(paste0(plot_dir, 'VarimaxPCA_scTrans_mergedNewSamples_scaled_-mt-removed.pdf'),width = 14,height=13) 
#pdf(paste0(plot_dir, 'VarimaxPCA_scTrans_mergedOldSamples_mt40_lib1500_MTremoved_1.pdf'),width = 14,height=13) 
for(i in PCs_to_check){ 
  pc_num = i
  rot_df <- data.frame(Varimax_1=embedd_df_rotated$Varimax_1,
                       emb_val=embedd_df_rotated[,pc_num],
                       cluster=as.character(merged_samples$res.0.6),
                       Library_size=merged_samples$nCount_RNA,
                       num_expressed_genes=merged_samples$nFeature_RNA,
                       strain=merged_samples$strain,
                       sample_name = merged_samples$sample_name)
  
  p1=ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=cluster))+geom_point()+
    theme_classic()+ylab(paste0('Varimax_',i))+scale_color_manual(values=colorPalatte)
  p2=ggplot(rot_df, aes(x=cluster, y=emb_val, fill=cluster))+geom_violin()+
    theme_classic()+scale_fill_manual(values=colorPalatte)
  p3=ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=strain))+geom_point()+
    theme_classic()+ylab(paste0('Varimax_',i))
  p4=ggplot(rot_df, aes(x=strain, y=emb_val, fill=strain))+geom_violin()+theme_classic()
  
  gridExtra::grid.arrange(p1,p2,p3,p4,ncol=2,nrow=2)
}
dev.off()                      



##### checking the sample composition of each cluster ##### 

clusters_pred = as.character(merged_samples$res.0.6)
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
                      clusters= as.character(merged_samples$res.0.6), 
                      sample_name=merged_samples$sample_name, 
                      strain=merged_samples$strain,
                      gene_expression=a_gene_expression,
                      detectedGenes_More250=merged_samples$nFeature_RNA>250
                      )

pdf(paste0(plot_dir, 'umap_plots_oldSamples_Mt_removed.pdf'))
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=library_size))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=n_expressed))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=mito_perc))+geom_point()+theme_classic()+scale_color_viridis(direction=1)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=gene_expression))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)+ggtitle(a_gene_name)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=sample_name))+geom_point(alpha=0.7)+theme_classic()+scale_color_manual(values = colorPalatte)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=strain))+geom_point(alpha=0.6)+theme_classic()+scale_color_manual(values = colorPalatte)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=clusters))+geom_point()+theme_classic()+scale_color_manual(values = colorPalatte)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=detectedGenes_More250))+geom_point(alpha=0.6)+theme_classic()+scale_color_manual(values = colorPalatte)

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
                      clusters=as.character(merged_samples$res.0.6), 
                      sample_name=merged_samples$sample_name,
                      strain=merged_samples$strain,
                      gene_expression=a_gene_expression)

pdf(paste0(plot_dir, 'tsne_plots_oldSamples_Mt_removed.pdf'))
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

for(pc_index in PCs_to_check){ 
  print(pc_index)
  loading_df <- data.frame(genes=rownames(rotatedLoadings),values=rotatedLoadings[,pc_index])
  loading_df_ord <- loading_df[order(loading_df$values, decreasing = T),]
  rownames(loading_df_ord) <- NULL
  print(head(loading_df_ord, 25))
  print(tail(loading_df_ord, 25))
  write.table(loading_df_ord, 
              paste0('~/RatLiver/Results/new_samples/rotated_loadings/rot_PC',pc_index,'_loadings.rnk'), 
              #paste0('~/RatLiver/Results/old_samples/rotated_loadings/rot_PC',pc_index,'_loadings.rnk'), 
            col.names = F, row.names = F, quote = F, sep = '\t')
}


#### running DE between DA and LEW cells in cluster 4 to check the contamination removal #### 

merged_samples_cluster4 <- merged_samples[,merged_samples$res.0.2 == 4]
Idents(merged_samples_cluster4) <- merged_samples_cluster4$strain
table(Idents(merged_samples_cluster4))
strain_gp <- c('rat_DA', 'rat_LEW')
Cluster4_strain_markers <- sapply(1:length(strain_gp), 
                          function(i) FindMarkers(merged_samples_cluster4, ident.1=strain_gp[i]), 
                          simplify = FALSE)
names(Cluster4_strain_markers) <- strain_gp
lapply(Cluster4_strain_markers, head,20)

#######################################################


merged_samples <- readRDS('~/RatLiver/Objects/merged_samples_newSamples.rds')

#saveRDS(all_merged_samples_seur, 'Results/preproc_rats/merged/all_merged_samples_seur.rds')
#merged_samples <- readRDS('Objects/merged_samples_newSamples.rds')


merged_samples <- readRDS('Objects/merged_samples_oldSamples.rds')
merged_samples <- readRDS('Objects/merged_samples_newSamples.rds')

merged_samples <- readRDS('Objects/merged_samples_oldSamples_mt40_lib1500_MTremoved.rds')
merged_samples <- readRDS('Objects/merged_samples_newSamples_MT-removed.rds')

### checking if cluster 10 is still a cluster after removal of it's MT genes
cluster_10_UMIs <- colnames(merged_samples)[as.character(merged_samples$res.1) == '1']
merged_samples_MTremoved$isCluster10 <- colnames(merged_samples_MTremoved) %in% cluster_10_UMIs

x = data.frame(getEmb(merged_samples_MTremoved, 'UMAP'), merged_samples_MTremoved$isCluster10)
head(x)
ggplot(x, aes(x=UMAP_1,y=UMAP_2,color=merged_samples_MTremoved.isCluster10))+geom_point()
table(merged_samples_MTremoved$res.1[merged_samples_MTremoved$isCluster10]) ## mainly is in the cluster 3 



### finding the equivalent clusters for the new version of the new samples
table(paste0('mt-removed-',as.character(merged_samples_mt_removed$res.1.2)),as.character(merged_samples$res.1))
data.frame(mt_removed=c('0', '1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', 
                        '2', '20', '21', '22', '23', '3', '4', '5', '6', '7', '8', '9'),
           mt_included=c('2', '0', '8', '11' , '12', '16 1' , '13 8', '14', '19 10', '15',
                         '18', '17', '1', '20', '21', '22', '23', '3', '4', '5', '10', '7', '9', '6'))




