### loading the required libraries
source('Codes/Functions.R')
source('Codes/FactorAnalysisUtils.R')
Initialize()

get_scaled_by_gene <- function(x){
  y <- scale(t(x), center = TRUE, scale = TRUE)
  y[is.nan(y)] <- 0
  y2 <- t(y)
  y2 <- as.data.frame(y2)
  return(y2)
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

qc_info = list(MIT_CUT_OFF=40, LIB_SIZE_CUT_OFF=2000)
sample_names = c('rat_DA_M_10WK_003', 'rat_DA_01_reseq', 'rat_Lew_01', 'rat_Lew_02' )
strain_info = c('DA', 'DA', 'LEW', 'LEW')

data_norm_list = list(NA, NA, NA, NA)
names(data_norm_list) = sample_names


for(i in 1:length(sample_names)){
  sample_name = sample_names[i]
  
  input_from_10x <- paste0("Data/", sample_name,'/')
  data <- CreateSeuratObject(counts=Read10X(input_from_10x, gene.column = 2),
                                 min.cells=0,min.features=1, 
                                 project = "snRNAseq")
  dim(data)
  
  
  MIT_PATTERN = '^Mt-'
  mito_genes_index <- grep(pattern = MIT_PATTERN, rownames(data) )
  data[["mito_perc"]] <- PercentageFeatureSet(data, features = mito_genes_index)
  
  to_drop_mito = data$mito_perc > qc_info$MIT_CUT_OFF
  to_drop_lib_size = data$nCount_RNA < qc_info$LIB_SIZE_CUT_OFF
  to_drop_numgenes = data$nFeature_RNA < 100 
  data$sample_name = sample_name
  
  summary(data$mito_perc )
  summary(data$nFeature_RNA )
  sum(!to_drop_mito & !to_drop_numgenes)
  sum(!to_drop_mito & !to_drop_lib_size)
  sum(to_drop_lib_size)
  sum(to_drop_numgenes)
  sum(to_drop_mito)
  
  data <- data[,!to_drop_mito & !to_drop_lib_size]
  data = SCTransform(data, vst.flavor = "v2", verbose = TRUE)
  
  data_norm_list[[i]] = data
}

rm(data)
gc()

merged_data <- merge(data_norm_list[[1]], c(data_norm_list[[2]], data_norm_list[[3]], data_norm_list[[4]] ), # 
                     add.cell.ids = names(data_norm_list), 
                     project = names(data_norm_list), 
                     merge.data = TRUE)

dim(merged_data)
Idents(merged_data) = merged_data$sample_name


merged_data <- SCTransform(merged_data,verbose=T)
merged_data = FindVariableFeatures(merged_data)

merged_data <- RunPCA(merged_data,verbose=T)
plot(100 * merged_data@reductions$pca@stdev^2 / merged_data@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")

#### UMAP on corrected PC components and clustering
merged_data <- RunHarmony(merged_data, group.by.vars = "sample_name")
merged_data <- RunUMAP(merged_data, reduction = "harmony", dims = 1:30, reduction.name = "umap_h")  %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

Resolution = 0.4
merged_data <- FindClusters(merged_data, resolution = Resolution, verbose = FALSE)
gene_name = 'Ecm1'
cluster_num = 5
is_cluster = ifelse(merged_data$seurat_clusters== cluster_num, paste0('cluster',cluster_num), 'other')

merged_data <- readRDS('Results/review_analysis/merged_samples_TLH_harmonized_qc_mt40_lib2000.rds')

df_umap <- data.frame(#UMAP_1=getEmb(merged_data, 'umap')[,1], 
  #UMAP_2=getEmb(merged_data, 'umap')[,2], 
  UMAP_h_1=getEmb(merged_data, 'umap_h')[,1], 
  UMAP_h_2=getEmb(merged_data, 'umap_h')[,2], 
  library_size= merged_data$nCount_RNA, 
  mito_perc=merged_data$mito_perc, 
  n_expressed=merged_data$nFeature_RNA,
  cluster=merged_data$seurat_clusters, 
  is_cluster=is_cluster,
  #cell_status = merged_data$cell_status,
  #nuclear_fraction=merged_data$nuclear_fraction, 
  Alb=GetAssayData(merged_data)['Alb',], 
  gene=GetAssayData(merged_data)[gene_name,],
  #strain = merged_data$strain
  sample_name = merged_data$sample_name)

ggplot(df_umap, aes(x=UMAP_h_1, y=UMAP_h_2, color=gene))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(gene_name)
ggplot(df_umap, aes(x=UMAP_h_1, y=UMAP_h_2, color=is_cluster))+geom_point(alpha=0.8)+theme_classic()#+ggtitle(paste0('res: ', Resolution))
ggplot(df_umap, aes(x=UMAP_h_1, y=UMAP_h_2, color=cluster))+geom_point(alpha=0.8)+theme_classic()+
  scale_color_manual(values = colorPalatte)#+ggtitle(paste0('res: ', Resolution))
ggplot(df_umap, aes(x=UMAP_h_1, y=UMAP_h_2, color=sample_name))+geom_point(alpha=0.3)+theme_classic()#+ggtitle(paste0('res: ', Resolution))

ggplot(df_umap, aes(x=UMAP_h_1, y=UMAP_h_2, color=Alb))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle('Alb')

########################################################################################################################



#######################################################################################################
##############  Running Varimax on the harmonized QC samples   ########################################
######################################################################################################

files_rds = sapply(1:length(data_norm_list), 
                   function(i) {
                     x = data_norm_list[[i]]
                     # x = x[,!x$SCT_snn_res.0.7 %in% hep_clusters[[i]]]
                     
                     DefaultAssay(x) <- 'RNA'
                     x@assays$SCT <- NULL
                     
                     x = SCTransform(x, conserve.memory=F, verbose=T, return.only.var.genes=F,
                                     variable.features.n = nrow(x[['RNA']]), 
                                     ## according to the paper scaling is not recommended prior to PCA:
                                     ## https://www.biorxiv.org/content/10.1101/576827v2
                                     do.scale = FALSE, ### default value 
                                     do.center = TRUE) ### default value ))
                     
                     x = CreateSeuratObject(GetAssayData(x, assay='SCT'))
                     x$sample_name = sample_names[i]
                     x$strain = strain_info[i]
                     
                     return(x)
                   },simplify = F)

names(files_rds) = sample_names


merged_samples <- merge(files_rds[[1]], c(files_rds[[2]], files_rds[[3]], files_rds[[4]]), 
                        add.cell.ids = names(files_rds), 
                        project = "harmonizedQC", 
                        merge.data = TRUE)

dim(merged_samples)
lapply(files_rds, dim)


merged_samples = SCTransform(merged_samples, conserve.memory=F, verbose=T, return.only.var.genes=F,
                variable.features.n = nrow(merged_samples[['RNA']]), 
                ## according to the paper scaling is not recommended prior to PCA:
                ## https://www.biorxiv.org/content/10.1101/576827v2
                do.scale = FALSE, ### default value 
                do.center = TRUE) ### default value ))

#merged_samples <- ScaleData(merged_samples, features = rownames(merged_samples))  
############ adding metadata calculated using the standard analysis
######## using the complete data
sum(colnames(merged_samples) != colnames(merged_data))
rownames(merged_samples@meta.data) == rownames(merged_data@meta.data)
merged_samples@meta.data = merged_data@meta.data

merged_samples$strain=unlist(lapply(str_split(string = merged_data@meta.data$sample_name,pattern = '_'), '[[',2))


####################################################################
##################    Peforming the varimax analysis ################
####################################################################
### needed data for the varimax PCA 
merged_samples <- RunPCA(merged_samples, features = rownames(merged_samples))  

loading_matrix = Loadings(merged_samples, 'pca')
# USING RNA is REALLY IMPORTANT. REMOVING THE SCT does not work and seurat obj seems to remember that 
gene_exp_matrix = GetAssayData(merged_samples, assay = 'RNA') 
dim(gene_exp_matrix)
dim(loading_matrix)

rot_data <- get_varimax_rotated(gene_exp_matrix, loading_matrix)
rotatedLoadings <- rot_data$rotLoadings
scores <- data.frame(rot_data$rotScores)
colnames(scores) = paste0('PC_', 1:ncol(scores))
head(scores)
head(rotatedLoadings)
plot(scores[,1], scores[,2]) ### checking if varimax factors look meaningful

embedd_df_rotated <- data.frame(scores)
PCs_to_check = 1:ncol(embedd_df_rotated)
embedd_df_rotated_2 <- embedd_df_rotated[,PCs_to_check]
colnames(embedd_df_rotated_2) <- paste0('Var_', PCs_to_check)


rot_data <- readRDS('Results/review_analysis/VarimaxPCA_merged_samples_allfeatures_TLH_harmonized_qc_mt40_lib2000_scTrans_res.rds')
#### adding meta data to the varimax data frame
sum(rownames(merged_samples@meta.data) != rownames(embedd_df_rotated_2))
embedd_df_rotated_2 = data.frame(cbind(embedd_df_rotated_2, merged_samples@meta.data)) ### run for the complete data

check_qc_cor(merged_samples, embedd_df_rotated_2[,PCs_to_check], 
             main='Varimax-PC embeddings correlation with technical covariates')


pdf('Results/review_analysis/VarimaxPCA_merged_samples_allfeatures_TLH_harmonized_qc_mt40_lib2000_scTrans_scale.pdf',width = 10,height=13) 

for(i in PCs_to_check){ 
  pc_num = i
  rot_df <- data.frame(Varimax_1=embedd_df_rotated_2$Var_1,
                       emb_val=embedd_df_rotated_2[,pc_num],
                       #cluster=as.character(merged_samples$final_cluster),
                       cluster=as.character(embedd_df_rotated_2$seurat_clusters),
                       #label = as.character(embedd_df_rotated_2$annot_TLH),
                       sample_name=embedd_df_rotated_2$sample_name,
                       strain=embedd_df_rotated_2$strain,
                       #nuclear_fraction = embedd_df_rotated_2$nuclear_fraction,
                       nCount_RNA = embedd_df_rotated_2$nCount_SCT,
                       UMAP_h_1=getEmb(merged_data, 'umap_h')[,1], 
                       UMAP_h_2=getEmb(merged_data, 'umap_h')[,2])
  
  p1=ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=cluster))+geom_point(alpha=0.8, size=1.1)+
    theme_classic()+ylab(paste0('Varimax_',i))+scale_color_manual(values=colorPalatte)

  p2=ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=sample_name))+geom_point(alpha=0.4, size=1.2)+
    theme_classic()+ylab(paste0('Varimax_',i))+scale_color_manual(values=colorPalatte)
  p3=ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=strain))+geom_point(alpha=0.4, size=1.2)+
    theme_classic()+ylab(paste0('Varimax_',i))+scale_color_manual(values=colorPalatte)
  
  p4=ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=nCount_RNA))+geom_point(size=1.2)+
    theme_classic()+ylab(paste0('Varimax_',i))+scale_color_viridis(direction = -1)
  
  
  #p5=ggplot(rot_df, aes(x=sample_name, y=emb_val, color=nCount_RNA))+geom_boxplot()+theme_classic()
  gridExtra::grid.arrange(p1,p2,p3,p4,ncol=2,nrow=2)
}
dev.off()  

strain_factor = 23
strain_factor = 7
i = strain_factor
head(loading_matrix[order(loading_matrix[,strain_factor], decreasing = T),],20)
head(loading_matrix[order(loading_matrix[,strain_factor], decreasing = F),],20)


ggplot(rot_df, aes(x=UMAP_h_1, y=UMAP_h_2, color=emb_val))+geom_point(alpha=0.4,size=1)+theme_classic()+
  scale_color_viridis(direction = -1, option = 'magma')+ggtitle(paste0('Varimax ', strain_factor))

ggplot(rot_df, aes(y=emb_val, x=strain, fill=strain))+geom_boxplot()+theme_classic()+
  ggtitle(paste0('Varimax ', strain_factor, ''))+ylab(paste0('Varimax-', strain_factor, ' Score'))

ggplot(rot_df[rot_df$cluster==5,], aes(y=emb_val, x=strain, fill=strain))+geom_boxplot()+theme_classic()+
  ggtitle(paste0('Varimax ', strain_factor, ' Macrophages'))+ylab(paste0('Varimax-', strain_factor, ' Score'))


######################################################
varimax_res_harmonized <- readRDS('Results/review_analysis/VarimaxPCA_merged_samples_allfeatures_TLH_harmonized_qc_mt40_lib2000_scTrans_res.rds')
varimax_res_TLH <- readRDS('Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds')

rotLoadings_harmonized = data.frame(sapply(1:ncol(varimax_res_harmonized$rotLoadings),function(i)varimax_res_harmonized$rotLoadings[,i],simplify = T ))
rotLoadings_harmonized = rotLoadings_harmonized[,1:25]
rotLoadings_TLH = data.frame(sapply(1:ncol(varimax_res_TLH$rotLoadings),function(i)varimax_res_TLH$rotLoadings[,i],simplify = T ))
rotLoadings_TLH = rotLoadings_TLH[,1:25]

colnames(rotLoadings_harmonized) = paste0('var_',1:ncol(rotLoadings_harmonized), '.hqc')
colnames(rotLoadings_TLH) = paste0('var_',1:ncol(rotLoadings_TLH), '.v0')

rotLoadings_harmonized$genes = rownames(rotLoadings_harmonized)
rotLoadings_TLH$genes = rownames(rotLoadings_TLH)


merged_var_loadings = merge(rotLoadings_harmonized, rotLoadings_TLH, by.x='genes', by.y='genes')
head(merged_var_loadings)

cor_mat = cor(merged_var_loadings[,-1])
colnames(cor_mat)
cor_mat_sub = cor_mat[1:(ncol(rotLoadings_harmonized)-1), ncol(rotLoadings_harmonized):ncol(cor_mat)]
pheatmap(cor_mat_sub)





########################################################################################################
##################     Performing correaltion between average expression of clusters  ##################
########################################################################################################
##################################################################
############## calculating the cluster average expression of the new data ############## 

nb.cols <- length(names(table(merged_data$seurat_clusters)))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols) #Pastel1

merged_data <- FindNeighbors(merged_data, reduction = "harmony", dims = 1:30)
merged_data <- FindClusters(merged_data, resolution = 0.7)
df_umap$cluster = merged_data$SCT_snn_res.0.7

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cluster))+geom_point(size=1,alpha=0.6)+
  theme_classic()+scale_color_manual(values = c(colorPalatte)) #mycolors
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cluster))+geom_point(size=1,alpha=0.6)+theme_classic()


merged_data$cluster = as.character(merged_data$seurat_clusters)
cluster_names_types = names(table(merged_data$cluster))
### dividing the expression matrix based on the clusters
cluster_expression <- sapply(1:length(cluster_names_types), function(i){
  a_cluster_name = cluster_names_types[i]
  GetAssayData(merged_data, 'data')[,merged_data$cluster == a_cluster_name] 
}, simplify = F)

names(cluster_expression) = cluster_names_types
lapply(cluster_expression, head)

## calculate the average expression of each gene in each cluster
cluster_average_exp <- lapply(cluster_expression, function(x){
  df = data.frame(average=rowSums(x)/ncol(x))
  return(df)
})
lapply(cluster_average_exp, dim)
lapply(cluster_average_exp, head)
## Concatenate all the clusters together to make a matrix
cluster_average_exp_df = do.call(cbind,cluster_average_exp)
#colnames(cluster_average_exp_df) = paste0('cluster_',names(cluster_average_exp))
colnames(cluster_average_exp_df) = paste0('c_',names(cluster_average_exp))
head(cluster_average_exp_df)

## scale and center all the genes in the matrix
cluster_average_exp_df <- get_scaled_by_gene(cluster_average_exp_df) ## scaling over the clusters of interest
cluster_average_exp_df$rat_ID=rownames(cluster_average_exp_df)
head(cluster_average_exp_df)

cluster_average_exp_df = cluster_average_exp_df[rownames(cluster_average_exp_df)%in%VariableFeatures(merged_samples),]
head(cluster_average_exp_df)
dim(cluster_average_exp_df)


##################################################################
###### importing cluster average means of reference maps 
ref_cluster_average_exp_df <- readRDS('~/RatLiver/Results/rat_old_cluster_average_exp_all.rds')

### refine the annotation names for each of the clusters
colnames_set1 = c("Hep (0)", "Hep (1)", "Hep (12)", "Hep (15)", "Hep (16)", "Hep (2)", "Hep (4)", 
                  "Hep (6)", "Hep (8)", "Lyz2/Cd74 Mo/Mac (9)", "Endothelial (11)", "Endothelial (3)", "Lymphocyte (13)",  
                  "Marco/Cd5l Mac (10)",  "Marco/Cd5l Mac (5)", "Mesenchymal (14)", "Mesenchymal (7)", "rat_ID" )

colnames(ref_cluster_average_exp_df) = colnames_set1
nFeatures = 2000

old_data_scClustViz_object <- "~/RatLiver/Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)

ref_data = your_scRNAseq_data_object
DefaultAssay(ref_data) <- 'RNA'
ref_data <- FindVariableFeatures(ref_data, nfeatures=nFeatures)
HVGs <- VariableFeatures(ref_data)

ref_cluster_average_exp_df = ref_cluster_average_exp_df[rownames(ref_cluster_average_exp_df)%in%HVGs,]


########################################################################
######## merging the Healthy Rat ref data with the new samples ######## 
average_exp_merged = merge(cluster_average_exp_df, ref_cluster_average_exp_df, by.x='rat_ID', by.y='rat_ID')

dim(average_exp_merged) ## additional layers to ne to be added to the filter -> 10126 one-2-one orthos
head(average_exp_merged)
rat_cor_mat = cor(average_exp_merged[,-1],method = 'pearson')
colnames(rat_cor_mat)

colnames(rat_cor_mat) <- gsub('Inflammatory', 'Inf', colnames(rat_cor_mat))
rownames(rat_cor_mat) = colnames(rat_cor_mat)
number_of_clusters = length(cluster_names_types)
rat_cor_mat = rat_cor_mat[1:number_of_clusters,(number_of_clusters+1):ncol(rat_cor_mat)] 
pheatmap(rat_cor_mat,color = inferno(20),  clustering_method='ward.D2')


get_matched_label <- function(index, cor.mat, thr){
  label = 'unknown'
  a_clust_value = cor.mat[index,]
  tmp_label =  names(which.max(cor.mat[index,]))
  if(a_clust_value[tmp_label]>thr) label = tmp_label
  return(label)
}

annotation.df = data.frame(cluster=rownames(rat_cor_mat),  
                           annotation=sapply(1:nrow(rat_cor_mat), 
                                             function(i) get_matched_label(index=i, rat_cor_mat, thr=0.3) ))

annotation.df$annotation_g = sub(" \\(.*", "", annotation.df$annotation)

dev.off()
gridExtra::grid.table(annotation.df[,-3])
dev.off()





