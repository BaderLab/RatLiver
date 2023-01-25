
## loading packages
source('Codes/Functions.R')
source('Codes/FactorAnalysisUtils.R')
Initialize()

merged_samples = readRDS('~/rat_sham_sn_data/sham_sn_merged_data.rds')
DefaultAssay(merged_samples) <- 'RNA'
merged_samples@assays$SCT <- NULL
merged_samples = SCTransform(merged_samples, conserve.memory=F, verbose=T, return.only.var.genes=F,
            variable.features.n = nrow(merged_samples[['RNA']]), 
            ## according to the paper scaling is not recommended prior to PCA:
            ## https://www.biorxiv.org/content/10.1101/576827v2
            do.scale = FALSE, ### default value 
            do.center = TRUE) 

merged_samples <- ScaleData(merged_samples)

################################################################
#### now sure if the samples need to be normalized individually 
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

####################################################################


### needed data for the varimax PCA 
merged_samples <- RunPCA(merged_samples, features = rownames(merged_samples))  

loading_matrix = Loadings(merged_samples, 'pca')
gene_exp_matrix = GetAssayData(merged_samples, assay = 'RNA') #SCT
dim(gene_exp_matrix)
dim(loading_matrix)

rot_data <- get_varimax_rotated(gene_exp_matrix, loading_matrix)
rotatedLoadings <- rot_data$rotLoadings
scores <- data.frame(rot_data$rotScores)
colnames(scores) = paste0('PC_', 1:ncol(scores))
head(scores)
head(rotatedLoadings)






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


embedd_df_rotated <- data.frame(scores)
PCs_to_check = 1:ncol(embedd_df_rotated)
embedd_df_rotated_2 <- embedd_df_rotated[,PCs_to_check]

colnames(embedd_df_rotated_2) <- paste0('Var.PC_', PCs_to_check)
check_qc_cor(merged_samples, embedd_df_rotated_2, 
             main='Varimax-PC embeddings correlation with technical covariates')


pdf('~/rat_sham_sn_data/VarimaxPCA_sham_sn_merged.pdf',width = 14,height=16) 

for(i in PCs_to_check){ 
  pc_num = i
  rot_df <- data.frame(Varimax_1=embedd_df_rotated_2$Var.PC_1,
                       emb_val=embedd_df_rotated_2[,pc_num],
                       #cluster=as.character(merged_samples$final_cluster),
                       cluster=as.character(merged_samples$cluster),
                       label = as.character(merged_samples$annot_TLH_g),
                       sample_name=merged_samples$sample_name,
                       strain=merged_samples$strain,
                       nuclear_fraction = merged_samples$nuclear_fraction,
                       nCount_RNA = merged_samples$nCount_RNA)
  
  p1=ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=cluster))+geom_point()+
    theme_classic()+ylab(paste0('Varimax_',i))+scale_color_manual(values=colorPalatte)
  p2=ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=label))+geom_point()+
    theme_classic()+ylab(paste0('Varimax_',i))+scale_color_manual(values=colorPalatte)
  p3=ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=sample_name))+geom_point()+
    theme_classic()+ylab(paste0('Varimax_',i))+scale_color_manual(values=colorPalatte)
  p6=ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=strain))+geom_point()+
    theme_classic()+ylab(paste0('Varimax_',i))+scale_color_manual(values=colorPalatte)
  
  p4=ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=nuclear_fraction))+geom_point()+
    theme_classic()+ylab(paste0('Varimax_',i))+scale_color_viridis(direction = -1)
  p5=ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=nCount_RNA))+geom_point()+
    theme_classic()+ylab(paste0('Varimax_',i))+scale_color_viridis(direction = -1)
  
  
  #p5=ggplot(rot_df, aes(x=sample_name, y=emb_val, color=nCount_RNA))+geom_boxplot()+theme_classic()
  gridExtra::grid.arrange(p1,p2,p3,p6,p4,p5,ncol=2,nrow=3)
}
dev.off()                      

