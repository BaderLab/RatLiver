### ToDo: try dynamic tree cutting on teh dendogram to define the clusters
source('Codes/Functions.R')
Initialize()
library(Seurat)
library(BiocParallel)
library(CoGAPS)


elbow_plot <- function(standardDev, title=''){
  ## Extract the standard deviation value for each Principal component and plot them
  ## INPUTS:
  ##  seurat_obj: seurat object
  ## RETURNS:
  ##  elbow (scree) plot for the dimension reduction
  
  ndims = 50
  percVar = (standardDev^2 * 100)/sum(standardDev^2)
  plot <- ggplot(data = data.frame(dims = 1:ndims, percVar = percVar[1:ndims])) +
    geom_point(mapping = aes_string(x = 'dims', y = 'percVar')) +
    labs( x= 'PCs', y = 'Percentage of explained variance')+ggtitle(title)
  theme_classic()
  print(plot)
}


new_data_scCLustViz_object_Immune <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_ImmuneSub.RData"
load(new_data_scCLustViz_object_Immune)

##### subseting the immune population of the new rat samples
### each complete sample has been normalized individually and subsetted for the immune clusters and merged and scaled 

merged_samples_sub <- readRDS('Results/new_samples/Immune_subclusters.rds') 
merged_samples_sub_raw <- readRDS('merged_new_samples_raw_immune.rds')
data <- as.matrix(GetAssayData(merged_samples_sub_raw))

######## PCA results on the immune clusters ########
PCA_embedding <- data.frame(Embeddings(your_scRNAseq_data_object, 'pca'))
head(PCA_embedding)

standardDev <- apply(PCA_embedding, 2, sd)
elbow_plot(standardDev, title = 'merged rat samples')
PCA_CUT_OFF = 10 # 11

PCA_embedding <- PCA_embedding[,1:PCA_CUT_OFF]
PCA_embedding$umi <- rownames(PCA_embedding)
head(PCA_embedding)



######## Varimax results on the immune clusters ########
rot_data <- readRDS('Results/new_samples/immune_varimax_results.rds') ### Varimax results on the separately scaled dataset

varimax_loadings <- data.frame(sapply(1:ncol(rot_data$rotLoadings), function(i) rot_data$rotLoadings[,i]))
colnames(varimax_loadings) <- paste0('Varimax_', 1:ncol(varimax_loadings))


varimax_embeddings <- data.frame(rot_data$rotScores)
colnames(varimax_embeddings) <- paste0('Varimax_', 1:ncol(varimax_embeddings))


rot_standardDev <- apply(varimax_embeddings, 2, sd)
elbow_plot(rot_standardDev, title = 'merged rat samples')
rot_percVar = (rot_standardDev^2 * 100)/sum(rot_standardDev^2)
names(rot_percVar) <- colnames(varimax_embeddings)
PCs_to_check <- names(rot_percVar[rot_percVar>quantile(rot_percVar, 0.5)])

varimax_embeddings <- varimax_embeddings[,colnames(varimax_embeddings) %in% PCs_to_check]
varimax_embeddings$umi <- rownames(varimax_embeddings)
head(varimax_embeddings)



######## scCoGAPs results ########
params <- new("CogapsParams")
params <- setDistributedParams(params, nSets=detectCores()-2)
getParam(params, "nPatterns")

# CoGAPS requires data to be genes x samples
Results.sc <- CoGAPS(data, params, distributed="single-cell", 
                     messages=T, transposeData=FALSE) #nIterations=1000
Results.sc <- readRDS('Results/new_samples/scCoGAPS_immuneCells.rds')

scCoGAPs_embeddings <- data.frame(Results.sc@sampleFactors)
scCoGAPs_embeddings$umi <- rownames(scCoGAPs_embeddings)
colnames(scCoGAPs_embeddings)[1:5] <- paste0('scCoGAPs_', 1:5)
head(scCoGAPs_embeddings)
head(Results.sc@featureLoadings)



######## cNMF results ########
cNMF_loading <- read.delim('Results/new_samples/cNMF/immune_ratLiver/cNMF/immune_ratLiver/immune_ratLiver.gene_spectra_score.k_7.dt_0_18.txt')
cNMF_loading.t <- t(cNMF_loading)
cNMF_loading <- data.frame(cNMF_loading.t[-1,])
colnames(cNMF_loading) <- paste0('cNMF', 1:ncol(cNMF_loading))
head(cNMF_loading)
read.delim('Results/new_samples/cNMF/immune_ratLiver/cNMF/immune_ratLiver/immune_ratLiver.gene_spectra_tpm.k_7.dt_0_18.txt')
x  = t(read.delim('Results/new_samples/cNMF/immune_ratLiver/cNMF/immune_ratLiver/immune_ratLiver.spectra.k_7.dt_0_18.consensus.txt'))
cNMF_embedding <- read.delim('Results/new_samples/cNMF/immune_ratLiver/cNMF/immune_ratLiver/immune_ratLiver.usages.k_7.dt_0_18.consensus.txt')
colnames(cNMF_embedding) <- c('umi', paste0('cNMF_', 1:(ncol(cNMF_embedding)-1)))
head(cNMF_embedding)



######## LDVAE results ########
LDVAE_loading <- read.csv('Results/new_samples/LDVAE/LDVAE_loadings.csv')
rownames(LDVAE_loading) <- LDVAE_loading$X
LDVAE_loading <- LDVAE_loading[,-1]
colnames(LDVAE_loading) <- paste0('LDVAE_', 0:(ncol(LDVAE_loading)-1))
head(LDVAE_loading)

LDVAE_embedding <- read.csv('Results/new_samples/LDVAE/LDVAE_embeddings.csv')
LDVAE_embedding <- LDVAE_embedding[,-1]

colnames(LDVAE_embedding) <- paste0('LDVAE_', 0:(ncol(LDVAE_embedding)-1))
rownames(LDVAE_embedding) <- colnames(merged_samples_sub)
LDVAE_embedding$umi <- rownames(LDVAE_embedding)
head(LDVAE_embedding)



######## fscLVM results ########
data <- readRDS('~/RatLiver/Results/new_samples/fsclvm_model_immuneCells.rds')

### accessing embedding and loading matrix
fscLVM_loading <- data.frame(rowData(data))
fscLVM_embedding <- data.frame(reducedDim(data, 'slalom'))
fscLVM_embedding_orig <- fscLVM_embedding
colnames(fscLVM_embedding) <- paste0('fscLVM_', 1:ncol(fscLVM_embedding))
fscLVM_embedding$umi <- colnames(merged_samples_sub)

colnames(fscLVM_embedding)

#fscLVM_embedding$cluster= merged_samples_sub$immune_clusters
#fscLVM_embedding$sample= merged_samples_sub$sample_name


######## NIFA results ########
#### default number of hidden factors are 6 

NIFAres <- readRDS('Results/new_samples/NIFA_immuneCells.rds')
NIFA_embedding <- data.frame(t(NIFAres$mu_S)) ## score matrix 
NIFA_loading <- data.frame(t(NIFAres$mu_A)) ## loading matrix

colnames(NIFA_embedding) <- paste0('NIFA_', 1:ncol(NIFA_embedding))
colnames(NIFA_loading) <- paste0('NIFA_', 1:ncol(NIFA_loading))
NIFA_embedding$umi <- colnames(merged_samples_sub)

head(NIFA_embedding)
head(NIFA_loading)
dim(merged_samples_sub)


################## making the umap  
sum(colnames(merged_samples_sub) != rownames(scCoGAPs_embeddings))
umap_df <- data.frame(Embeddings(your_scRNAseq_data_object, 'umap'))
umap_df$umi <- rownames(umap_df)
umap_df <- merge(umap_df, scCoGAPs_embeddings, 'umi', 'umi', all.x=T)
umap_df <- merge(umap_df, varimax_embeddings, 'umi', 'umi', all.x=T)
umap_df <- merge(umap_df, cNMF_embedding, 'umi', 'umi', all.x=T)
umap_df <- merge(umap_df, LDVAE_embedding, 'umi', 'umi', all.x=T)
#umap_df <- merge(umap_df, PCA_embedding, 'umi', 'umi', all.x=T)
umap_df <- merge(umap_df, fscLVM_embedding, 'umi', 'umi', all.x=T)
umap_df <- merge(umap_df, NIFA_embedding, 'umi', 'umi', all.x=T)

umap_df$strain = unlist(lapply(str_split(umap_df$umi,pattern = '_' ), function(x) x[[2]]))
umap_df$immune_clusters = merged_samples_sub$immune_clusters
umap_df$init_clusters = merged_samples_sub$init_clusters
head(umap_df)

ggplot(umap_df, aes(UMAP_1,UMAP_2, color=LDVAE_0 ))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)
ggplot(umap_df, aes(x=strain, y=LDVAE_0, fill=strain))+geom_boxplot()+theme_bw()

ggplot(umap_df, aes(x=LDVAE_0, fill=strain))+geom_density(alpha=0.4)+theme_bw() #scCoGAPs
ggplot(umap_df, aes(x=Varimax_4, fill=strain))+geom_density(alpha=0.4)+theme_bw()
ggplot(umap_df, aes(x=cNMF_1, fill=strain))+geom_density(alpha=0.4)+theme_bw()
ggplot(umap_df, aes(x=LDVAE_7, fill=strain))+geom_density(alpha=0.4)+theme_bw()


total_factors_df <- umap_df[,4:(ncol(umap_df)-3)]
head(total_factors_df)
pheatmap(cor(scale(total_factors_df)))
total_factors_df <- readRDS('Results/total_factors_df.rds')

pdf('Plots/compareAllFactors_immuneNewSamples.pdf', height = 20, width = 20)
pheatmap(cor(total_factors_df,method = 'pearson'), main = 'pearson')
pheatmap(cor(total_factors_df,method = 'spearman'), main = 'spearman')
dev.off()




##### interesting clusters to investigate ####
gp_1 <- c('LDVAE_0', 'LDVAE_8', 'cNMF_2', 'scCoGAPs_3', 'NIFA_5')
gp_2 <- c('fscLVM_19', 'fscLVM_17', 'fscLVM_34', 'fscLVM_18', 'Varimax_1', 'fscLVM_24' )
total_factors_df.cor <- cor(total_factors_df[,colnames(total_factors_df) %in% c(gp_1, gp_2)],method = 'pearson')
out <- pheatmap(total_factors_df.cor, inferno(10,direction = -1))

### checking varimax_1 form gp-2
varimax_loadings.ord <- varimax_loadings[order(as.numeric(varimax_loadings$Varimax_1), decreasing = T),]
head(varimax_loadings.ord,10) ### inflammatory responses, strain-specific
rownames(head(varimax_loadings.ord,20))

### checking LDVAE_0 form gp-1
LDVAE_loading.ord <- LDVAE_loading[order(LDVAE_loading$LDVAE_0,decreasing = T),]
head(LDVAE_loading.ord,20) ### non-Inf?
tail(LDVAE_loading.ord)
rownames(head(LDVAE_loading.ord,20))

dendogram_thr = 2
total_factors_df.cor <- cor(total_factors_df,method = 'pearson')
out <- pheatmap(total_factors_df.cor, main = 'pearson')
#Re-order original data to match ordering in heatmap (top-to-bottom)
sorted_tree.df = data.frame(sorted_factor=rownames(total_factors_df.cor[out$tree_row[["order"]],]))

plot(out$tree_row)
abline(h=dendogram_thr, col="red", lty=2, lwd=2)


# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/WORKSHOP/2014/Langfelder-NetworkDay-clustering.pdf
library(dynamicTreeCut)
help("cutreeDynamic")
# Input:
# clustering tree
# dissimilarity matrix that was used to produce the tree
# multiple options to fine-tune cluster criteria and PAM stage
# Most important options:
# DeepSplit (0-4): controls how finely clusters will be split
# pamStage (FALSE or TRUE): turns PAM stage off/on

cutreeDynamicTree(out$tree_row, maxTreeHeight=2, deepSplit=F)
cutreeHybrid(out$tree_row, total_factors_df.cor,deepSplit = 3)
#Cut the row  dendrogram at a Euclidean distance dis-similarity of 8
dendrogram.cut.df <- data.frame(clusters=sort(cutree(out$tree_row, h=dendogram_thr)))
dendrogram.cut.df$factors = rownames(dendrogram.cut.df)
dendrogram.cut.df_sorted <- dendrogram.cut.df[match(sorted_tree.df$sorted_factor,dendrogram.cut.df$factors),]

table(dendrogram.cut.df_sorted$clusters[dendrogram.cut.df_sorted$factors %in% strain_specific_factors_names])
table(dendrogram.cut.df_sorted$clusters[dendrogram.cut.df_sorted$factors %in% c('scCoGAPs_2', 'Varimax_4', 'cNMF_4', 'LDVAE_2',
                                                                                'fscLVM_2', 'fscLVM_12', 'fscLVM_14', 'NIFA_6')])

head(dendrogram.cut.df$factors[1])
cluster_df <- data.frame(table(dendrogram.cut.df$clusters))
colnames(cluster_df) = c('cluster', 'num_elements')
head(cluster_df)



#### visualizing factors in each cluster in the UMAP ####

a_cluster = 25
pdf('Plots/test.pdf', width = 15, height = 14)
for(a_cluster in as.numeric(cluster_df$cluster)){
  if(a_cluster == 16) next
  
  dendrogram.cut.df <- data.frame(clusters=sort(cutree(out$tree_row, h=dendogram_thr)))
  dendrogram.cut.df$factors = rownames(dendrogram.cut.df)
  
  dendrogram.cut.df = dendrogram.cut.df[dendrogram.cut.df$clusters==a_cluster,]
  my_plots = sapply(1:nrow(dendrogram.cut.df), function(i){
    a_factor <- dendrogram.cut.df$factors[i]
    a_cluster <- dendrogram.cut.df$clusters[i]
    df = umap_df[,c('UMAP_1', 'UMAP_2', a_factor)]
    colnames(df) = c('UMAP_1', 'UMAP_2', 'Factor')
    p=ggplot(df,aes(UMAP_1,UMAP_2, color=Factor))+geom_point()+theme_classic()+
      scale_color_viridis(option = 'plasma', direction = -1)+
      ggtitle(paste0(a_factor))
    return(p)
  }, simplify = F)
  
  num_cluster_element = cluster_df$num_elements[cluster_df$cluster==a_cluster]
  n_row_val = 1
  n_col_val = 1
  if(num_cluster_element==2) {n_row_val=1; n_col_val=2}
  if(num_cluster_element==3) {n_row_val=1; n_col_val=3}
  if(num_cluster_element==3) {n_row_val=1; n_col_val=3}
  if(num_cluster_element==4) {n_row_val=2; n_col_val=2}
  if(num_cluster_element %in% c(5,6)) {n_row_val=2; n_col_val=3}
  
  p=plot_grid(plotlist = my_plots, nrow = n_row_val, ncol = n_col_val)
  print(p)
}

dev.off()



### mainly in cluster 2, a few in 8 and 24
strain_specific_factors_names_2 <- dendrogram.cut.df_sorted$factors[dendrogram.cut.df_sorted$clusters == 2]

data.frame(strain_specific_factors_names_2)
total_factors_df.pca <- prcomp(total_factors_df, scale = TRUE, center = FALSE)
total_factors_df.pca <- data.frame(total_factors_df.pca$rotation)
total_factors_df.pca$method <- rownames(total_factors_df.pca)
ggplot(total_factors_df.pca, aes(PC1, PC2, color=method))+geom_point()




############ ruuning wilcoxon test on the factors to identify the strain-specific ones #######

is.DA <- umap_df$strain=="DA"
head(total_factors_df)

p_value_list <- rep(-1, ncol(total_factors_df))
FC_list = rep(NA, ncol(total_factors_df))

for(factor_i in 1:ncol(total_factors_df)){
  wilcoxon.res = wilcox.test(x = total_factors_df[is.DA,factor_i], 
                             y = total_factors_df[!is.DA,factor_i],
                             alternative = "two.sided",
                             mu = 0, paired = FALSE)
  
  p_value_list[factor_i] <- wilcoxon.res$p.value
  FC_list[factor_i] <- mean(total_factors_df[is.DA,factor_i])/mean(total_factors_df[!is.DA,factor_i]) # DA/LEW
}
is_significant <- p_value_list < 0.0001/ncol(total_factors_df)
strain_specific_factors <- data.frame(factor=colnames(total_factors_df)[is_significant], 
                                      FC_list = FC_list[is_significant],
                                      p_val=p_value_list[is_significant])
head(strain_specific_factors)
dim(strain_specific_factors)
strain_specific_factors_ord <- strain_specific_factors[order(abs(strain_specific_factors$FC_list), decreasing = T),]
strain_specific_factors_names <- strain_specific_factors_ord$factor


columns_to_check <- which(colnames(umap_df) %in% strain_specific_factors_names)
columns_to_check <- which(colnames(umap_df) %in% strain_specific_factors_names_2)
pdf('~/RatLiver/Plots/strainSpecificFactors_newSamples_immuneCells_2.pdf',width = 21, height = 14)
for(i in columns_to_check){
  
  p1=ggplot(umap_df, aes(y=umap_df[,i], x=strain,fill=strain))+geom_violin()+theme_classic()+ylab(colnames(umap_df)[i])
  p2=ggplot(umap_df, aes(umap_df[,i],fill=strain))+geom_density(alpha=0.6)+theme_classic()+xlab(colnames(umap_df)[i])
  
  p3=ggplot(umap_df, aes(y=umap_df[,i], x=scCoGAPs_1,color=strain))+
    geom_point()+theme_classic()+ylab(colnames(umap_df)[i])+xlab('scCoGAPs_1')
  
  p4=ggplot(umap_df, aes(y=umap_df[,i], x=scCoGAPs_1,color=immune_clusters))+
    geom_point()+theme_classic()+ylab(colnames(umap_df)[i])+xlab('scCoGAPs_1')
  
  p5=ggplot(umap_df, aes(y=umap_df[,i], x=scCoGAPs_1,color=init_clusters))+
    geom_point()+theme_classic()+ylab(colnames(umap_df)[i])+xlab('scCoGAPs_1')
  
  p6=ggplot(umap_df, aes(y=umap_df[,i], x=init_clusters,fill=init_clusters))+geom_violin()+theme_classic()+ylab(colnames(umap_df)[i])
  
  gridExtra::grid.arrange(p1,p2,p6, p3,p4,p5,nrow=2,ncol=3)
  
}
dev.off()
