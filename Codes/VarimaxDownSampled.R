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


old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
new_data_scCLustViz_object <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved.RData"

load(old_data_scClustViz_object)
load(new_data_scCLustViz_object)

merged_samples <- your_scRNAseq_data_object
merged_samples$clusters <- paste(as.character(sCVdata_list$res.0.6@Clusters))
merged_samples$strain <- sapply(strsplit(colnames(merged_samples), '_'), '[[', 2)
merged_samples_orig <- merged_samples
### make freq table for clusters of each map ###
freq.df = data.frame(table(merged_samples$clusters ))
colnames(freq.df) = c('set1-cluster', '#cell')
gridExtra::grid.table(freq.df)


annot_df <- readRDS('cluster_labels/new_labels.rds')
annot_df <- annot_df[,c('cluster','labels')]


##### Finding the UMIs which need to be removed based on the downsampling results ####
#### generating the meta dataframe ####
meta = data.frame(labels=sCVdata_list$res.0.6@Clusters, 
                  umi=names(sCVdata_list$res.0.6@Clusters),
                  strain=sapply(str_split(names(sCVdata_list$res.0.6@Clusters),'_'), '[[', 2))
head(meta)
meta.merged <- merge(meta, annot_df, by.x='labels',by.y='cluster', all.x=T)
colnames(meta.merged) <- c('cluster', 'umi', 'strain', 'labels')
head(meta.merged)


#### down-sampling #### 

#### based on Hep cells in total
meta_split <- split.data.frame(meta.merged, meta.merged$labels)
Hep_df <- meta_split$Hep
Hep_downsamp <- downSample(Hep_df, as.factor(Hep_df$strain), list = FALSE, yname = "Class")
umi_2bRemoved = Hep_df$umi[! Hep_df$umi %in% Hep_downsamp$umi] 

#### based on individual Hep clusters
meta_split_cl <- split.data.frame(meta.merged, meta.merged$cluster)
Hep_clusters <- annot_df$cluster[annot_df$labels=='Hep'] ### including only the Hep clusters
Hep_clusters <- c(Hep_clusters, '3') ## including the Heps + cluster 3 - LSECs in set-2 map
cl_df <- meta_split_cl[Hep_clusters]
cl_downsamp_list <- lapply(cl_df, function(df) downSample(df, as.factor(df$strain), 
                                                     list = FALSE, yname = "Class"))
cl_downsamp <- do.call(rbind, cl_downsamp_list)
umi_2bRemoved= meta.merged$umi[meta.merged$cluster %in% Hep_clusters][! meta.merged$umi[meta.merged$cluster %in% Hep_clusters] %in% cl_downsamp$umi] 

####  down-sampling all the clusters to evaluate effect of #cells
meta_split_cl <- split.data.frame(meta.merged, meta.merged$cluster)
cl_downsamp_list <- lapply(meta_split_cl, function(df) downSample(df, as.factor(df$strain), 
                                                          list = FALSE, yname = "Class"))
cl_downsamp <- do.call(rbind, cl_downsamp_list)
umi_2bRemoved= meta.merged$umi[! meta.merged$umi %in% cl_downsamp$umi] 



############ subsetting the input matrix, re-scaling and PCA ####
merged_samples <- merged_samples[,!colnames(merged_samples) %in% umi_2bRemoved] 
merged_samples <- ScaleData(merged_samples)
merged_samples <- RunPCA(merged_samples, features = rownames(merged_samples), reduction.name = 'pca_sub')  

loading_matrix = Loadings(merged_samples, 'pca_sub')
gene_exp_matrix = GetAssayData(merged_samples, assay = 'RNA')

dim(loading_matrix)
dim(gene_exp_matrix)

### varimax PCA
rot_data <- get_varimax_rotated(gene_exp_matrix, loading_matrix)
rotatedLoadings <- rot_data$rotLoadings
scores <- data.frame(rot_data$rotScores)
colnames(scores) = paste0('Varimax_', 1:ncol(scores))
dim(scores)
dim(rotatedLoadings)
dim(merged_samples) ### running PCA removes around 20 genes which have low variance
head(rotatedLoadings)


rot_standardDev <- apply(rot_data$rotScores, 2, sd)
elbow_plot(rot_standardDev, title = 'Varimax on last 2 rat samples')
ndims = 50
rot_percVar = (rot_standardDev^2 * 100)/sum(rot_standardDev^2)
names(rot_percVar) <- paste0(1:(length(rot_percVar)))
perc_variance_threshold = 0.1
PCs_to_check <- as.numeric(names(rot_percVar[rot_percVar>perc_variance_threshold  ]))

embedd_df_rotated <- data.frame(scores)

##### Visualizing the PCA plots after the varimax rotation #####
pdf('Plots/VarimaxPCA_scTrans_mergedNewSamples_allclusters_DownSampled.pdf',width = 14,height=13) 

for(i in PCs_to_check){ 
  pc_num = i
  rot_df <- data.frame(Varimax_1=embedd_df_rotated$Varimax_1,
                       emb_val=embedd_df_rotated[,pc_num],
                       cluster=as.character(merged_samples$clusters),
                       Library_size=merged_samples$nCount_RNA,
                       num_expressed_genes=merged_samples$nFeature_RNA,
                       strain = merged_samples$strain)
  
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

table(strain=merged_samples$strain, cluster=as.character(merged_samples$clusters))
df <- data.frame(table(strain=merged_samples$strain, cluster=as.character(merged_samples$clusters)))
ggplot(df, aes(fill=strain, x=cluster, y=Freq)) + 
  geom_bar(position="stack", stat="identity")+scale_fill_manual(values=colorPalatte)+ggtitle('all clusters down-sampled')


###############################################################
#################################################
######### Progressive Down-sampling of cluster-13, of set-1, matched with varimax-10 ##############

#### based on individual Hep clusters
meta_split_cl <- split.data.frame(meta.merged, meta.merged$cluster)
a_cluster = '13'
cl_df <- meta_split_cl[a_cluster][[1]]
N = nrow(cl_df)
n_vector = seq(from=10, to=N, by=40)

for (i in 2:length(n_vector)){
  
  ### returning the merged samples back to its original state 
  merged_samples <- merged_samples_orig
  
  #n = n_vector[i]
  n = 120
  sampled_umi_indices = round(runif(n, min=1, max=N))
  sampled_umi = cl_df$umi[sampled_umi_indices]
  
  umi_2bRemoved= cl_df$umi[!cl_df$umi %in% sampled_umi]
  
  
  ############ subsetting the input matrix, re-scaling and PCA ####
  merged_samples <- merged_samples[,!colnames(merged_samples) %in% umi_2bRemoved] 
  merged_samples <- ScaleData(merged_samples)
  merged_samples <- RunPCA(merged_samples, features = rownames(merged_samples), reduction.name = 'pca_sub')  
  
  loading_matrix = Loadings(merged_samples, 'pca_sub')
  gene_exp_matrix = GetAssayData(merged_samples, assay = 'RNA')
  
  dim(loading_matrix)
  dim(gene_exp_matrix)
  
  ### varimax PCA
  rot_data <- get_varimax_rotated(gene_exp_matrix, loading_matrix)
  rotatedLoadings <- rot_data$rotLoadings
  scores <- data.frame(rot_data$rotScores)
  colnames(scores) = paste0('Varimax_', 1:ncol(scores))
  dim(scores)
  dim(rotatedLoadings)
  dim(merged_samples) ### running PCA removes around 20 genes which have low variance
  head(rotatedLoadings)
  
  
  rot_standardDev <- apply(rot_data$rotScores, 2, sd)
  elbow_plot(rot_standardDev, title = 'Varimax on last 2 rat samples')
  ndims = 50
  rot_percVar = (rot_standardDev^2 * 100)/sum(rot_standardDev^2)
  names(rot_percVar) <- paste0(1:(length(rot_percVar)))
  perc_variance_threshold = 0.1
  PCs_to_check <- as.numeric(names(rot_percVar[rot_percVar>perc_variance_threshold  ]))
  
  embedd_df_rotated <- data.frame(scores)
  dev.off()
  ##### Visualizing the PCA plots after the varimax rotation #####
  #pdf('Plots/VarimaxPCA_scTrans_mergedNewSamples_allclusters_DownSampled.pdf',width = 14,height=13) 
  #pdf(paste0('Plots/DownSampling/VarimaxPCA_progressive_DownSampled_size:',n,'.pdf', collapse = ''),width = 14,height=13) 
  pdf('Plots/DownSampling/set1_varimax7_LEW_DownSampled.pdf',width = 14,height=13) 
  pdf('Plots/DownSampling/set1_varimax7_DA_DownSampled.pdf',width = 14,height=13) 
  pdf('Plots/DownSampling/set1_varimax7_balanced.pdf',width = 14,height=13) 
  
  table(strain=merged_samples$strain, cluster=as.character(merged_samples$clusters))
  df <- data.frame(table(strain=merged_samples$strain, cluster=as.character(merged_samples$clusters)))
  ggplot(df, aes(fill=strain, x=cluster, y=Freq)) + 
    geom_bar(position="stack", stat="identity")+scale_fill_manual(values=colorPalatte)+ggtitle('all clusters down-sampled')
  
  for(i in PCs_to_check){ 
    pc_num = i
    rot_df <- data.frame(Varimax_1=embedd_df_rotated$Varimax_1,
                         emb_val=embedd_df_rotated[,pc_num],
                         cluster=as.character(merged_samples$clusters),
                         Library_size=merged_samples$nCount_RNA,
                         num_expressed_genes=merged_samples$nFeature_RNA,
                         strain = merged_samples$strain)
    
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
  
}

##### making varimax-7, representing cluster 7 imbalanced and evaluating if it effects it in any way


varimax_res <- readRDS('Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds')
varimax_df <- data.frame(varimax_res$rotScores)
plot(varimax_df$X1, varimax_df$X7)
varimax_df$var7_cells <- varimax_df$X7 < (-2)
head(varimax_df)

df <- data.frame(umi=rownames(varimax_df), 
                 strain = merged_samples$strain, 
                 cluster = merged_samples$clusters, 
                 is_var7 = varimax_df$var7_cells)

df_var7 <- df[df$is_var7,]
table(df_var7$cluster)
umi_2bRemoved =  df_var7$umi[df_var7$strain=='DA'] # df_var7$umi[df_var7$strain=='DA']


df_var7_ds <- downSample(df_var7, as.factor(df_var7$strain), list = FALSE, yname = "Class")
umi_2bRemoved = df_var7$umi[!df_var7$umi %in% df_var7_ds$umi]
table(df_var7_ds$strain)



